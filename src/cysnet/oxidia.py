from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


DEFAULT_ID_COLS = ["Protein", "Protein.Names", "Gene.Names", "Residue", "Site", "Sequence"]
DEFAULT_KEY_COLS = ["Protein", "Residue", "Site", "Sequence"]
DEFAULT_META_COLS = ["Protein.Names", "Gene.Names"]


def read_table(path: str | Path, sep: str = "\t") -> pd.DataFrame:
    """
    Read a tabular Oxi-DIA input file.

    DIA-NN site output is usually tab-delimited, but this function is kept
    simple so that comma-delimited files can also be supported later.
    """
    return pd.read_csv(path, sep=sep)


def infer_measurement_columns(
    df: pd.DataFrame,
    id_cols: Iterable[str] = DEFAULT_ID_COLS,
) -> list[str]:
    """
    Return columns that are treated as quantitative run/sample columns.

    All columns not listed as site identity or metadata columns are treated
    as measurement columns.
    """
    id_set = set(id_cols)
    return [c for c in df.columns if c not in id_set]


def validate_oxidia_site_tables(
    light: pd.DataFrame,
    heavy: pd.DataFrame,
    key_cols: Iterable[str] = DEFAULT_KEY_COLS,
) -> None:
    """
    Validate reduced/light and oxidised/heavy site tables before redox calculation.
    """
    key_cols = list(key_cols)

    for name, df in {"light": light, "heavy": heavy}.items():
        missing = [c for c in key_cols if c not in df.columns]
        if missing:
            raise ValueError(f"{name} table is missing key columns: {missing}")

        duplicated = df.duplicated(key_cols).sum()
        if duplicated:
            raise ValueError(
                f"{name} table contains {duplicated} duplicated site keys. "
                "Each Protein/Residue/Site/Sequence combination must be unique."
            )


def compute_site_redox_from_light_heavy(
    light: pd.DataFrame,
    heavy: pd.DataFrame,
    sample_columns: Optional[list[str]] = None,
    key_cols: Iterable[str] = DEFAULT_KEY_COLS,
    meta_cols: Iterable[str] = DEFAULT_META_COLS,
    id_cols: Iterable[str] = DEFAULT_ID_COLS,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute site-level cysteine redox from reduced/light and oxidised/heavy matrices.

    Redox logic
    -----------
    Only light and heavy site intensities define the redox marginal:

        oxidised fraction = H / (H + L)
        percent oxidised = H / (H + L) * 100

    Single-channel sites are retained:

        L only -> 0%
        H only -> 100%
        L + H  -> H / (H + L)
        neither -> missing

    Returns
    -------
    site_percent:
        Wide table with site metadata and percent oxidised per sample/run.

    site_coverage:
        Wide table with site metadata and a 0/1 detection flag per sample/run.

    redox_long:
        Long CysNet-ready table with fractional marginals.
    """
    key_cols = list(key_cols)
    meta_cols = list(meta_cols)

    validate_oxidia_site_tables(light, heavy, key_cols=key_cols)

    if sample_columns is None:
        light_samples = set(infer_measurement_columns(light, id_cols=id_cols))
        heavy_samples = set(infer_measurement_columns(heavy, id_cols=id_cols))
        sample_columns = sorted(light_samples.union(heavy_samples))

    light_i = light.set_index(key_cols)
    heavy_i = heavy.set_index(key_cols)

    keys = light_i.index.union(heavy_i.index)

    light_i = light_i.reindex(keys)
    heavy_i = heavy_i.reindex(keys)

    meta = pd.DataFrame(index=keys)

    for c in meta_cols:
        if c in heavy_i.columns and c in light_i.columns:
            meta[c] = heavy_i[c].fillna(light_i[c])
        elif c in heavy_i.columns:
            meta[c] = heavy_i[c]
        elif c in light_i.columns:
            meta[c] = light_i[c]

    percent = {}
    detected = {}

    for c in sample_columns:
        if c in heavy_i.columns:
            h = pd.to_numeric(heavy_i[c], errors="coerce")
        else:
            h = pd.Series(np.nan, index=keys)

        if c in light_i.columns:
            l = pd.to_numeric(light_i[c], errors="coerce")
        else:
            l = pd.Series(np.nan, index=keys)

        if (h.dropna() < 0).any() or (l.dropna() < 0).any():
            raise ValueError(f"Negative intensity detected in sample column: {c}")

        total = h.fillna(0.0) + l.fillna(0.0)
        is_detected = total > 0

        percent[c] = np.where(is_detected, h.fillna(0.0) / total * 100.0, np.nan)
        detected[c] = is_detected.astype(int)

    percent_df = pd.DataFrame(percent, index=keys)
    detected_df = pd.DataFrame(detected, index=keys)

    site_percent = meta.join(percent_df).reset_index()
    site_coverage = meta.join(detected_df.add_suffix("_detected")).reset_index()

    long_rows = []

    for c in sample_columns:
        tmp = percent_df[[c]].rename(columns={c: "percent_oxidised"}).reset_index()
        tmp["sample_id"] = c
        tmp["marginal"] = tmp["percent_oxidised"] / 100.0
        tmp["detected"] = detected_df[c].to_numpy().astype(bool)

        tmp = tmp[tmp["detected"]].copy()

        tmp = tmp.rename(
            columns={
                "Protein": "protein_id",
                "Site": "site_id",
                "Residue": "residue",
                "Sequence": "sequence",
            }
        )

        long_rows.append(
            tmp[
                [
                    "protein_id",
                    "sample_id",
                    "site_id",
                    "residue",
                    "sequence",
                    "marginal",
                    "percent_oxidised",
                    "detected",
                ]
            ]
        )

    redox_long = pd.concat(long_rows, ignore_index=True)

    return site_percent, site_coverage, redox_long


def summarise_site_redox(
    site_percent: pd.DataFrame,
    site_coverage: pd.DataFrame,
    key_cols: Iterable[str] = DEFAULT_KEY_COLS,
    meta_cols: Iterable[str] = DEFAULT_META_COLS,
) -> pd.DataFrame:
    """
    Summarise site-level redox per sample/run.

    Reports mean oxidation, median oxidation, detected sites and degenerate
    0%/100% sites.
    """
    non_sample_cols = set(key_cols).union(meta_cols)
    sample_cols = [c for c in site_percent.columns if c not in non_sample_cols]

    rows = []

    for c in sample_cols:
        values = pd.to_numeric(site_percent[c], errors="coerce")
        coverage_col = f"{c}_detected"

        if coverage_col in site_coverage.columns:
            detected = pd.to_numeric(site_coverage[coverage_col], errors="coerce") > 0
        else:
            detected = values.notna()

        degenerate = ((values == 0.0) | (values == 100.0)) & detected

        rows.append(
            {
                "sample_id": c,
                "global_mean_percent_ox": float(values.mean(skipna=True)),
                "median_percent_ox": float(values.median(skipna=True)),
                "n_sites_detected": int(detected.sum()),
                "n_sites_degenerate_0_or_100": int(degenerate.sum()),
            }
        )

    return pd.DataFrame(rows)


def write_oxidia_outputs(
    light_path: str | Path,
    heavy_path: str | Path,
    outdir: str | Path,
    study_name: str,
    sep: str = "\t",
) -> dict[str, Path]:
    """
    Load L/H Oxi-DIA site tables, calculate redox marginals and write outputs.

    This is the first CysNet import layer for Oxi-DIA style inputs.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    light = read_table(light_path, sep=sep)
    heavy = read_table(heavy_path, sep=sep)

    site_percent, site_coverage, redox_long = compute_site_redox_from_light_heavy(
        light=light,
        heavy=heavy,
    )

    summary = summarise_site_redox(site_percent, site_coverage)

    paths = {
        "site_percent": outdir / f"{study_name}_site_percent_oxidised.tsv",
        "site_coverage": outdir / f"{study_name}_site_coverage_nfiles.tsv",
        "summary": outdir / f"{study_name}_sample_summary.tsv",
        "redox_marginals": outdir / f"{study_name}_redox_marginals.tsv",
    }

    site_percent.to_csv(paths["site_percent"], sep="\t", index=False)
    site_coverage.to_csv(paths["site_coverage"], sep="\t", index=False)
    summary.to_csv(paths["summary"], sep="\t", index=False)
    redox_long.to_csv(paths["redox_marginals"], sep="\t", index=False)

    return paths
