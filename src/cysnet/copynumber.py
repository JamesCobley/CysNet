from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


AVOGADRO = 6.02214076e23
DEFAULT_INJECTED_MASS_G = 500e-9

AA_MW = {
    "G": 57.0519,
    "A": 71.0788,
    "S": 87.0782,
    "P": 97.1167,
    "V": 99.1326,
    "T": 101.1051,
    "C": 103.1388,
    "L": 113.1594,
    "I": 113.1594,
    "N": 114.1038,
    "D": 115.0886,
    "Q": 128.1307,
    "K": 128.1741,
    "E": 129.1155,
    "M": 131.1926,
    "H": 137.1411,
    "F": 147.1766,
    "R": 156.1875,
    "Y": 163.1760,
    "W": 186.2132,
    "U": 150.0388,
}

WATER_MW = 18.01528


def parse_fasta_cysteines_and_mw(fasta_path: str | Path) -> dict[str, dict[str, float]]:
    """
    Parse a FASTA file and return accession -> cysteine count and molecular weight.

    Molecular weight is estimated from residue masses plus water.

    UniProt-style headers are handled as:

        >sp|P12345|NAME_HUMAN ...

    where P12345 is used as the accession.

    For non-UniProt headers, the first whitespace-delimited token is used.
    """
    fasta_path = Path(fasta_path)

    records: dict[str, dict[str, float]] = {}
    accession = None
    sequence_parts: list[str] = []

    def flush(acc: str | None, parts: list[str]) -> None:
        if acc is None:
            return

        sequence = "".join(parts).upper()

        if not sequence:
            return

        records[acc] = {
            "fasta_cysteines": int(sequence.count("C")),
            "molecular_weight_da": float(
                sum(AA_MW.get(residue, 110.0) for residue in sequence) + WATER_MW
            ),
        }

    with fasta_path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                flush(accession, sequence_parts)
                sequence_parts = []

                header = line[1:]

                if "|" in header:
                    fields = header.split("|")
                    accession = fields[1] if len(fields) > 1 else fields[0].split()[0]
                else:
                    accession = header.split()[0]

            else:
                sequence_parts.append(line)

        flush(accession, sequence_parts)

    return records


def resolve_protein_group(
    protein_group: str,
    fasta_records: dict[str, dict[str, float]],
) -> str | None:
    """
    Resolve a semicolon-delimited protein group to the first accession found in FASTA.

    Tries exact accession first, then canonical accession with any '-isoform'
    suffix removed.
    """
    for accession in str(protein_group).split(";"):
        accession = accession.strip()

        if accession in fasta_records:
            return accession

        canonical = accession.split("-")[0]
        if canonical in fasta_records:
            return canonical

    return None


def infer_protein_id_column(pg_matrix: pd.DataFrame) -> str:
    """
    Infer the protein-group identifier column from common DIA-NN-style names.
    """
    candidates = [
        "Protein.Group",
        "Protein.Ids",
        "Protein.Group.Ids",
        "Protein",
        "protein_id",
    ]

    for col in candidates:
        if col in pg_matrix.columns:
            return col

    raise ValueError(
        "Could not identify a protein ID column in the protein matrix. "
        f"Tried: {candidates}"
    )


def infer_quantity_columns(
    pg_matrix: pd.DataFrame,
    protein_id_col: str,
) -> list[str]:
    """
    Infer sample/quantity columns in a protein matrix.

    All columns except the protein identifier are treated as possible quantity
    columns. Non-numeric values are coerced later.
    """
    return [c for c in pg_matrix.columns if c != protein_id_col]


def _safe_log2(x: float) -> float:
    if not np.isfinite(x) or x <= 0:
        return np.nan
    return float(np.log2(x))


def copy_constrained_substates(
    detected_cysteines: int,
    protein_copies: float,
) -> tuple[float, bool]:
    """
    Compute the copy-constrained number of realised observed-coordinate substates.

    observed state space = 2^R_detected

    realised substates = min(2^R_detected, protein copies)

    Returns
    -------
    realised_substates, copy_limited
    """
    if detected_cysteines <= 0:
        return np.nan, False

    if not np.isfinite(protein_copies) or protein_copies <= 0:
        return np.nan, False

    log2_copies = np.log2(protein_copies)

    if detected_cysteines <= log2_copies:
        return float(2 ** int(detected_cysteines)), False

    return float(protein_copies), True


def build_copy_number_tables(
    redox_marginals: pd.DataFrame,
    protein_matrix: pd.DataFrame,
    fasta_path: str | Path,
    injected_mass_g: float = DEFAULT_INJECTED_MASS_G,
    protein_col: str = "protein_id",
    sample_col: str = "sample_id",
    site_col: str = "site_id",
    marginal_col: str = "marginal",
    protein_matrix_id_col: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Scale protein-group intensities to molecular copies and calculate
    copy-constrained observed-coordinate substate capacity.

    Parameters
    ----------
    redox_marginals:
        Long CysNet redox marginal table, usually produced by oxidia-sites.

    protein_matrix:
        Wide protein LFQ/intensity matrix. One row per protein group and one
        column per sample/run.

    fasta_path:
        FASTA used for DIA-NN search. Used to estimate molecular weight.

    injected_mass_g:
        Total injected protein mass in grams. Default is 500 ng.

    Returns
    -------
    protein_copy_number:
        One row per sample/protein with scaled mass, copies, detected cysteine
        count and copy-constrained observed-coordinate substate capacity.

    summary:
        One row per sample with copy-number and state-space summaries.
    """
    required = {protein_col, sample_col, site_col, marginal_col}
    missing = required.difference(redox_marginals.columns)

    if missing:
        raise ValueError(f"Redox marginals table is missing required columns: {sorted(missing)}")

    if injected_mass_g <= 0:
        raise ValueError("Injected mass must be positive.")

    fasta_records = parse_fasta_cysteines_and_mw(fasta_path)

    if protein_matrix_id_col is None:
        protein_matrix_id_col = infer_protein_id_column(protein_matrix)

    sample_ids = sorted(redox_marginals[sample_col].dropna().astype(str).unique())

    available_quantity_cols = infer_quantity_columns(
        protein_matrix,
        protein_id_col=protein_matrix_id_col,
    )

    sample_columns = [s for s in sample_ids if s in available_quantity_cols]

    if not sample_columns:
        raise ValueError(
            "No sample columns in the protein matrix matched sample_id values in "
            "the redox marginals table."
        )

    pg = protein_matrix.copy()
    pg[protein_matrix_id_col] = pg[protein_matrix_id_col].astype(str)
    pg = pg.set_index(protein_matrix_id_col)

    detected = (
        redox_marginals
        .dropna(subset=[protein_col, sample_col, site_col])
        .groupby([sample_col, protein_col], as_index=False)[site_col]
        .nunique()
        .rename(columns={site_col: "detected_cysteines"})
    )

    protein_rows = []

    for sample_id in sample_columns:
        intensity = pd.to_numeric(pg[sample_id], errors="coerce")
        total_intensity = float(intensity.sum(skipna=True))

        if not np.isfinite(total_intensity) or total_intensity <= 0:
            continue

        sample_detected = detected[detected[sample_col].astype(str) == str(sample_id)]

        for _, row in sample_detected.iterrows():
            protein_id = str(row[protein_col])
            detected_cys = int(row["detected_cysteines"])

            if protein_id not in pg.index:
                continue

            resolved = resolve_protein_group(protein_id, fasta_records)

            if resolved is None:
                continue

            mw = fasta_records[resolved]["molecular_weight_da"]
            fasta_cys = int(fasta_records[resolved]["fasta_cysteines"])

            raw_intensity = pd.to_numeric(pd.Series([pg.loc[protein_id, sample_id]]), errors="coerce").iloc[0]

            if not np.isfinite(raw_intensity) or raw_intensity <= 0:
                continue

            scaled_mass_g = raw_intensity / total_intensity * injected_mass_g
            protein_copies = scaled_mass_g / mw * AVOGADRO

            realised, copy_limited = copy_constrained_substates(
                detected_cysteines=detected_cys,
                protein_copies=protein_copies,
            )

            protein_rows.append(
                {
                    "sample_id": sample_id,
                    "protein_id": protein_id,
                    "resolved_accession": resolved,
                    "raw_intensity": float(raw_intensity),
                    "total_sample_intensity": total_intensity,
                    "scaled_mass_g": float(scaled_mass_g),
                    "molecular_weight_da": float(mw),
                    "protein_copies": float(protein_copies),
                    "fasta_cysteines": fasta_cys,
                    "detected_cysteines": detected_cys,
                    "coverage_percent": (
                        detected_cys / fasta_cys * 100.0
                        if fasta_cys > 0
                        else np.nan
                    ),
                    "complete": detected_cys == fasta_cys,
                    "log10_full_state_space": fasta_cys * np.log10(2.0),
                    "log10_observed_state_space": detected_cys * np.log10(2.0),
                    "realised_substates": realised,
                    "copy_limited": bool(copy_limited),
                }
            )

    protein_copy_number = pd.DataFrame(protein_rows)

    if protein_copy_number.empty:
        raise ValueError(
            "No protein copy-number rows could be generated. Check that protein IDs, "
            "sample IDs and FASTA accessions match."
        )

    redox_for_weighting = redox_marginals.copy()
    redox_for_weighting[sample_col] = redox_for_weighting[sample_col].astype(str)
    redox_for_weighting[protein_col] = redox_for_weighting[protein_col].astype(str)

    copy_lookup = protein_copy_number[
        ["sample_id", "protein_id", "protein_copies"]
    ].copy()

    weighted = redox_for_weighting.merge(
        copy_lookup,
        left_on=[sample_col, protein_col],
        right_on=["sample_id", "protein_id"],
        how="inner",
    )

    if "percent_oxidised" in weighted.columns:
        weighted["oxidation_percent_for_weighting"] = pd.to_numeric(
            weighted["percent_oxidised"],
            errors="coerce",
        )
    else:
        weighted["oxidation_percent_for_weighting"] = (
            pd.to_numeric(weighted[marginal_col], errors="coerce") * 100.0
        )

    summary_rows = []

    for sample_id, block in protein_copy_number.groupby("sample_id", sort=False):
        n_valid = int(len(block))
        sum_realised = float(block["realised_substates"].sum(skipna=True))
        floor = n_valid

        w_block = weighted[weighted["sample_id"] == sample_id]
        mask = (
            np.isfinite(w_block["oxidation_percent_for_weighting"])
            & np.isfinite(w_block["protein_copies"])
            & (w_block["protein_copies"] > 0)
        )

        if mask.any():
            copy_weighted_mean_ox_pct = float(
                np.sum(
                    w_block.loc[mask, "oxidation_percent_for_weighting"]
                    * w_block.loc[mask, "protein_copies"]
                )
                / np.sum(w_block.loc[mask, "protein_copies"])
            )
        else:
            copy_weighted_mean_ox_pct = np.nan

        summary_rows.append(
            {
                "sample_id": sample_id,
                "n_cys_proteins": n_valid,
                "total_scaled_protein_copies": float(block["protein_copies"].sum(skipna=True)),
                "sum_realised_substates": sum_realised,
                "floor_one_per_protein": floor,
                "realised_above_one_floor": sum_realised - floor,
                "pct_copy_limited": (
                    float(block["copy_limited"].mean() * 100.0)
                    if n_valid
                    else np.nan
                ),
                "copy_weighted_mean_ox_pct": copy_weighted_mean_ox_pct,
            }
        )

    summary = pd.DataFrame(summary_rows)

    return protein_copy_number, summary


def write_copy_number_outputs(
    redox_marginals_path: str | Path,
    protein_matrix_path: str | Path,
    fasta_path: str | Path,
    outdir: str | Path,
    study_name: str,
    injected_mass_g: float = DEFAULT_INJECTED_MASS_G,
    sep: str = "\t",
    protein_matrix_sep: str = "\t",
) -> dict[str, Path]:
    """
    Write copy-number scaling outputs from CysNet redox marginals, a protein
    LFQ/intensity matrix and FASTA.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    redox_marginals = pd.read_csv(redox_marginals_path, sep=sep)
    protein_matrix = pd.read_csv(protein_matrix_path, sep=protein_matrix_sep)

    protein_copy_number, summary = build_copy_number_tables(
        redox_marginals=redox_marginals,
        protein_matrix=protein_matrix,
        fasta_path=fasta_path,
        injected_mass_g=injected_mass_g,
    )

    paths = {
        "protein_copy_number": outdir / f"{study_name}_protein_copy_number.tsv",
        "copy_substate_summary": outdir / f"{study_name}_copy_substate_summary.tsv",
    }

    protein_copy_number.to_csv(paths["protein_copy_number"], sep="\t", index=False)
    summary.to_csv(paths["copy_substate_summary"], sep="\t", index=False)

    return paths
