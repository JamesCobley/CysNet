from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def parse_fasta_cysteine_counts(fasta_path: str | Path) -> dict[str, int]:
    """
    Parse a FASTA file and return accession -> cysteine count.

    UniProt-style headers are handled as:

        >sp|P12345|NAME_HUMAN ...

    where P12345 is used as the accession.

    For non-UniProt headers, the first whitespace-delimited token is used.
    """
    fasta_path = Path(fasta_path)

    cys_counts: dict[str, int] = {}
    accession = None
    sequence_parts: list[str] = []

    def flush(acc: str | None, parts: list[str]) -> None:
        if acc is not None:
            sequence = "".join(parts).upper()
            cys_counts[acc] = sequence.count("C")

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

    return cys_counts


def resolve_protein_group_to_fasta(
    protein_group: str,
    cysteine_counts: dict[str, int],
) -> tuple[str | None, int | float]:
    """
    Resolve a semicolon-delimited protein group to the first accession found in FASTA.

    The function first tries the accession exactly as written. If that fails,
    it tries stripping an isoform suffix after '-'.

    Returns
    -------
    resolved_accession, fasta_cysteine_count
    """
    for accession in str(protein_group).split(";"):
        accession = accession.strip()

        if accession in cysteine_counts:
            return accession, cysteine_counts[accession]

        canonical = accession.split("-")[0]
        if canonical in cysteine_counts:
            return canonical, cysteine_counts[canonical]

    return None, np.nan


def _log10_sum_power2(r_values: Iterable[float]) -> float:
    """
    Compute log10(sum(2^R)) without requiring scipy.

    This is used for theoretical proteoform-space summaries.
    """
    r = np.asarray(list(r_values), dtype=float)
    r = r[np.isfinite(r)]

    if len(r) == 0:
        return np.nan

    x = r * np.log(2.0)
    xmax = np.max(x)

    return float((xmax + np.log(np.sum(np.exp(x - xmax)))) / np.log(10.0))


def build_protein_topology_table(
    redox_marginals: pd.DataFrame,
    fasta_path: str | Path,
    protein_col: str = "protein_id",
    sample_col: str = "sample_id",
    site_col: str = "site_id",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build per-protein/sample cysteine topology bookkeeping from CysNet redox marginals.

    Parameters
    ----------
    redox_marginals:
        Long CysNet-ready redox table. Must contain protein_id, sample_id and site_id.

    fasta_path:
        FASTA file used to define total cysteine topology.

    Returns
    -------
    protein_topology:
        One row per sample/protein with FASTA cysteine count, detected cysteine
        count, coverage percentage and completeness.

    summary:
        One row per sample with coverage and theoretical state-space summaries.
    """
    required = {protein_col, sample_col, site_col}
    missing = required.difference(redox_marginals.columns)

    if missing:
        raise ValueError(f"Redox marginals table is missing required columns: {sorted(missing)}")

    cysteine_counts = parse_fasta_cysteine_counts(fasta_path)

    detected = (
        redox_marginals
        .dropna(subset=[protein_col, sample_col, site_col])
        .groupby([sample_col, protein_col], as_index=False)[site_col]
        .nunique()
        .rename(columns={site_col: "detected_cysteines"})
    )

    protein_groups = detected[protein_col].drop_duplicates()

    resolved_rows = []
    for protein_group in protein_groups:
        resolved_accession, fasta_cysteines = resolve_protein_group_to_fasta(
            protein_group,
            cysteine_counts,
        )

        resolved_rows.append(
            {
                protein_col: protein_group,
                "resolved_accession": resolved_accession,
                "fasta_cysteines": fasta_cysteines,
            }
        )

    resolved = pd.DataFrame(resolved_rows)

    protein_topology = detected.merge(resolved, on=protein_col, how="left")

    protein_topology["matched_fasta"] = protein_topology["resolved_accession"].notna()

    protein_topology = protein_topology[protein_topology["matched_fasta"]].copy()

    protein_topology["fasta_cysteines"] = protein_topology["fasta_cysteines"].astype(int)

    protein_topology["coverage_percent"] = (
        protein_topology["detected_cysteines"]
        / protein_topology["fasta_cysteines"]
        * 100.0
    )

    protein_topology["complete"] = (
        protein_topology["detected_cysteines"] == protein_topology["fasta_cysteines"]
    )

    protein_topology["log10_full_state_space"] = (
        protein_topology["fasta_cysteines"] * np.log10(2.0)
    )

    protein_topology["log10_observed_state_space"] = (
        protein_topology["detected_cysteines"] * np.log10(2.0)
    )

    protein_topology = protein_topology.rename(
        columns={
            protein_col: "protein_id",
            sample_col: "sample_id",
        }
    )

    summary_rows = []

    for sample_id, block in protein_topology.groupby("sample_id", sort=False):
        total_fasta_cys = int(block["fasta_cysteines"].sum())
        total_detected_cys = int(block["detected_cysteines"].sum())

        summary_rows.append(
            {
                "sample_id": sample_id,
                "n_protein_groups_detected": int(len(block)),
                "n_complete_protein_groups": int(block["complete"].sum()),
                "n_incomplete_protein_groups": int((~block["complete"]).sum()),
                "n_fasta_cysteines": total_fasta_cys,
                "n_detected_cysteines": total_detected_cys,
                "cysteine_coverage_percent": (
                    total_detected_cys / total_fasta_cys * 100.0
                    if total_fasta_cys > 0
                    else np.nan
                ),
                "log10_sum_full_state_space": _log10_sum_power2(block["fasta_cysteines"]),
                "log10_sum_observed_state_space": _log10_sum_power2(block["detected_cysteines"]),
            }
        )

    summary = pd.DataFrame(summary_rows)

    return protein_topology, summary


def write_topology_outputs(
    redox_marginals_path: str | Path,
    fasta_path: str | Path,
    outdir: str | Path,
    study_name: str,
    sep: str = "\t",
) -> dict[str, Path]:
    """
    Write CysNet FASTA topology outputs from a redox_marginals table.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    redox_marginals = pd.read_csv(redox_marginals_path, sep=sep)

    protein_topology, summary = build_protein_topology_table(
        redox_marginals=redox_marginals,
        fasta_path=fasta_path,
    )

    paths = {
        "protein_topology": outdir / f"{study_name}_protein_topology.tsv",
        "topology_summary": outdir / f"{study_name}_topology_summary.tsv",
    }

    protein_topology.to_csv(paths["protein_topology"], sep="\t", index=False)
    summary.to_csv(paths["topology_summary"], sep="\t", index=False)

    return paths
