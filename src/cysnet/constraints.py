from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def classify_redox_value(value: float) -> str:
    """
    Classify a cysteine redox marginal or percent value as reduced, partial or oxidised.

    Accepts either fractional values in [0, 1] or percent values in [0, 100].
    """
    if pd.isna(value):
        return "missing"

    value = float(value)

    if value > 1.0:
        value = value / 100.0

    if value == 0.0:
        return "reduced(0)"

    if value == 1.0:
        return "oxidised(1)"

    if 0.0 < value < 1.0:
        return "partial"

    raise ValueError(f"Redox value outside expected range: {value}")


def polytope_dimension_from_intermediates(n_intermediate: int) -> int:
    """
    Return the affine dimension of the feasible complete-coverage polytope
    after boundary coordinates have collapsed.

    For j intermediate coordinates:

        j = 0 -> exact singleton
        j = 1 -> exact two-state distribution, dimension 0
        j >= 2 -> 2^j - j - 1

    This is the dimension of the distribution over the non-degenerate
    observed-coordinate state space under normalisation plus j marginal
    constraints.
    """
    j = int(n_intermediate)

    if j < 0:
        raise ValueError("n_intermediate must be non-negative.")

    if j <= 1:
        return 0

    return int((2**j) - j - 1)


def classify_solution_type(complete: bool, n_intermediate: int) -> str:
    """
    Classify the observed-coordinate constraint class for a protein/sample.
    """
    if not complete:
        return "incomplete_observed_coordinate_constraints"

    if n_intermediate == 0:
        return "exact_singleton"

    if n_intermediate == 1:
        return "exact_two_state"

    return "inexact_bounded"


def build_constraint_tables(
    redox_marginals: pd.DataFrame,
    protein_topology: pd.DataFrame,
    protein_col: str = "protein_id",
    sample_col: str = "sample_id",
    site_col: str = "site_id",
    marginal_col: str = "marginal",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build CysNet coverage and theorem-constraint summary tables.

    Parameters
    ----------
    redox_marginals:
        Long-format CysNet redox marginal table, usually produced by
        `cysnet oxidia-sites`.

        Required columns by default:
            protein_id, sample_id, site_id, marginal

    protein_topology:
        Protein/sample topology table, usually produced by `cysnet topology`.

        Required columns by default:
            protein_id, sample_id, fasta_cysteines, detected_cysteines, complete

    Returns
    -------
    per_protein_constraints:
        One row per protein/sample with complete/incomplete coverage,
        redox-coordinate composition and theorem constraint class.

    coverage_classes:
        One row per sample summarising redox composition stratified by
        complete versus incomplete protein coverage.

    constraint_summary:
        One row per sample summarising exact, bounded and incomplete
        observed-coordinate constraint classes.
    """
    required_redox = {protein_col, sample_col, site_col, marginal_col}
    missing_redox = required_redox.difference(redox_marginals.columns)

    if missing_redox:
        raise ValueError(
            "Redox marginals table is missing required columns: "
            f"{sorted(missing_redox)}"
        )

    required_topology = {
        protein_col,
        sample_col,
        "fasta_cysteines",
        "detected_cysteines",
        "complete",
    }
    missing_topology = required_topology.difference(protein_topology.columns)

    if missing_topology:
        raise ValueError(
            "Protein topology table is missing required columns: "
            f"{sorted(missing_topology)}"
        )

    redox = redox_marginals.copy()
    topology = protein_topology.copy()

    redox[protein_col] = redox[protein_col].astype(str)
    redox[sample_col] = redox[sample_col].astype(str)
    redox[site_col] = redox[site_col].astype(str)
    redox[marginal_col] = pd.to_numeric(redox[marginal_col], errors="coerce")

    topology[protein_col] = topology[protein_col].astype(str)
    topology[sample_col] = topology[sample_col].astype(str)
    topology["fasta_cysteines"] = pd.to_numeric(
        topology["fasta_cysteines"],
        errors="coerce",
    )
    topology["detected_cysteines"] = pd.to_numeric(
        topology["detected_cysteines"],
        errors="coerce",
    )

    if topology["complete"].dtype != bool:
        topology["complete"] = topology["complete"].astype(str).str.lower().isin(
            ["true", "1", "yes"]
        )

    # ------------------------------------------------------------------
    # Site-level redox composition, then attach coverage class.
    # ------------------------------------------------------------------
    site = redox.dropna(
        subset=[protein_col, sample_col, site_col, marginal_col]
    ).copy()

    site["redox_class"] = site[marginal_col].map(classify_redox_value)

    site_with_topology = site.merge(
        topology[
            [
                sample_col,
                protein_col,
                "fasta_cysteines",
                "detected_cysteines",
                "complete",
            ]
        ],
        on=[sample_col, protein_col],
        how="left",
    )

    site_with_topology = site_with_topology.dropna(
        subset=["fasta_cysteines", "detected_cysteines"]
    ).copy()

    site_with_topology["coverage"] = np.where(
        site_with_topology["complete"],
        "complete",
        "incomplete",
    )

    # ------------------------------------------------------------------
    # Per-protein constraint table.
    # ------------------------------------------------------------------
    grouped = site.groupby([sample_col, protein_col], sort=False)

    redox_counts = grouped[marginal_col].agg(
        n_detected_sites="count",
        n_fixed_reduced=lambda x: int((x == 0.0).sum()),
        n_fixed_oxidised=lambda x: int((x == 1.0).sum()),
        n_intermediate=lambda x: int(((x > 0.0) & (x < 1.0)).sum()),
        at_least_one_oxiform=lambda x: bool((x > 0.0).any()),
    ).reset_index()

    per_protein = topology.merge(
        redox_counts,
        on=[sample_col, protein_col],
        how="inner",
    )

    per_protein = per_protein.rename(
        columns={
            "fasta_cysteines": "R_total",
            "detected_cysteines": "R_detected",
        }
    )

    per_protein["R_total"] = per_protein["R_total"].astype(int)
    per_protein["R_detected"] = per_protein["R_detected"].astype(int)
    per_protein["n_detected_sites"] = per_protein["n_detected_sites"].astype(int)
    per_protein["n_fixed_reduced"] = per_protein["n_fixed_reduced"].astype(int)
    per_protein["n_fixed_oxidised"] = per_protein["n_fixed_oxidised"].astype(int)
    per_protein["n_intermediate"] = per_protein["n_intermediate"].astype(int)

    per_protein["coverage"] = np.where(
        per_protein["complete"],
        "complete",
        "incomplete",
    )

    per_protein["n_degenerate"] = (
        per_protein["n_fixed_reduced"] + per_protein["n_fixed_oxidised"]
    )

    per_protein["observed_state_space_log2"] = per_protein["R_detected"]
    per_protein["full_state_space_log2"] = per_protein["R_total"]

    # Boundary coordinates collapse one binary dimension each in the observed
    # state space. Unobserved cysteines are not inferred.
    per_protein["collapsed_space_log2"] = (
        per_protein["R_detected"] - per_protein["n_degenerate"]
    )
    per_protein["collapse_extent_log2"] = per_protein["n_degenerate"]

    per_protein["multi_intermediate"] = per_protein["n_intermediate"] >= 2

    per_protein["polytope_dim"] = per_protein["n_intermediate"].map(
        polytope_dimension_from_intermediates
    )

    per_protein.loc[
        ~per_protein["complete"],
        "polytope_dim",
    ] = np.nan

    per_protein["solution_type"] = [
        classify_solution_type(bool(complete), int(n_intermediate))
        for complete, n_intermediate in zip(
            per_protein["complete"],
            per_protein["n_intermediate"],
        )
    ]

    keep_cols = [
        sample_col,
        protein_col,
        "R_total",
        "R_detected",
        "coverage",
        "n_degenerate",
        "n_fixed_reduced",
        "n_fixed_oxidised",
        "n_intermediate",
        "observed_state_space_log2",
        "full_state_space_log2",
        "collapsed_space_log2",
        "collapse_extent_log2",
        "at_least_one_oxiform",
        "multi_intermediate",
        "polytope_dim",
        "solution_type",
    ]

    optional_cols = [
        c
        for c in [
            "resolved_accession",
            "coverage_percent",
            "log10_full_state_space",
            "log10_observed_state_space",
        ]
        if c in per_protein.columns
    ]

    per_protein_constraints = per_protein[keep_cols + optional_cols].copy()

    # ------------------------------------------------------------------
    # Coverage class redox composition.
    # ------------------------------------------------------------------
    coverage_records = []

    samples = sorted(site_with_topology[sample_col].dropna().astype(str).unique())

    for sample_id in samples:
        sample_sites = site_with_topology[
            site_with_topology[sample_col].astype(str) == str(sample_id)
        ]

        ctab = pd.crosstab(
            sample_sites["coverage"],
            sample_sites["redox_class"],
        ).reindex(
            index=["complete", "incomplete"],
            columns=["reduced(0)", "partial", "oxidised(1)"],
            fill_value=0,
        )

        sample_proteins = per_protein_constraints[
            per_protein_constraints[sample_col].astype(str) == str(sample_id)
        ]

        complete = sample_proteins[sample_proteins["coverage"] == "complete"]
        incomplete = sample_proteins[sample_proteins["coverage"] == "incomplete"]

        complete_polytope_ge2partial = int(
            (
                (complete["n_intermediate"] >= 2)
                & (complete["solution_type"] == "inexact_bounded")
            ).sum()
        )

        coverage_records.append(
            {
                sample_col: sample_id,
                "proteins_complete": int(len(complete)),
                "proteins_incomplete": int(len(incomplete)),
                "cys_complete": int(ctab.loc["complete"].sum()),
                "cys_incomplete": int(ctab.loc["incomplete"].sum()),
                "complete_reduced_0": int(ctab.loc["complete", "reduced(0)"]),
                "complete_partial": int(ctab.loc["complete", "partial"]),
                "complete_oxidised_1": int(ctab.loc["complete", "oxidised(1)"]),
                "incomplete_reduced_0": int(ctab.loc["incomplete", "reduced(0)"]),
                "incomplete_partial": int(ctab.loc["incomplete", "partial"]),
                "incomplete_oxidised_1": int(ctab.loc["incomplete", "oxidised(1)"]),
                "complete_polytope_ge2partial": complete_polytope_ge2partial,
            }
        )

    coverage_classes = pd.DataFrame(coverage_records)

    # ------------------------------------------------------------------
    # Constraint summary.
    # ------------------------------------------------------------------
    summary_rows = []

    for sample_id, block in per_protein_constraints.groupby(sample_col, sort=False):
        complete = block[block["coverage"] == "complete"]
        incomplete = block[block["coverage"] == "incomplete"]

        summary_rows.append(
            {
                sample_col: sample_id,
                "proteins_total": int(len(block)),
                "proteins_complete": int(len(complete)),
                "proteins_incomplete": int(len(incomplete)),
                "cys_complete": int(complete["R_detected"].sum()),
                "cys_incomplete": int(incomplete["R_detected"].sum()),
                "exact_singletons": int(
                    (complete["solution_type"] == "exact_singleton").sum()
                ),
                "exact_two_state": int(
                    (complete["solution_type"] == "exact_two_state").sum()
                ),
                "inexact_bounded": int(
                    (complete["solution_type"] == "inexact_bounded").sum()
                ),
                "incomplete_with_oxiform": int(
                    incomplete["at_least_one_oxiform"].sum()
                ),
                "complete_with_oxiform": int(
                    complete["at_least_one_oxiform"].sum()
                ),
                "proteins_with_oxiform_total": int(
                    block["at_least_one_oxiform"].sum()
                ),
                "proteins_multi_intermediate": int(
                    block["multi_intermediate"].sum()
                ),
                "total_collapse_extent_log2": int(
                    block["collapse_extent_log2"].sum()
                ),
            }
        )

    constraint_summary = pd.DataFrame(summary_rows)

    return per_protein_constraints, coverage_classes, constraint_summary


def write_constraint_outputs(
    redox_marginals_path: str | Path,
    protein_topology_path: str | Path,
    outdir: str | Path,
    study_name: str,
    sep: str = "\t",
) -> dict[str, Path]:
    """
    Write CysNet coverage and theorem-constraint outputs.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    redox_marginals = pd.read_csv(redox_marginals_path, sep=sep)
    protein_topology = pd.read_csv(protein_topology_path, sep=sep)

    per_protein_constraints, coverage_classes, constraint_summary = build_constraint_tables(
        redox_marginals=redox_marginals,
        protein_topology=protein_topology,
    )

    paths = {
        "per_protein_constraints": outdir / f"{study_name}_per_protein_constraints.tsv",
        "coverage_classes": outdir / f"{study_name}_coverage_classes.tsv",
        "constraint_summary": outdir / f"{study_name}_constraint_summary.tsv",
    }

    per_protein_constraints.to_csv(
        paths["per_protein_constraints"],
        sep="\t",
        index=False,
    )

    coverage_classes.to_csv(
        paths["coverage_classes"],
        sep="\t",
        index=False,
    )

    constraint_summary.to_csv(
        paths["constraint_summary"],
        sep="\t",
        index=False,
    )

    return paths
