from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd


def _read_if_exists(path: Path) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return None


def _section(title: str) -> str:
    line = "=" * len(title)
    return f"\n{title}\n{line}\n"


def _table_preview(df: pd.DataFrame | None, max_rows: int = 12) -> str:
    if df is None:
        return "Not available.\n"

    if df.empty:
        return "Available, but empty.\n"

    return df.head(max_rows).to_string(index=False) + "\n"


def _metric_table(df: pd.DataFrame | None) -> str:
    if df is None:
        return "Not available.\n"

    if df.empty:
        return "Available, but empty.\n"

    return df.to_string(index=False) + "\n"


def _safe_value(df: pd.DataFrame | None, column: str, default="not available"):
    if df is None or df.empty or column not in df.columns:
        return default

    value = df[column].iloc[0]

    if pd.isna(value):
        return default

    return value


def write_meaning_report(
    outdir: str | Path,
    study_name: str,
) -> Path:
    """
    Write a human-readable CysNet meaning report.

    The report automatically detects which outputs are present.

    Required/core outputs:
      - sample summary
      - topology summary
      - constraint summary
      - coverage classes
      - sample proteoform summary
      - cohort proteoform totals

    Optional copy-number outputs:
      - protein copy-number summary
      - copy constraints
      - exact substate copies
      - copy constraint summaries
    """
    outdir = Path(outdir)

    paths = {
        "site_percent": outdir / f"{study_name}_site_percent_oxidised.tsv",
        "site_coverage": outdir / f"{study_name}_site_coverage_nfiles.tsv",
        "sample_summary": outdir / f"{study_name}_sample_summary.tsv",
        "redox_marginals": outdir / f"{study_name}_redox_marginals.tsv",
        "protein_topology": outdir / f"{study_name}_protein_topology.tsv",
        "topology_summary": outdir / f"{study_name}_topology_summary.tsv",
        "per_protein_constraints": outdir / f"{study_name}_per_protein_constraints.tsv",
        "coverage_classes": outdir / f"{study_name}_coverage_classes.tsv",
        "constraint_summary": outdir / f"{study_name}_constraint_summary.tsv",
        "sample_proteoform_summary": outdir / f"{study_name}_sample_proteoform_summary.tsv",
        "cohort_proteoform_totals": outdir / f"{study_name}_cohort_proteoform_totals.tsv",
        "resolved_distributions": outdir / f"{study_name}_resolved_distributions.tsv",
        "protein_copy_number": outdir / f"{study_name}_protein_copy_number.tsv",
        "copy_substate_summary": outdir / f"{study_name}_copy_substate_summary.tsv",
        "copy_constraints": outdir / f"{study_name}_copy_constraints.tsv",
        "exact_substate_copies": outdir / f"{study_name}_exact_substate_copies.tsv",
        "copy_constraint_summary": outdir / f"{study_name}_copy_constraint_summary.tsv",
        "copy_constraint_cohort_summary": outdir / f"{study_name}_copy_constraint_cohort_summary.tsv",
    }

    data = {name: _read_if_exists(path) for name, path in paths.items()}

    has_copy_number = data["protein_copy_number"] is not None
    has_copy_constraints = data["copy_constraints"] is not None

    report = []

    report.append("CysNet meaning report\n")
    report.append("=====================\n\n")
    report.append(f"Study: {study_name}\n")
    report.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n")
    report.append(f"Output folder: {outdir}\n")

    report.append(_section("1. Workflow mode"))

    if has_copy_number:
        report.append(
            "CysNet detected protein copy-number outputs. This run used the copy-aware workflow.\n\n"
        )
        report.append(
            "Mode: redox marginals + FASTA topology + theorem constraints + protein copy-number scaling.\n"
        )
    else:
        report.append(
            "CysNet did not detect protein copy-number outputs. This run used the fraction-scale workflow.\n\n"
        )
        report.append(
            "Mode: redox marginals + FASTA topology + theorem constraints only.\n"
        )

    report.append(
        "\nThe redox calculation is always based only on the L and H cysteine-site channels:\n"
    )
    report.append("  oxidised fraction = H / (H + L)\n")
    report.append("  L-only sites are 0% oxidised.\n")
    report.append("  H-only sites are 100% oxidised.\n")
    report.append("  sites with neither channel are missing.\n")

    report.append(_section("2. Files detected"))

    for name, path in paths.items():
        status = "present" if path.exists() else "missing"
        report.append(f"{name}: {status} :: {path}\n")

    report.append(_section("3. Oxi-DIA redox summary"))
    report.append(
        "This table summarises the site-level redox import from the L/H matrices.\n\n"
    )
    report.append(_table_preview(data["sample_summary"]))

    report.append(_section("4. FASTA topology summary"))
    report.append(
        "This table maps detected cysteine sites onto the FASTA-derived cysteine topology.\n"
    )
    report.append(
        "Complete proteins have all FASTA cysteines detected. Incomplete proteins are interpreted only over the observed cysteine-coordinate projection.\n\n"
    )
    report.append(_table_preview(data["topology_summary"]))

    report.append(_section("5. Constraint summary"))
    report.append(
        "This table summarises exact, bounded and incomplete observed-coordinate constraint classes.\n\n"
    )
    report.append(_table_preview(data["constraint_summary"]))

    report.append("\nInterpretation of solution_type:\n")
    report.append("  exact_singleton: complete coverage and one exact observed-coordinate state.\n")
    report.append("  exact_two_state: complete coverage and one intermediate cysteine, giving an exact two-state distribution.\n")
    report.append("  inexact_bounded: complete coverage but multiple feasible distributions remain.\n")
    report.append("  incomplete_observed_coordinate_constraints: incomplete coverage; CysNet constrains only measured coordinates.\n")

    report.append(_section("6. Coverage classes"))
    report.append(
        "This table separates detected cysteine redox classes by complete versus incomplete protein coverage.\n\n"
    )
    report.append(_table_preview(data["coverage_classes"]))

    report.append(_section("7. Sample-level proteoform and oxiform summary"))
    report.append(
        "Resolved proteoforms are complete-coverage exact observed-coordinate distributions.\n"
    )
    report.append(
        "Constrained proteoforms are complete bounded cases plus incomplete observed-coordinate constraints.\n\n"
    )
    report.append(_table_preview(data["sample_proteoform_summary"]))

    report.append(_section("8. Cohort proteoform totals"))
    report.append(
        "Unique-by-distribution counts treat the same protein with the same state probabilities as the same resolved proteoform distribution.\n"
    )
    report.append(
        "Unique-by-structure counts ignore probabilities and compare only the occupied substate support.\n\n"
    )
    report.append(_metric_table(data["cohort_proteoform_totals"]))

    report.append(_section("9. Resolved exact distributions"))
    resolved = data["resolved_distributions"]

    if resolved is None:
        report.append("No resolved distribution table was found.\n")
    elif resolved.empty:
        report.append("No exact resolved distributions were recorded.\n")
    else:
        report.append(
            f"Resolved exact distribution rows: {len(resolved):,}\n\n"
        )
        report.append(_table_preview(resolved))

    report.append(_section("10. Copy-number scaling"))

    if not has_copy_number:
        report.append(
            "No protein copy-number table was detected. Copy-number-scaled interpretation was not run.\n"
        )
    else:
        report.append(
            "Protein-group intensities were scaled to the injected protein mass and converted to molecular copies using FASTA-derived molecular weight.\n\n"
        )
        report.append("Copy substate summary:\n\n")
        report.append(_table_preview(data["copy_substate_summary"]))

    report.append(_section("11. Copy-number-scaled constraints"))

    if not has_copy_constraints:
        report.append(
            "No copy-constraint table was detected. This is expected if no PG / LFQ matrix was supplied.\n"
        )
    else:
        report.append(
            "CysNet scaled observed-coordinate oxiform and fully reduced bounds into molecular copy-number bounds.\n\n"
        )
        report.append("Per-sample copy constraint summary:\n\n")
        report.append(_table_preview(data["copy_constraint_summary"]))

        report.append("\nCohort copy constraint summary:\n\n")
        report.append(_metric_table(data["copy_constraint_cohort_summary"]))

        exact_copies = data["exact_substate_copies"]

        if exact_copies is not None and not exact_copies.empty:
            report.append("\nExact substate copy-number rows:\n\n")
            report.append(_table_preview(exact_copies))
        else:
            report.append("\nNo exact substate copy-number rows were recorded.\n")

    report.append(_section("12. Mathematical interpretation"))

    report.append(
        "CysNet treats cysteine oxidation values as first-order marginals over binary observed-coordinate substates.\n"
    )
    report.append(
        "For complete-coverage proteins, CysNet may resolve an exact observed-coordinate oxiform distribution.\n"
    )
    report.append(
        "For incomplete proteins, CysNet does not infer unmeasured cysteine states. It reports exact constraints over the measured coordinate projection.\n\n"
    )

    report.append("For a binary state s and cysteine marginals m_j:\n")
    report.append("  q_j(s) = m_j if s_j = 1, else 1 - m_j\n")
    report.append("  lower_s = max(0, sum_j q_j(s) - (R - 1))\n")
    report.append("  upper_s = min_j q_j(s)\n\n")

    report.append(
        "For copy-aware runs, the observed-coordinate oxiform fraction is bounded by sharp union bounds:\n"
    )
    report.append("  oxiform_min_fraction = max(p_i)\n")
    report.append("  oxiform_max_fraction = min(1, sum(p_i))\n\n")
    report.append("The fully reduced observed-coordinate fraction is bounded by:\n")
    report.append("  fully_reduced_min_fraction = max(0, 1 - sum(p_i))\n")
    report.append("  fully_reduced_max_fraction = 1 - max(p_i)\n\n")

    report.append(
        "All copy-number-scaled quantities are bounds over the observed cysteine-coordinate projection.\n"
    )

    report.append(_section("13. Reviewer-safe interpretation"))

    report.append(
        "CysNet does not claim to directly identify full molecular proteoforms from bottom-up data.\n"
    )
    report.append(
        "Instead, it reports theorem-constrained oxiform/proteoform subspaces compatible with measured cysteine redox marginals.\n"
    )
    report.append(
        "Complete coverage can yield exact observed-coordinate distributions. Incomplete coverage yields exact observed-coordinate constraints, not full unmeasured proteoform resolution.\n"
    )

    report_path = outdir / f"{study_name}_CysNet_meaning_report.txt"
    report_path.write_text("".join(report), encoding="utf-8")

    return report_path
