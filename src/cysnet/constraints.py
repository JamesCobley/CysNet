from **future** import annotations

from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

ROUND_DP_DEFAULT = 2
EFF_R_CAP_DEFAULT = 20

def classify_redox_value(value: float) -> str:
"""
Classify a cysteine redox marginal or percent value as reduced, partial or oxidised.

```
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
```

def polytope_dimension_from_intermediates(n_intermediate: int) -> int:
"""
Return the affine dimension of the feasible complete-coverage polytope
after boundary coordinates have collapsed.

```
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
```

def classify_solution_type(complete: bool, n_intermediate: int) -> str:
"""
Classify the observed-coordinate constraint class for a protein/sample.
"""
if not complete:
return "incomplete_observed_coordinate_constraints"

```
if n_intermediate == 0:
    return "exact_singleton"

if n_intermediate == 1:
    return "exact_two_state"

return "inexact_bounded"
```

def resolve_exact_distribution(
marginals,
round_dp: int = ROUND_DP_DEFAULT,
eff_r_cap: int = EFF_R_CAP_DEFAULT,
) -> dict[tuple[int, ...], float] | None:
"""
Resolve an exact observed-coordinate distribution when the feasible
polytope dimension is zero.

```
Returns None if the distribution is bounded or too large to enumerate.

Boundary-only and one-intermediate cases are exact. The more general
rank check is retained for safety.
"""
m = np.asarray(marginals, dtype=float)

fixed = (m == 0.0) | (m == 1.0)
intermediate = (m > 0.0) & (m < 1.0)
j = int(intermediate.sum())

if j == 0:
    return {tuple(int(round(x)) for x in m): 1.0}

if j > eff_r_cap:
    return None

sub = np.array(list(product([0, 1], repeat=j)))

A = np.vstack([sub.T, np.ones((1, sub.shape[0]))])
dim = int(sub.shape[0] - np.linalg.matrix_rank(A))

if dim != 0:
    return None

b = np.concatenate([m[intermediate], [1.0]])
weights, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

intermediate_idx = np.where(intermediate)[0]
dist: dict[tuple[int, ...], float] = {}

for k, substate in enumerate(sub):
    if weights[k] <= 1e-9:
        continue

    vector = np.zeros(len(m), dtype=int)
    vector[fixed] = m[fixed].round().astype(int)

    for ii, coordinate_idx in enumerate(intermediate_idx):
        vector[coordinate_idx] = int(substate[ii])

    fraction = round(float(weights[k]), round_dp)

    if fraction > 0:
        key = tuple(int(x) for x in vector)
        dist[key] = dist.get(key, 0.0) + fraction

return dist
```

def count_oxiform_states(distribution: dict[tuple[int, ...], float]) -> int:
"""
Count substates containing at least one oxidised observed coordinate.
"""
return int(sum(1 for state in distribution if any(bit == 1 for bit in state)))

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

```
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
```

def build_sample_proteoform_summary(
redox_marginals: pd.DataFrame,
per_protein_constraints: pd.DataFrame,
protein_col: str = "protein_id",
sample_col: str = "sample_id",
site_col: str = "site_id",
marginal_col: str = "marginal",
round_dp: int = ROUND_DP_DEFAULT,
eff_r_cap: int = EFF_R_CAP_DEFAULT,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
"""
Build sample-level resolved/constrained proteoform and oxiform summaries.

```
Resolved means complete-coverage exact distributions:
  - exact_singleton
  - exact_two_state

Constrained means:
  - complete but bounded
  - incomplete observed-coordinate constraints

Cohort totals are reported both including repeated sample observations
and as unique resolved proteoforms.

Unique-by-distribution:
  (protein, frozenset((substate, rounded_fraction)))

Unique-by-structure:
  (protein, frozenset(substate_support))
"""
redox = redox_marginals.copy()
constraints = per_protein_constraints.copy()

redox[protein_col] = redox[protein_col].astype(str)
redox[sample_col] = redox[sample_col].astype(str)
redox[site_col] = redox[site_col].astype(str)
redox[marginal_col] = pd.to_numeric(redox[marginal_col], errors="coerce")

redox = redox.sort_values([sample_col, protein_col, site_col]).copy()

constraints[protein_col] = constraints[protein_col].astype(str)
constraints[sample_col] = constraints[sample_col].astype(str)

exact_solution_types = {"exact_singleton", "exact_two_state"}

all_dist_keys = set()
all_struct_keys = set()
all_oxi_dist_keys = set()
all_oxi_struct_keys = set()

total_resolved_multiple = 0
total_resolved_oxi_multiple = 0
total_constrained_multiple = 0
total_constrained_oxi_multiple = 0

per_sample_rows = []
resolved_rows = []

for sample_id, block in constraints.groupby(sample_col, sort=False):
    n_resolved = 0
    n_resolved_oxiform_containing = 0
    n_constrained = 0
    n_constrained_with_oxiform = 0

    redox_sample = redox[redox[sample_col] == sample_id]

    for _, row in block.iterrows():
        protein_id = row[protein_col]
        solution_type = row["solution_type"]

        protein_redox = (
            redox_sample[redox_sample[protein_col] == protein_id]
            .sort_values(site_col)[marginal_col]
            .dropna()
            .to_numpy(dtype=float)
        )

        if len(protein_redox) == 0:
            continue

        if solution_type in exact_solution_types:
            dist = resolve_exact_distribution(
                protein_redox,
                round_dp=round_dp,
                eff_r_cap=eff_r_cap,
            )

            if dist is None:
                n_constrained += 1

                if bool(row["at_least_one_oxiform"]):
                    n_constrained_with_oxiform += 1

                continue

            n_resolved += 1
            n_oxi_states = count_oxiform_states(dist)

            if n_oxi_states > 0:
                n_resolved_oxiform_containing += 1

            distribution_key = (
                protein_id,
                frozenset(dist.items()),
            )

            structure_key = (
                protein_id,
                frozenset(dist.keys()),
            )

            all_dist_keys.add(distribution_key)
            all_struct_keys.add(structure_key)

            if n_oxi_states > 0:
                all_oxi_dist_keys.add(distribution_key)
                all_oxi_struct_keys.add(structure_key)

            resolved_rows.append(
                {
                    sample_col: sample_id,
                    protein_col: protein_id,
                    "solution_type": solution_type,
                    "n_resolved_substates": int(len(dist)),
                    "n_oxiform_substates": int(n_oxi_states),
                    "distribution": ";".join(
                        f"{''.join(map(str, state))}:{fraction}"
                        for state, fraction in sorted(dist.items())
                    ),
                    "support_structure": ";".join(
                        "".join(map(str, state))
                        for state in sorted(dist.keys())
                    ),
                }
            )

        else:
            n_constrained += 1

            if bool(row["at_least_one_oxiform"]):
                n_constrained_with_oxiform += 1

    per_sample_rows.append(
        {
            sample_col: sample_id,
            "resolved_proteoforms": int(n_resolved),
            "resolved_oxiform_containing": int(n_resolved_oxiform_containing),
            "constrained_proteoforms": int(n_constrained),
            "constrained_with_oxiform": int(n_constrained_with_oxiform),
        }
    )

    total_resolved_multiple += n_resolved
    total_resolved_oxi_multiple += n_resolved_oxiform_containing
    total_constrained_multiple += n_constrained
    total_constrained_oxi_multiple += n_constrained_with_oxiform

sample_proteoform_summary = pd.DataFrame(per_sample_rows)

cohort_totals = pd.DataFrame(
    {
        "metric": [
            "resolved_proteoforms_incl_multiples",
            "resolved_proteoforms_unique_by_distribution",
            "resolved_proteoforms_unique_by_structure",
            "resolved_oxiform_incl_multiples",
            "resolved_oxiform_unique_by_distribution",
            "resolved_oxiform_unique_by_structure",
            "constrained_incl_multiples",
            "constrained_with_oxiform_incl_multiples",
        ],
        "value": [
            int(total_resolved_multiple),
            int(len(all_dist_keys)),
            int(len(all_struct_keys)),
            int(total_resolved_oxi_multiple),
            int(len(all_oxi_dist_keys)),
            int(len(all_oxi_struct_keys)),
            int(total_constrained_multiple),
            int(total_constrained_oxi_multiple),
        ],
    }
)

resolved_distributions = pd.DataFrame(resolved_rows)

return sample_proteoform_summary, cohort_totals, resolved_distributions
```

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

```
redox_marginals = pd.read_csv(redox_marginals_path, sep=sep)
protein_topology = pd.read_csv(protein_topology_path, sep=sep)

per_protein_constraints, coverage_classes, constraint_summary = build_constraint_tables(
    redox_marginals=redox_marginals,
    protein_topology=protein_topology,
)

sample_proteoform_summary, cohort_totals, resolved_distributions = (
    build_sample_proteoform_summary(
        redox_marginals=redox_marginals,
        per_protein_constraints=per_protein_constraints,
    )
)

paths = {
    "per_protein_constraints": outdir / f"{study_name}_per_protein_constraints.tsv",
    "coverage_classes": outdir / f"{study_name}_coverage_classes.tsv",
    "constraint_summary": outdir / f"{study_name}_constraint_summary.tsv",
    "sample_proteoform_summary": outdir / f"{study_name}_sample_proteoform_summary.tsv",
    "cohort_proteoform_totals": outdir / f"{study_name}_cohort_proteoform_totals.tsv",
    "resolved_distributions": outdir / f"{study_name}_resolved_distributions.tsv",
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

sample_proteoform_summary.to_csv(
    paths["sample_proteoform_summary"],
    sep="\t",
    index=False,
)

cohort_totals.to_csv(
    paths["cohort_proteoform_totals"],
    sep="\t",
    index=False,
)

resolved_distributions.to_csv(
    paths["resolved_distributions"],
    sep="\t",
    index=False,
)

return paths
```

