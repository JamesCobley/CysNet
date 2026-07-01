from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

def frechet_oxiform_bounds(marginals) -> tuple[float, float]:
"""
Bounds for the probability of at least one oxidised observed cysteine.

```
For events A_i = cysteine i oxidised:

    lower P(union A_i) = max_i P(A_i)
    upper P(union A_i) = min(1, sum_i P(A_i))

These are sharp Fréchet bounds without assuming independence.
"""
marg = np.asarray(marginals, dtype=float)

if len(marg) == 0:
    return np.nan, np.nan

return float(np.max(marg)), float(min(1.0, np.sum(marg)))
```

def fully_reduced_bounds(marginals) -> tuple[float, float]:
"""
Bounds for the probability that all observed cysteines are reduced.

```
This is the complement of the union of oxidised events:

    P(all reduced) = 1 - P(any oxidised)

Therefore:

    lower = max(0, 1 - sum(p_i))
    upper = 1 - max(p_i)
"""
marg = np.asarray(marginals, dtype=float)

if len(marg) == 0:
    return np.nan, np.nan

return float(max(0.0, 1.0 - np.sum(marg))), float(1.0 - np.max(marg))
```

def observed_coordinate_string(marginals) -> str:
"""
Represent observed marginals as a compact coordinate string.

```
Boundary values are 0/1; partial values are X.
"""
chars = []

for p in np.asarray(marginals, dtype=float):
    if np.isclose(p, 0.0):
        chars.append("0")
    elif np.isclose(p, 1.0):
        chars.append("1")
    else:
        chars.append("X")

return "".join(chars)
```

def exact_distribution_if_dim0(marginals) -> list[dict[str, float]]:
"""
Return exact observed-coordinate distribution for dimension-zero cases.

```
Cases:
  - no intermediate marginals: one deterministic observed-coordinate state
  - one intermediate marginal: exact two-state distribution

For >=2 intermediate marginals, the feasible set is generally bounded,
so this returns an empty list.
"""
marg = np.asarray(marginals, dtype=float)

inter_idx = np.where((marg > 0.0) & (marg < 1.0))[0]

if len(inter_idx) == 0:
    state = "".join("1" if np.isclose(p, 1.0) else "0" for p in marg)
    return [{"state": state, "probability": 1.0}]

if len(inter_idx) == 1:
    idx = int(inter_idx[0])
    p = float(marg[idx])

    state0 = []
    state1 = []

    for i, q in enumerate(marg):
        if i == idx:
            state0.append("0")
            state1.append("1")
        elif np.isclose(q, 1.0):
            state0.append("1")
            state1.append("1")
        else:
            state0.append("0")
            state1.append("0")

    return [
        {"state": "".join(state0), "probability": 1.0 - p},
        {"state": "".join(state1), "probability": p},
    ]

return []
```

def build_copy_constraint_tables(
redox_marginals: pd.DataFrame,
per_protein_constraints: pd.DataFrame,
protein_copy_number: pd.DataFrame,
protein_col: str = "protein_id",
sample_col: str = "sample_id",
site_col: str = "site_id",
marginal_col: str = "marginal",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
"""
Build copy-number-scaled CysNet constraint tables.

```
Inputs
------
redox_marginals:
    Long-format CysNet redox marginal table.

per_protein_constraints:
    Output from `cysnet.constraints.write_constraint_outputs`.

protein_copy_number:
    Output from `cysnet.copynumber.write_copy_number_outputs`.

Outputs
-------
copy_constraints:
    One row per protein/sample with copy-scaled oxiform and fully reduced
    bounds.

exact_substate_copies:
    One row per exact observed-coordinate substate for exact singleton and
    exact two-state cases.

copy_constraint_summary:
    One row per sample.

copy_constraint_cohort_summary:
    One-row cohort summary.
"""
required_redox = {protein_col, sample_col, site_col, marginal_col}
missing_redox = required_redox.difference(redox_marginals.columns)

if missing_redox:
    raise ValueError(
        "Redox marginals table is missing required columns: "
        f"{sorted(missing_redox)}"
    )

required_constraints = {
    protein_col,
    sample_col,
    "R_total",
    "R_detected",
    "coverage",
    "n_degenerate",
    "n_intermediate",
    "solution_type",
    "at_least_one_oxiform",
}
missing_constraints = required_constraints.difference(per_protein_constraints.columns)

if missing_constraints:
    raise ValueError(
        "Per-protein constraints table is missing required columns: "
        f"{sorted(missing_constraints)}"
    )

required_copy = {protein_col, sample_col, "protein_copies"}
missing_copy = required_copy.difference(protein_copy_number.columns)

if missing_copy:
    raise ValueError(
        "Protein copy-number table is missing required columns: "
        f"{sorted(missing_copy)}"
    )

redox = redox_marginals.copy()
constraints = per_protein_constraints.copy()
copies = protein_copy_number.copy()

for df in [redox, constraints, copies]:
    df[protein_col] = df[protein_col].astype(str)
    df[sample_col] = df[sample_col].astype(str)

redox[site_col] = redox[site_col].astype(str)
redox[marginal_col] = pd.to_numeric(redox[marginal_col], errors="coerce")
redox = redox.dropna(subset=[protein_col, sample_col, site_col, marginal_col])
redox = redox.sort_values([sample_col, protein_col, site_col]).copy()

copies["protein_copies"] = pd.to_numeric(copies["protein_copies"], errors="coerce")

copy_lookup = copies[
    [sample_col, protein_col, "protein_copies"]
].drop_duplicates(
    subset=[sample_col, protein_col],
    keep="first",
)

merged = constraints.merge(
    copy_lookup,
    on=[sample_col, protein_col],
    how="left",
)

redox_grouped = {
    key: block.sort_values(site_col)[marginal_col].to_numpy(dtype=float)
    for key, block in redox.groupby([sample_col, protein_col], sort=False)
}

constraint_rows = []
exact_rows = []

for _, rec in merged.iterrows():
    sample_id = str(rec[sample_col])
    protein_id = str(rec[protein_col])

    marg = redox_grouped.get((sample_id, protein_id), np.array([], dtype=float))
    cp = rec["protein_copies"]

    has_copies = np.isfinite(cp) and cp > 0

    Rdet = int(rec["R_detected"])
    Rtot = int(rec["R_total"])
    k_fixed = int(rec["n_degenerate"])
    n_inter = int(rec["n_intermediate"])

    ox_min_frac, ox_max_frac = frechet_oxiform_bounds(marg)
    red_min_frac, red_max_frac = fully_reduced_bounds(marg)

    ox_min_copies = ox_min_frac * cp if has_copies and np.isfinite(ox_min_frac) else np.nan
    ox_max_copies = ox_max_frac * cp if has_copies and np.isfinite(ox_max_frac) else np.nan

    red_min_copies = red_min_frac * cp if has_copies and np.isfinite(red_min_frac) else np.nan
    red_max_copies = red_max_frac * cp if has_copies and np.isfinite(red_max_frac) else np.nan

    compatible_fibre_copies = cp if has_copies else np.nan

    observed_space_log2 = Rdet
    pruned_space_log2 = Rdet - k_fixed
    excluded_fraction = 1.0 - (2.0 ** (-k_fixed)) if k_fixed >= 0 else np.nan

    min_required_observed_substates = 2 if n_inter > 0 else 1

    copy_limited_realised_observed_substates = np.nan

    if has_copies:
        if Rdet <= np.log2(cp):
            copy_limited_realised_observed_substates = float(2**Rdet)
        else:
            copy_limited_realised_observed_substates = float(cp)

    constraint_rows.append(
        {
            sample_col: sample_id,
            protein_col: protein_id,
            "coverage": rec["coverage"],
            "solution_type": rec["solution_type"],
            "R_total": Rtot,
            "R_detected": Rdet,
            "observed_coordinate": observed_coordinate_string(marg),
            "protein_copy_number": float(cp) if has_copies else np.nan,
            "observed_space_log2": observed_space_log2,
            "pruned_space_log2": pruned_space_log2,
            "collapse_extent_log2": k_fixed,
            "excluded_fraction_by_boundary_priors": excluded_fraction,
            "copy_limited_realised_observed_substates": copy_limited_realised_observed_substates,
            "min_required_observed_substates": min_required_observed_substates,
            "at_least_one_oxiform": bool(rec["at_least_one_oxiform"]),
            "oxiform_min_fraction": ox_min_frac,
            "oxiform_max_fraction": ox_max_frac,
            "oxiform_min_copies": ox_min_copies,
            "oxiform_max_copies": ox_max_copies,
            "fully_reduced_min_fraction": red_min_frac,
            "fully_reduced_max_fraction": red_max_frac,
            "fully_reduced_min_copies": red_min_copies,
            "fully_reduced_max_copies": red_max_copies,
            "compatible_fibre_copies": compatible_fibre_copies,
        }
    )

    if (
        rec["coverage"] == "complete"
        and rec["solution_type"] in {"exact_singleton", "exact_two_state"}
        and has_copies
    ):
        dist = exact_distribution_if_dim0(marg)

        for state_rec in dist:
            state = state_rec["state"]
            probability = float(state_rec["probability"])

            exact_rows.append(
                {
                    sample_col: sample_id,
                    protein_col: protein_id,
                    "state": state,
                    "probability": probability,
                    "substate_copies": probability * cp,
                    "is_oxiform": "1" in state,
                    "k_oxidised_observed": state.count("1"),
                    "protein_copy_number": float(cp),
                }
            )

copy_constraints = pd.DataFrame(constraint_rows)
exact_substate_copies = pd.DataFrame(exact_rows)

summary_rows = []

for sample_id, block in copy_constraints.groupby(sample_col, sort=False):
    summary_rows.append(
        {
            sample_col: sample_id,
            "protein_records": int(len(block)),
            "records_with_copy_number": int(block["protein_copy_number"].notna().sum()),
            "total_protein_copies_in_cysnet_records": float(
                block["protein_copy_number"].sum(skipna=True)
            ),
            "records_with_oxiform_evidence": int(block["at_least_one_oxiform"].sum()),
            "oxiform_min_copies_total": float(block["oxiform_min_copies"].sum(skipna=True)),
            "oxiform_max_copies_total": float(block["oxiform_max_copies"].sum(skipna=True)),
            "fully_reduced_min_copies_total": float(
                block["fully_reduced_min_copies"].sum(skipna=True)
            ),
            "fully_reduced_max_copies_total": float(
                block["fully_reduced_max_copies"].sum(skipna=True)
            ),
            "compatible_fibre_copies_total": float(
                block["compatible_fibre_copies"].sum(skipna=True)
            ),
            "mean_oxiform_min_fraction": float(
                block["oxiform_min_fraction"].mean(skipna=True)
            ),
            "mean_oxiform_max_fraction": float(
                block["oxiform_max_fraction"].mean(skipna=True)
            ),
            "total_collapse_extent_log2": float(
                block["collapse_extent_log2"].sum(skipna=True)
            ),
        }
    )

copy_constraint_summary = pd.DataFrame(summary_rows)

if len(exact_substate_copies):
    exact_oxiform_substate_rows = int(exact_substate_copies["is_oxiform"].sum())
    exact_substate_copies_total = float(
        exact_substate_copies["substate_copies"].sum(skipna=True)
    )
    exact_oxiform_substate_copies_total = float(
        exact_substate_copies.loc[
            exact_substate_copies["is_oxiform"],
            "substate_copies",
        ].sum(skipna=True)
    )
else:
    exact_oxiform_substate_rows = 0
    exact_substate_copies_total = 0.0
    exact_oxiform_substate_copies_total = 0.0

copy_constraint_cohort_summary = pd.DataFrame(
    [
        {
            "protein_sample_records": int(len(copy_constraints)),
            "records_with_copy_number": int(
                copy_constraints["protein_copy_number"].notna().sum()
            ),
            "total_protein_copies_in_cysnet_records": float(
                copy_constraints["protein_copy_number"].sum(skipna=True)
            ),
            "records_with_oxiform_evidence": int(
                copy_constraints["at_least_one_oxiform"].sum()
            ),
            "oxiform_min_copies_total": float(
                copy_constraints["oxiform_min_copies"].sum(skipna=True)
            ),
            "oxiform_max_copies_total": float(
                copy_constraints["oxiform_max_copies"].sum(skipna=True)
            ),
            "fully_reduced_min_copies_total": float(
                copy_constraints["fully_reduced_min_copies"].sum(skipna=True)
            ),
            "fully_reduced_max_copies_total": float(
                copy_constraints["fully_reduced_max_copies"].sum(skipna=True)
            ),
            "exact_substate_rows": int(len(exact_substate_copies)),
            "exact_oxiform_substate_rows": exact_oxiform_substate_rows,
            "exact_substate_copies_total": exact_substate_copies_total,
            "exact_oxiform_substate_copies_total": exact_oxiform_substate_copies_total,
        }
    ]
)

return (
    copy_constraints,
    exact_substate_copies,
    copy_constraint_summary,
    copy_constraint_cohort_summary,
)
```

def write_copy_constraint_outputs(
redox_marginals_path: str | Path,
per_protein_constraints_path: str | Path,
protein_copy_number_path: str | Path,
outdir: str | Path,
study_name: str,
sep: str = "\t",
) -> dict[str, Path]:
"""
Write copy-number-scaled CysNet constraint outputs.
"""
outdir = Path(outdir)
outdir.mkdir(parents=True, exist_ok=True)

```
redox_marginals = pd.read_csv(redox_marginals_path, sep=sep)
per_protein_constraints = pd.read_csv(per_protein_constraints_path, sep=sep)
protein_copy_number = pd.read_csv(protein_copy_number_path, sep=sep)

(
    copy_constraints,
    exact_substate_copies,
    copy_constraint_summary,
    copy_constraint_cohort_summary,
) = build_copy_constraint_tables(
    redox_marginals=redox_marginals,
    per_protein_constraints=per_protein_constraints,
    protein_copy_number=protein_copy_number,
)

paths = {
    "copy_constraints": outdir / f"{study_name}_copy_constraints.tsv",
    "exact_substate_copies": outdir / f"{study_name}_exact_substate_copies.tsv",
    "copy_constraint_summary": outdir / f"{study_name}_copy_constraint_summary.tsv",
    "copy_constraint_cohort_summary": outdir / f"{study_name}_copy_constraint_cohort_summary.tsv",
}

copy_constraints.to_csv(paths["copy_constraints"], sep="\t", index=False)
exact_substate_copies.to_csv(paths["exact_substate_copies"], sep="\t", index=False)
copy_constraint_summary.to_csv(paths["copy_constraint_summary"], sep="\t", index=False)
copy_constraint_cohort_summary.to_csv(
    paths["copy_constraint_cohort_summary"],
    sep="\t",
    index=False,
)

return paths
```
