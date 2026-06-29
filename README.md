# CysNet

CysNet is a theorem-constrained workflow for inferring oxiform constraints from bottom-up cysteine redox proteomics.

CysNet treats site-level cysteine oxidation values as marginals over an unobserved distribution of cysteine-redox proteoform substates. Given site-level redox measurements, FASTA-derived cysteine topology and optional protein copy-number estimates, CysNet reports compatible, excluded, required, exact and bounded oxiform substates.

CysNet does not over-resolve incomplete bottom-up data. It only returns what the mathematics supports. Complete-coverage proteins may collapse to exact oxiform ensemble distributions; incomplete proteins return the exact constraint class supported by the measured cysteine coordinates.

## What CysNet does

CysNet converts bottom-up cysteine redox proteomics into a protein-level constraint report.

The core workflow is:

```text
L / UniMod_108 reduced site matrix
H / UniMod_776 oxidised site matrix
FASTA used for the DIA-NN search
optional protein-group LFQ / PG matrix
        ↓
site-level redox marginals
        ↓
FASTA-derived cysteine topology
        ↓
complete / incomplete coverage classes
        ↓
theorem-constrained proteoform and oxiform summaries
        ↓
optional copy-number-scaled oxiform bounds
        ↓
downloadable CysNet output bundle
```

## Current status

CysNet currently implements:

* theorem-constrained state-bound inference from cysteine redox marginals;
* Oxi-DIA light/heavy site matrix import;
* redox marginal calculation from reduced and oxidised cysteine channels;
* FASTA-derived cysteine topology bookkeeping;
* complete versus incomplete coverage classification;
* sample-level proteoform and oxiform summaries;
* cohort-level resolved and constrained proteoform totals;
* optional protein copy-number scaling from protein-group intensities;
* optional copy-number-scaled oxiform and fully reduced bounds;
* exact substate copy numbers for exact observed-coordinate solutions;
* command-line workflows for theorem testing, Oxi-DIA import and topology;
* a Colab-native upload helper for non-command-line use;
* optional experimental widget and Streamlit routes;
* unit tests for the core theorem.

## Redox logic

For Oxi-DIA inputs, CysNet assumes:

```text
L = Light / UniMod_108 = reduced cysteine signal
H = Heavy / UniMod_776 = reversibly oxidised cysteine signal
```

The cysteine redox marginal is calculated only from the L and H site channels:

```text
oxidised fraction = H / (H + L)
percent oxidised = H / (H + L) * 100
```

Single-channel sites are retained:

```text
L only      -> 0% oxidised
H only      -> 100% oxidised
L + H       -> H / (H + L)
neither L/H -> missing
```

The protein LFQ / PG matrix is not used to calculate redox marginals. Protein abundance or copy-number information is used only after redox calculation for copy-number-scaled oxiform constraints.

## Installation

CysNet can be installed directly from the GitHub repository into a local Python environment or a Google Colab notebook.

### Clone the repository

```bash
git clone https://github.com/JamesCobley/CysNet.git
cd CysNet
```

### Install CysNet with development dependencies

```bash
python -m pip install -e ".[dev]" --upgrade --upgrade-strategy only-if-needed
```

This installs CysNet together with the dependencies required to run the tests.

### Optional extras

For the experimental notebook widget:

```bash
python -m pip install -e ".[dev,widget]" --upgrade --upgrade-strategy only-if-needed
```

For the optional Streamlit app:

```bash
python -m pip install -e ".[dev,app]" --upgrade --upgrade-strategy only-if-needed
```

## Run the tests

```bash
pytest -q
```

Expected output:

```text
7 passed
```

These tests check the core theorem behaviour, including boundary exclusion, exact oxiform resolution, bounded feasible ensembles, oxiform union bounds, fully reduced bounds and invalid input handling.

## Recommended Colab workflow

The recommended non-command-line route is the Colab-native upload helper.

This avoids third-party widget issues and uses Colab's built-in file upload/download system.

In a fresh Colab notebook:

```python
%cd /content
!rm -rf CysNet
!git clone https://github.com/JamesCobley/CysNet.git
%cd /content/CysNet

!python -m pip install -e ".[dev]" --upgrade --upgrade-strategy only-if-needed
!pytest -q
```

Then run:

```python
import sys
sys.path.insert(0, "/content/CysNet/src")

from cysnet.colab import run_colab_upload

run_colab_upload()
```

The Colab helper asks for:

```text
study name
L/H delimiter
PG delimiter
injected protein mass in ng
L / Light / UniMod_108 site matrix
H / Heavy / UniMod_776 site matrix
FASTA file used for the DIA-NN search
optional PG / protein LFQ matrix
```

If no PG matrix is provided, CysNet runs the fraction-scale workflow:

```text
Oxi-DIA redox import
FASTA topology
constraint classification
sample-level proteoform summary
cohort proteoform totals
resolved exact distributions
ZIP download
```

If a PG matrix is provided, CysNet additionally runs:

```text
protein copy-number scaling
copy-constrained observed-coordinate substate capacity
copy-number-scaled oxiform bounds
copy-number-scaled fully reduced bounds
exact substate copy numbers where the solution is exact
copy-aware summary tables
```

## Main outputs

A typical CysNet run writes:

```text
<study>_site_percent_oxidised.tsv
<study>_site_coverage_nfiles.tsv
<study>_sample_summary.tsv
<study>_redox_marginals.tsv
<study>_protein_topology.tsv
<study>_topology_summary.tsv
<study>_per_protein_constraints.tsv
<study>_coverage_classes.tsv
<study>_constraint_summary.tsv
<study>_sample_proteoform_summary.tsv
<study>_cohort_proteoform_totals.tsv
<study>_resolved_distributions.tsv
```

If a PG / protein LFQ matrix is supplied, CysNet also writes:

```text
<study>_protein_copy_number.tsv
<study>_copy_substate_summary.tsv
<study>_copy_constraints.tsv
<study>_exact_substate_copies.tsv
<study>_copy_constraint_summary.tsv
<study>_copy_constraint_cohort_summary.tsv
```

## Command-line usage

CysNet currently exposes command-line workflows for theorem checks, Oxi-DIA import and FASTA topology.

```text
cysnet theorem
cysnet oxidia-sites
cysnet topology
```

The constraint and copy-constraint modules can currently be run through the Python API or through the Colab helper.

## 1. Theorem examples

The theorem command takes cysteine oxidation marginals as fractions.

### Exact observed-coordinate solution

```bash
cysnet theorem 0 0 0.25
```

Expected result:

```text
solution_type    exact_singleton
```

This corresponds to three observed cysteine coordinates with marginals:

```text
C1 = 0
C2 = 0
C3 = 0.25
```

The only compatible distribution is:

```text
000 = 0.75
001 = 0.25
```

All states containing oxidation at C1 or C2 are excluded.

### Bounded observed-coordinate solution

```bash
cysnet theorem 0 0.25 0.25
```

Expected result:

```text
solution_type    bounded
```

This corresponds to three observed cysteine coordinates with marginals:

```text
C1 = 0
C2 = 0.25
C3 = 0.25
```

CysNet excludes all states containing oxidation at C1, but the remaining observed-coordinate oxiform ensemble is not uniquely resolved. The feasible distribution is bounded rather than exact.

## 2. Oxi-DIA site import

The Oxi-DIA importer takes two site-level matrices:

```text
L / Light / UniMod_108 site matrix
H / Heavy / UniMod_776 site matrix
```

Run:

```bash
cysnet oxidia-sites \
  --light UniMod_108_sites.tsv \
  --heavy UniMod_776_sites.tsv \
  --study MY_STUDY \
  --out results
```

This writes:

```text
results/MY_STUDY_site_percent_oxidised.tsv
results/MY_STUDY_site_coverage_nfiles.tsv
results/MY_STUDY_sample_summary.tsv
results/MY_STUDY_redox_marginals.tsv
```

The key CysNet-ready output is:

```text
MY_STUDY_redox_marginals.tsv
```

This is a long-format table containing fractional redox marginals for each detected protein, sample and cysteine site.

## 3. FASTA topology bookkeeping

The topology command maps detected cysteine sites onto the FASTA used for the DIA-NN search.

Run:

```bash
cysnet topology \
  --redox-marginals results/MY_STUDY_redox_marginals.tsv \
  --fasta search_database.fasta \
  --study MY_STUDY \
  --out results
```

This writes:

```text
results/MY_STUDY_protein_topology.tsv
results/MY_STUDY_topology_summary.tsv
```

The protein topology table reports, for each protein and sample:

```text
sample_id
protein_id
resolved_accession
fasta_cysteines
detected_cysteines
coverage_percent
complete
log10_full_state_space
log10_observed_state_space
```

The `complete` column indicates whether all FASTA cysteines were detected for that protein/sample.

```text
complete = True  -> detected cysteines equal FASTA cysteine count
complete = False -> observed-coordinate constraint only
```

## 4. Constraint classification through Python

After generating redox marginals and topology outputs, run:

```python
from cysnet.constraints import write_constraint_outputs

constraint_paths = write_constraint_outputs(
    redox_marginals_path="results/MY_STUDY_redox_marginals.tsv",
    protein_topology_path="results/MY_STUDY_protein_topology.tsv",
    outdir="results",
    study_name="MY_STUDY",
)
```

This writes:

```text
results/MY_STUDY_per_protein_constraints.tsv
results/MY_STUDY_coverage_classes.tsv
results/MY_STUDY_constraint_summary.tsv
results/MY_STUDY_sample_proteoform_summary.tsv
results/MY_STUDY_cohort_proteoform_totals.tsv
results/MY_STUDY_resolved_distributions.tsv
```

### Per-protein constraints

`per_protein_constraints.tsv` reports, for each protein/sample:

```text
sample_id
protein_id
R_total
R_detected
coverage
n_degenerate
n_fixed_reduced
n_fixed_oxidised
n_intermediate
observed_state_space_log2
full_state_space_log2
collapsed_space_log2
collapse_extent_log2
at_least_one_oxiform
multi_intermediate
polytope_dim
solution_type
```

The `solution_type` column can be:

```text
exact_singleton
exact_two_state
inexact_bounded
incomplete_observed_coordinate_constraints
```

### Coverage classes

`coverage_classes.tsv` summarises cysteine redox classes within complete and incomplete proteins:

```text
complete_reduced_0
complete_partial
complete_oxidised_1
incomplete_reduced_0
incomplete_partial
incomplete_oxidised_1
complete_polytope_ge2partial
```

### Proteoform and oxiform summaries

`sample_proteoform_summary.tsv` reports:

```text
resolved_proteoforms
resolved_oxiform_containing
constrained_proteoforms
constrained_with_oxiform
```

`cohort_proteoform_totals.tsv` reports resolved and constrained totals, including:

```text
resolved_proteoforms_incl_multiples
resolved_proteoforms_unique_by_distribution
resolved_proteoforms_unique_by_structure
resolved_oxiform_incl_multiples
resolved_oxiform_unique_by_distribution
resolved_oxiform_unique_by_structure
constrained_incl_multiples
constrained_with_oxiform_incl_multiples
```

`resolved_distributions.tsv` records exact observed-coordinate distributions for exact complete-coverage solutions.

## 5. Protein copy-number scaling through Python

If a PG / protein LFQ matrix is available, CysNet can scale protein-group intensities to molecular copy number.

Run:

```python
from cysnet.copynumber import write_copy_number_outputs

copy_paths = write_copy_number_outputs(
    redox_marginals_path="results/MY_STUDY_redox_marginals.tsv",
    protein_matrix_path="pg_matrix.tsv",
    fasta_path="search_database.fasta",
    outdir="results",
    study_name="MY_STUDY",
    injected_mass_g=500e-9,
)
```

This writes:

```text
results/MY_STUDY_protein_copy_number.tsv
results/MY_STUDY_copy_substate_summary.tsv
```

The copy-number table reports, for each protein/sample:

```text
raw_intensity
scaled_mass_g
molecular_weight_da
protein_copies
fasta_cysteines
detected_cysteines
realised_substates
copy_limited
```

The copy-substate summary reports sample-level copy-number and state-space capacity summaries.

## 6. Copy-number-scaled constraints through Python

After running `constraints` and `copynumber`, run:

```python
from cysnet.copyconstraints import write_copy_constraint_outputs

copy_constraint_paths = write_copy_constraint_outputs(
    redox_marginals_path="results/MY_STUDY_redox_marginals.tsv",
    per_protein_constraints_path="results/MY_STUDY_per_protein_constraints.tsv",
    protein_copy_number_path="results/MY_STUDY_protein_copy_number.tsv",
    outdir="results",
    study_name="MY_STUDY",
)
```

This writes:

```text
results/MY_STUDY_copy_constraints.tsv
results/MY_STUDY_exact_substate_copies.tsv
results/MY_STUDY_copy_constraint_summary.tsv
results/MY_STUDY_copy_constraint_cohort_summary.tsv
```

The copy-constraint table reports:

```text
protein_copy_number
observed_coordinate
observed_space_log2
pruned_space_log2
collapse_extent_log2
excluded_fraction_by_boundary_priors
copy_limited_realised_observed_substates
min_required_observed_substates
oxiform_min_fraction
oxiform_max_fraction
oxiform_min_copies
oxiform_max_copies
fully_reduced_min_fraction
fully_reduced_max_fraction
fully_reduced_min_copies
fully_reduced_max_copies
compatible_fibre_copies
```

The oxiform copy bounds are sharp Fréchet union bounds over the observed cysteine coordinates:

```text
oxiform_min_fraction = max(p_i)
oxiform_max_fraction = min(1, sum(p_i))
```

where `p_i` are the observed cysteine oxidation marginals.

The fully reduced bounds are the complement:

```text
fully_reduced_min_fraction = max(0, 1 - sum(p_i))
fully_reduced_max_fraction = 1 - max(p_i)
```

When a complete-coverage protein has an exact observed-coordinate solution, `exact_substate_copies.tsv` reports the exact substate probability and copy number.

## Full Python workflow

A complete scripted workflow is:

```python
from cysnet.oxidia import write_oxidia_outputs
from cysnet.topology import write_topology_outputs
from cysnet.constraints import write_constraint_outputs
from cysnet.copynumber import write_copy_number_outputs
from cysnet.copyconstraints import write_copy_constraint_outputs

study = "MY_STUDY"
outdir = "results"

oxidia_paths = write_oxidia_outputs(
    light_path="UniMod_108_sites.tsv",
    heavy_path="UniMod_776_sites.tsv",
    outdir=outdir,
    study_name=study,
    sep="\t",
)

topology_paths = write_topology_outputs(
    redox_marginals_path=oxidia_paths["redox_marginals"],
    fasta_path="search_database.fasta",
    outdir=outdir,
    study_name=study,
    sep="\t",
)

constraint_paths = write_constraint_outputs(
    redox_marginals_path=oxidia_paths["redox_marginals"],
    protein_topology_path=topology_paths["protein_topology"],
    outdir=outdir,
    study_name=study,
    sep="\t",
)

copy_paths = write_copy_number_outputs(
    redox_marginals_path=oxidia_paths["redox_marginals"],
    protein_matrix_path="pg_matrix.tsv",
    fasta_path="search_database.fasta",
    outdir=outdir,
    study_name=study,
    injected_mass_g=500e-9,
)

copy_constraint_paths = write_copy_constraint_outputs(
    redox_marginals_path=oxidia_paths["redox_marginals"],
    per_protein_constraints_path=constraint_paths["per_protein_constraints"],
    protein_copy_number_path=copy_paths["protein_copy_number"],
    outdir=outdir,
    study_name=study,
)
```

Without a PG matrix, omit the final two blocks:

```text
write_copy_number_outputs
write_copy_constraint_outputs
```

CysNet will still produce valid redox, topology, constraint, sample-summary and cohort-summary outputs.

## Required Oxi-DIA site matrix structure

The reduced and oxidised site matrices should contain site identity columns and sample/run intensity columns.

Required site identity columns:

```text
Protein
Residue
Site
Sequence
```

Recommended metadata columns:

```text
Protein.Names
Gene.Names
```

Example:

```csv
Protein,Protein.Names,Gene.Names,Residue,Site,Sequence,S1,S2
P12345,Example protein,GENE1,C,45,ACDEFG,1000,1200
P12345,Example protein,GENE1,C,102,ACDCFG,500,600
```

The L and H files should use the same site identity structure. CysNet aligns the two channels by:

```text
Protein
Residue
Site
Sequence
```

## Required protein matrix structure

The optional PG / protein LFQ matrix should contain one protein-group identifier column and one or more sample intensity columns.

CysNet attempts to infer the protein-group identifier column from common names:

```text
Protein.Group
Protein.Ids
Protein.Group.Ids
Protein
protein_id
```

Sample columns should match the `sample_id` values in the CysNet redox marginal table.

Example:

```csv
Protein.Group,S1,S2
P12345,1000000,1200000
Q99999,500000,750000
```

The protein matrix is used only for copy-number scaling. It is not used to calculate cysteine oxidation.

## FASTA requirements

The FASTA should be the same FASTA used for the DIA-NN search.

CysNet parses UniProt-style headers such as:

```text
>sp|P12345|PROTEIN_HUMAN ...
```

and uses the accession:

```text
P12345
```

For non-UniProt headers, CysNet uses the first whitespace-delimited token.

For semicolon-delimited protein groups, CysNet resolves the first accession that matches the FASTA. If an isoform accession is not found, CysNet also tries the canonical accession by removing the suffix after `-`.

## Mathematical scope

CysNet operates on binary cysteine-redox substates.

For a protein with `R` cysteine coordinates, the binary state space has:

```text
2^R states
```

For a state `s` and site marginals `m_j`, CysNet defines:

```text
q_j(s) = m_j       if s_j = 1
q_j(s) = 1 - m_j   if s_j = 0
```

The sharp bounds for each state are:

```text
lower_s = max(0, sum_j q_j(s) - (R - 1))
upper_s = min_j q_j(s)
```

CysNet uses these bounds to classify each substate as:

```text
compatible
excluded
required
fixed-positive
```

If every state has identical lower and upper bounds, the observed-coordinate ensemble is exact. Otherwise, the feasible set is bounded.

For complete-coverage proteins, exact solutions can be reported as resolved observed-coordinate oxiform distributions.

For incomplete proteins, CysNet reports constraints over the measured cysteine-coordinate projection. It does not infer unmeasured cysteine states.

## Copy-number mathematical scope

For observed cysteine oxidation marginals `p_i`, CysNet bounds the fraction of molecules carrying at least one oxidised observed cysteine using sharp union bounds:

```text
lower oxiform fraction = max(p_i)
upper oxiform fraction = min(1, sum(p_i))
```

It bounds the fully reduced observed-coordinate fraction as:

```text
lower fully reduced fraction = max(0, 1 - sum(p_i))
upper fully reduced fraction = 1 - max(p_i)
```

When protein copy number is available, these fractions are multiplied by the estimated protein copy number:

```text
oxiform_min_copies = protein_copies * max(p_i)
oxiform_max_copies = protein_copies * min(1, sum(p_i))
```

These are bounds over the observed cysteine-coordinate projection. For incomplete proteins, they are not claims about unmeasured cysteine coordinates.

## Current milestone

The current repository implements the first tested CysNet v1 skeleton:

* validation of cysteine oxidation marginals;
* enumeration of binary observed-coordinate oxiform substates;
* sharp lower and upper bounds for each substate;
* compatible, excluded, required and fixed-positive state calls;
* exact versus bounded solution classification;
* oxiform union bounds;
* fully reduced bounds;
* Oxi-DIA L/H site import;
* redox marginal calculation;
* FASTA-derived cysteine topology;
* complete versus incomplete protein coverage;
* coverage-class and constraint summaries;
* sample-level resolved/constrained proteoform summaries;
* cohort-level unique resolved proteoform totals;
* resolved exact distributions;
* optional protein copy-number scaling;
* optional copy-number-scaled oxiform and fully reduced bounds;
* optional exact substate copy-number outputs;
* command-line workflows;
* Colab-native upload workflow;
* optional experimental widget;
* optional Streamlit app;
* unit tests for the core theorem cases.

This establishes the software foundation for extending CysNet from direct marginal inputs to full protein-level redox tables, FASTA-derived cysteine topology, copy-number-scaled oxiform constraints and later full oxiform identity/weight interpretation.
