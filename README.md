# CysNet

CysNet is a theorem-constrained workflow for inferring oxiform constraints from bottom-up cysteine redox proteomics.

CysNet treats site-level cysteine oxidation values as marginals over an unobserved distribution of cysteine-redox proteoform substates. Given site-level redox measurements, FASTA-derived cysteine topology and optional protein copy-number estimates, CysNet reports compatible, excluded, required, exact and bounded oxiform substates.

CysNet does not over-resolve incomplete bottom-up data. It only returns what the mathematics supports. Complete-coverage proteins may collapse to exact oxiform ensemble distributions; incomplete proteins return the exact constraint class supported by the measured cysteine coordinates.

## Current status

CysNet currently implements:

* theorem-constrained state-bound inference from cysteine redox marginals;
* Oxi-DIA light/heavy site matrix import;
* redox marginal calculation from reduced and oxidised cysteine channels;
* FASTA-derived cysteine topology bookkeeping;
* complete versus incomplete coverage classification;
* command-line workflows;
* a simple Jupyter/Colab upload widget;
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

The protein LFQ matrix is not used to calculate redox marginals. Protein abundance or copy-number information is used later for abundance-scaled or copy-number-scaled oxiform constraints.

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

### Install CysNet with notebook widget support

```bash
python -m pip install -e ".[dev,widget]" --upgrade --upgrade-strategy only-if-needed
```

### Install CysNet with Streamlit app support

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

## Command-line usage

CysNet currently exposes three command-line workflows:

```text
cysnet theorem
cysnet oxidia-sites
cysnet topology
```

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

The topology summary reports sample-level coverage and state-space bookkeeping, including:

```text
number of detected protein groups
number of complete protein groups
number of incomplete protein groups
number of FASTA cysteines
number of detected cysteines
cysteine coverage percentage
log10 summed full state space
log10 summed observed-coordinate state space
```

## Full command-line workflow

A typical CysNet v1 workflow is:

```bash
git clone https://github.com/JamesCobley/CysNet.git
cd CysNet

python -m pip install -e ".[dev]" --upgrade --upgrade-strategy only-if-needed

cysnet oxidia-sites \
  --light UniMod_108_sites.tsv \
  --heavy UniMod_776_sites.tsv \
  --study MY_STUDY \
  --out results

cysnet topology \
  --redox-marginals results/MY_STUDY_redox_marginals.tsv \
  --fasta search_database.fasta \
  --study MY_STUDY \
  --out results
```

## Notebook upload widget

CysNet includes a simple Jupyter/Colab upload widget for users who do not want to use command-line file paths.

Install with widget support:

```bash
python -m pip install -e ".[dev,widget]" --upgrade --upgrade-strategy only-if-needed
```

Launch in a notebook:

```python
from cysnet.widget import launch_oxidia_widget

launch_oxidia_widget()
```

The widget asks for:

```text
study name
L / Light / UniMod_108 site matrix
H / Heavy / UniMod_776 site matrix
FASTA file used for the DIA-NN search
```

The widget then computes:

```text
site-level percent oxidised
detected-site coverage
sample-level redox summary
CysNet-ready redox marginals
FASTA-derived cysteine topology
complete/incomplete protein coverage
topology summary
```

It also creates a ZIP archive containing all generated outputs.

### Colab note

If Colab cannot find the package after cloning, add the source path manually:

```python
import sys
sys.path.insert(0, "/content/CysNet/src")

from cysnet.widget import launch_oxidia_widget
launch_oxidia_widget()
```

## Streamlit app

CysNet also includes an optional Streamlit upload app.

Install with app support:

```bash
python -m pip install -e ".[dev,app]" --upgrade --upgrade-strategy only-if-needed
```

Run locally:

```bash
streamlit run src/cysnet/app.py
```

The Streamlit app provides the same basic upload workflow:

```text
upload L
upload H
upload FASTA
enter study name
run CysNet
download outputs
```

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

For incomplete proteins, CysNet reports constraints over the measured cysteine-coordinate projection. It does not infer unmeasured cysteine states.

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
* command-line workflows;
* notebook upload widget;
* optional Streamlit app;
* unit tests for the core theorem cases.

This establishes the software foundation for extending CysNet from direct marginal inputs to full protein-level redox tables, FASTA-derived cysteine topology and copy-number-scaled oxiform constraints.

