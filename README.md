# CysNet

CysNet is a theorem-constrained workflow for inferring oxiform constraints from bottom-up cysteine redox proteomics.

CysNet treats site-level cysteine oxidation values as marginals over an unobserved distribution of cysteine-redox proteoform substates. Given site-level redox measurements, FASTA-derived cysteine topology and optional protein copy-number estimates, CysNet reports compatible, excluded, required, exact and bounded oxiform substates.

CysNet does not over-resolve incomplete bottom-up data (i.e., it only does what the mathematics allows). Complete-coverage proteins may collapse to exact oxiform ensemble distributions; incomplete proteins return the exact constraint class supported by the measured cysteine coordinates.
## Installation and first test

CysNet can be installed directly from the GitHub repository into a local Python environment or a Google Colab notebook.

### Clone the repository

```bash
git clone https://github.com/JamesCobley/CysNet.git
cd CysNet
```

### Install CysNet in editable mode

```bash
python -m pip install -e ".[dev]" --upgrade --upgrade-strategy only-if-needed
```

This installs CysNet together with the development dependencies required to run the tests.

### Run the theorem tests

```bash
pytest -q
```

Expected output:

```text
7 passed
```

These tests check the core theorem behaviour, including boundary exclusion, exact oxiform resolution, bounded feasible ensembles, oxiform union bounds, fully reduced bounds and invalid input handling.

## Command-line examples

CysNet currently exposes a minimal command-line interface for testing cysteine-redox marginals directly.

### Exact observed-coordinate solution

```bash
cysnet 0 0 0.25
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
cysnet 0 0.25 0.25
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

## Current milestone

The current repository implements the first tested CysNet theorem core:

* validation of cysteine oxidation marginals;
* enumeration of binary observed-coordinate oxiform substates;
* sharp lower and upper bounds for each substate;
* compatible, excluded, required and fixed-positive state calls;
* exact versus bounded solution classification;
* oxiform union bounds;
* fully reduced bounds;
* a minimal command-line interface;
* unit tests for the core theorem cases.

This establishes the software foundation for extending CysNet from direct marginal inputs to protein-level redox tables, FASTA-derived cysteine topology and copy-number-scaled oxiform constraints.
