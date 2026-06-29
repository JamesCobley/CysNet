# CysNet

CysNet is a theorem-constrained workflow for inferring oxiform constraints from bottom-up cysteine redox proteomics.

CysNet treats site-level cysteine oxidation values as marginals over an unobserved distribution of cysteine-redox proteoform substates. Given site-level redox measurements, FASTA-derived cysteine topology and optional protein copy-number estimates, CysNet reports compatible, excluded, required, exact and bounded oxiform substates.

CysNet does not over-resolve incomplete bottom-up data (i.e., it only does what the mathematics allows). Complete-coverage proteins may collapse to exact oxiform ensemble distributions; incomplete proteins return the exact constraint class supported by the measured cysteine coordinates.
