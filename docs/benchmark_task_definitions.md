# Benchmark Task Definitions

## Task A: Breast Subtype (Fine-grained biology)

**Contrast**: Luminal vs Basal

**Purpose**: Test whether harmonization preserves fine-grained molecular subtype
differences that are biologically important within breast cancer.

**CPTAC**: PDC000120 with mixture-balanced subset to reduce subtype-plex confounding.
PAM50 subtypes: LumA and LumB grouped as Luminal; Basal as Basal.

**CCLE**: Fixed cell lines.
- Basal: HCC70, HCC1806, HCC1143, MDA-MB-468
- Luminal: CAMA-1, MCF7, T-47D, ZR-75-1

**Markers**: ESR1, PGR, GATA3, FOXA1 (luminal); EGFR, KRT5, KRT17, FOXC1 (basal)

**Block order** (interleaved for harmonization comparison):
CPTAC Luminal | CCLE Luminal | CPTAC Basal | CCLE Basal

**Config**: `configs/tasks/breast_subtype.yaml`

---

## Task B: Breast vs Lung (Coarse biology)

**Contrast**: Breast vs Lung

**Purpose**: Test whether harmonization preserves cancer-type level differences,
which are a stronger biological signal than subtype.

**CPTAC**: PDC000120 (Breast) + PDC000153 (Lung)

**CCLE**: Breast and Lung tissue-of-origin cell lines from sample_info_ccle.csv

**Markers**: NKX2-1, SFTPB, NAPSA (lung); GATA3, FOXA1, ESR1 (breast);
EGFR, ERBB2, CDH1, KRT7, MUC1 (shared/variable)

**Block order** (interleaved):
CPTAC Breast | CCLE Breast | CPTAC Lung | CCLE Lung

**Config**: `configs/tasks/breast_vs_lung.yaml`

---

## Subset Strategies

See existing `reports/benchmark_master/benchmark_subsampling_strategy.md` for
detailed statistical rationale on mixture-balanced CPTAC subsets and CCLE
design limitations.
