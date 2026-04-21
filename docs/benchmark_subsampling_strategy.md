# Benchmark Subsampling Strategy

## Guiding Principles

1. Do not randomly subsample unless absolutely necessary
2. If randomization is used, make it seeded and reproducible
3. Reduce confounding from study / mixture / plex / imbalance where possible
4. Document rationale for every subset decision

## Task A: Breast Subtype

### CPTAC (PDC000120)

**Strategy**: Mixture-balanced subset

CPTAC breast TMT experiments are organized into mixtures (plexes). PAM50
subtypes are not uniformly distributed across mixtures, creating
subtype-mixture confounding. The mixture-balanced strategy selects samples
to minimize this confound.

**Implementation**: `scripts/benchmark/subset_strategies.R::mixture_balanced_subset()`

**Current subset**:
- Luminal (LumA + LumB): ~46 samples
- Basal: ~16 samples
- Documented in `reports/benchmark_master/subtype_subset_summary.tsv`

### CCLE

**Strategy**: Fixed cell line list

CCLE design is constrained — limited breast cell lines with known PAM50-like
subtype annotations. Cell lines selected based on published subtype classification.

- Basal: HCC70, HCC1806, HCC1143, MDA-MB-468 (4 lines)
- Luminal: CAMA-1, MCF7, T-47D, ZR-75-1 (4 lines)

**Limitation**: Small sample size; any statistical inference on CCLE subtype
alone has very low power. Cross-domain comparison relies on directional
agreement rather than significance.

## Task B: Breast vs Lung

### CPTAC

**Strategy**: All available samples from PDC000120 (Breast) and PDC000153 (Lung)

No subsampling needed for the CPTAC side — the two studies are naturally
separate with no shared plexes, so there is no cross-study confounding within
CPTAC. The biological signal (cancer type) is strong.

### CCLE

**Strategy**: All cell lines annotated as Breast or Lung tissue of origin

Uses `data/ccle_peptide/sample_info_ccle.csv` for tissue mapping.
No subsampling; includes all available lines for each tissue.

## Reproducibility

All subset decisions are:
- Encoded in YAML configs (`configs/tasks/*.yaml`)
- Deterministic (no random seed needed for current strategies)
- Documented in this file and in task-specific metadata CSVs
