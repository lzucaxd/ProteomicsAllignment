# Diagnostics Master Summary

## Purpose

Diagnostics are first-class benchmark outputs that validate assumptions
and identify issues before drawing conclusions from DA or metrics.

## Diagnostic Categories

### 1. Distribution Checks
- Raw vs log boxplots by domain
- Per-sample intensity distributions
- **What to expect**: CPTAC and CCLE on comparable log2 scales
- **What worries us**: Bimodal distributions, extreme outlier samples

### 2. Bridge / Reference Channel QC
- Bridge coverage per protein
- Bridge offset distribution
- Extreme offset proteins
- **What to expect**: Unimodal offset distribution centered around platform difference
- **What worries us**: Large fraction of proteins without bridge data

### 3. Effect-Size and Spread Summaries
- FC scatter plots (CPTAC logFC vs CCLE logFC)
- FC agreement statistics (correlation, same-direction fraction)
- **What to expect**: Positive correlation, majority same direction
- **What worries us**: Low correlation, method flattening biological effects

### 4. PCA / Structure Checks
- PCA colored by domain, condition, plex
- Domain R² on top PCs
- **What to expect**: Domain effect reduced after harmonization
- **What worries us**: Domain still dominates PC1, or biology lost

### 5. Method-Specific QC
- Bridge-aware: offset/scale factor distributions, QC flags
- Celligner: DE genes, pre/post PCA, gene set reduction
- **Location**: `reports/benchmark_master/methods/<method>/`

## Implementation

- R diagnostics: `scripts/benchmark/diagnostics.R`
- Python structure plots: `src/harmonize/benchmark/plots/structure.py`
- Python FC scatter: `src/harmonize/benchmark/plots/diagnostics.py`

## Output Location

`reports/benchmark_master/diagnostics/`
