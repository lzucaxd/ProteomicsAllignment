# Benchmark Overview

## Philosophy

The benchmark answers: **Can a harmonization method improve cross-domain
comparability while preserving biology?**

A good method should:
- Reduce residual domain effect
- Preserve lineage / cancer-type structure
- Preserve subtype structure
- Preserve marker behavior
- Improve cross-domain effect-size agreement
- Improve tumor-cell line matching quality

## Evaluation Levels

### Level 1 — Native-domain baseline inference

- Uses original, non-aligned CPTAC and CCLE data
- MSstatsTMT for CPTAC (native TMT inference)
- limma on gene matrix for CCLE
- Produces per-domain DA tables as directional/size anchors
- **Separate from representation-level comparison**

### Level 2 — Representation-level feature comparison

- For each transformed matrix (Raw, Bridge, Celligner, ...)
- limma-based DA within each domain
- Cross-domain FC agreement (correlation, same-direction fraction)
- Marker gene summaries and direction checks

### Level 3 — Marker profile evaluation

- Method-agnostic polished profile plots
- Fixed sample order and y-axis across methods
- Faint individual points + bold block medians
- 3-4 markers per page, meeting/paper quality

### Level 4 — Structure and geometry evaluation

- PCA variance decomposition (domain R², condition R²)
- kNN purity by domain and condition
- Silhouette scores
- Domain/condition classification accuracy
- PCA/UMAP visualizations

### Level 5 — Matching / retrieval evaluation

- Top-k same-lineage match rate
- Top-k same-subtype match rate
- Nearest-neighbor consistency across methods
- **Currently scaffolded, partially implemented**

## Critical Statistical Split

| Layer | When to use | Statistical tool | Notes |
|-------|------------|-----------------|-------|
| Native-domain | Original TMT design intact | MSstatsTMT | CPTAC only |
| Representation-level | Transformed/aligned matrices | limma | All methods |

This distinction is encoded in: code structure, config, output directories, naming.
