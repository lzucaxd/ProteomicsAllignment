# Structure Metrics Overview (Level 4)

## Purpose

Quantify whether a harmonized representation preserves biological structure
while reducing technical domain effects.

## Metrics

### PCA Variance Decomposition

For each label (domain, condition), compute ANOVA R² on top PCs:
- `domain_r2_pc1`, `domain_r2_pc2`: How much of PC1/PC2 variance is
  explained by domain (CPTAC vs CCLE). Lower is better for harmonization.
- `condition_r2_pc1`, `condition_r2_pc2`: How much is explained by biology.
  Higher is better.
- `*_r2_top5`: Average R² across top 5 PCs.

### kNN Purity

For each sample, fraction of k-nearest neighbors sharing the same label.
- `knn_purity_domain`: High means domains cluster separately (bad for mixing).
- `knn_purity_condition`: High means biology is preserved (good).

### Silhouette Score

- `silhouette_domain`: Positive means domain clusters are well-separated (bad).
- `silhouette_condition`: Positive means biological groups are well-separated (good).

### Classification Accuracy

Logistic regression on top PCs:
- `classification_acc_domain`: How well can a classifier predict domain?
  Near 0.5 = good mixing. Near 1.0 = domains still separate.
- `classification_acc_condition`: How well can it predict biology?
  Higher = biological structure preserved.

## Ideal Outcome

| Metric | Raw (baseline) | Good harmonization |
|--------|---------------|-------------------|
| domain_r2_pc1 | high | low |
| condition_r2_pc1 | moderate | high (or maintained) |
| knn_purity_domain | high | low |
| knn_purity_condition | moderate | maintained or improved |
| classification_acc_domain | ~1.0 | closer to 0.5 |
| classification_acc_condition | moderate | maintained or improved |

## Implementation

`src/harmonize/benchmark/metrics/structure.py`

Uses: sklearn (PCA, LogisticRegression, KNeighborsClassifier, silhouette_score)
