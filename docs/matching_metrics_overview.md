# Matching Metrics Overview (Level 5)

## Purpose

Evaluate whether harmonization improves the ability to correctly match
CPTAC tumor samples to biologically similar CCLE cell lines.

## Status

**Partially implemented / scaffold.** Core distance computation and top-k
matching are implemented. Cross-method NN consistency is implemented.
More sophisticated matching (marker consistency in matched pairs, rank-based
metrics) can be added.

## Metrics

### Top-k Same-Condition Match Rate

For each CPTAC tumor, find its k nearest CCLE neighbors (by correlation
distance in the harmonized space). Fraction of tumors where at least one
of the top-k neighbors shares the same condition (e.g., Breast tumor
matched to Breast cell line).

### Nearest-Neighbor Consistency

Compare NN assignments across methods: for each tumor, what fraction of
its top-k CCLE neighbors are the same across two method representations?
High overlap = methods agree on which cell lines are most similar.

## What Is Needed

- Reliable tissue/subtype annotations for both CPTAC and CCLE samples
- Sufficient samples per condition in both domains
- Meaningful distance metric in the shared feature space

## Implementation

`src/harmonize/benchmark/metrics/matching.py`

Uses: scipy (cdist), numpy, pandas
