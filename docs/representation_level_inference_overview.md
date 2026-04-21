# Representation-Level Inference Overview (Level 2)

## Purpose

Evaluate each harmonization method's transformed matrix using the **same
statistical framework** — limma-based differential abundance — so methods
are comparable on a level playing field.

## Why limma (Not MSstatsTMT)

Once data has been transformed by a harmonization method (bridge correction,
Celligner alignment, etc.), the original TMT experimental design is no longer
intact in the same way. Using MSstatsTMT on transformed data would be
statistically inappropriate because the variance model assumes TMT-specific
structure that no longer holds.

Instead, we use limma on the gene-level matrix, which makes fewer assumptions
about the data generation process and treats samples as independent observations.

## What It Produces

For each method × task × domain:
- Per-domain limma DA results (logFC, t-stat, p-value, adj.p.value)
- Cross-domain FC agreement (correlation, same-direction fraction)
- Marker direction checks

## Metrics Computed

| Metric | Description | Ideal |
|--------|-------------|-------|
| FC correlation | Pearson r between CPTAC and CCLE logFCs | High positive |
| Same-direction fraction | % of genes with concordant FC sign | High |
| Median FC difference | CPTAC logFC minus CCLE logFC, median | Near 0 |
| Marker concordance | Fraction of markers with matching direction | High |

## Running

```bash
python scripts/run_benchmark.py --methods raw bridge_shift celligner
```

## Output Location

`reports/benchmark_master/benchmark_results/<method>/<task>/representation_da/`
