# Native-Domain Inference Overview (Level 1)

## Purpose

Establish the statistically grounded baseline DA results using the **original
TMT experimental design** before any cross-domain harmonization is applied.

## Why Separate

The original TMT design (mixtures, reference channels, within-plex normalization)
is only meaningful in the native data. Once data is transformed by harmonization
methods, the TMT design structure is no longer intact. Therefore:

- **Native-domain inference** uses MSstatsTMT (protein-level, TMT-aware)
- **Representation-level inference** uses limma (gene-level, matrix-based)

These must never be mixed or compared as if they were the same thing.

## What It Produces

For each domain and task:

| Domain | Tool | Input | Output |
|--------|------|-------|--------|
| CPTAC | MSstatsTMT `groupComparisonTMT` | msstats_input.tsv | Protein-level DA tables |
| CCLE | limma on gene matrix | gene_matrix.csv | Gene-level DA tables |

## Use in Benchmark

Native-domain DA results serve as:
1. **Directional anchors**: Expected fold-change directions for markers
2. **Size anchors**: Approximate magnitude of biological effects
3. **Baseline comparison**: Are representation-level results consistent with
   what native inference shows?

## Running

```bash
python scripts/run_native_baselines.py
```

## Output Location

`reports/benchmark_master/native_domain_da/`
