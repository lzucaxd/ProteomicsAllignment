# Preprocessing Overview

## Purpose

The preprocessing layer prepares benchmark-ready inputs from raw CPTAC and
CCLE gene matrices. It is fully separated from method implementations and
benchmark evaluation.

## What It Does

1. **Load data**: CPTAC gene matrices (per study), CCLE corrected gene matrix
2. **Build metadata**: Standardized sample_meta (sample_id, domain, condition, study_id)
3. **Build shared feature space**: Gene intersection, prevalence filtering, SD filtering
4. **Task-specific subsets**: Subtype-balanced CPTAC, fixed CCLE cell lines, BvL tissue mapping
5. **Save benchmark-ready matrices**: Combined matrices + metadata CSVs

## Preprocessing Rules by Method

| Method | Feature level | Imputation | Standardization | Notes |
|--------|--------------|------------|-----------------|-------|
| Raw | gene | none | none | Intersection only |
| Bridge-aware | gene | none | none | Bridge extraction from msstats_input.tsv |
| Celligner | gene | within-domain median | z-score | 70% prevalence, SD > 0.01 |

## Protein vs Gene

- **Native-domain inference**: May use protein-level MSstatsTMT summaries
- **Representation methods**: Operate on gene-level matrices
- **Never mixed silently**: Config specifies `feature_level: gene` or `feature_level: protein`

## Running

```bash
python scripts/run_preprocessing.py --config configs/preprocessing/default.yaml
```

Outputs go to `data/processed/`.
