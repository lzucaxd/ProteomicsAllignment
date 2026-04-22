# Data provenance

This repository’s **clone name** is `ProteomicsAllignment` (historical spelling on disk). Paths below are **relative to the repo root**.

## Raw data (not tracked in git)

Large downloads and intermediates are **gitignored** (see root `.gitignore`).

| Location | Source | Description |
|----------|--------|-------------|
| `data/pdc_psm/{study_id}/` | PDC | Downloaded CPTAC PSM tables per study |
| `data/manifests/` | PDC | **Local only:** file manifest CSVs (not in git — URLs expire). See **`data/manifests/README.md`** + **`example_pdc_file_manifest.csv`**. |
| `data/ccle_peptide/` | DepMap / Broad-style inputs | CCLE peptide-level reporter-ion tables and helpers |
| `data/results/{study_id}/` | This pipeline | Per-study MSstatsTMT outputs (`gene_matrix.csv`, etc.) — **ignored** by default due to size |

**Note:** Some documentation refers to an idealized `data/pdc/PDC000120/` layout; in **this** repo, CPTAC raw downloads live under **`data/pdc_psm/`** (see [`PIPELINE_README.md`](PIPELINE_README.md)).

## Processed data (tracked selectively)

| Location | Produced by | Description |
|----------|-------------|-------------|
| `data/processed/union/*.csv` | `scripts/run_preprocessing.py` | Per-task **shared gene matrices** + `sample_meta_*.csv` |
| `data/processed/methods/{method}/transformed_*.csv` | `scripts/benchmark/regenerate_methods_union.py` (Step 2 of overnight) | Harmonized representations for **raw**, **bridge_shift**, **bridge_scale**, **celligner** |
| `data/processed/*annotation*` | `scripts/benchmark/process_ccle_annotations_v2.R` | CCLE subtype / BvL labels used in tasks |

Some older duplicate matrices may exist under `data/processed/` **and** `data/processed/union/`; **`data/processed/union/`** is what **`run_overnight_v2.sh`** uses.

## How to regenerate

### 1) PSM → gene matrices (once per study; long-running)

From repo root, follow **[`../pipeline/psm_to_gene_matrix/README.md`](../pipeline/psm_to_gene_matrix/README.md)** (short map) then **[`data/PIPELINE_README.md`](PIPELINE_README.md)** (full narrative). Typical entrypoints:

```bash
cd data
# Option A: per-manifest pipeline (see PIPELINE_README for prerequisites)
./run_pipeline_per_manifest.sh

# Option B: batch studies with a local CPTAC mirror (requires CPTAC_LOCAL_MIRROR)
./run_batch_studies.sh
```

The main R driver is **`data/pdc_psm_to_msstatsTMT_protein_matrix.R`**.

### 2) Union matrices + benchmark (after `gene_matrix.csv` exist)

```bash
# From repo root — full benchmark (hours; permutations + ceilings are slow)
bash scripts/benchmark/run_overnight_v2.sh

# Inspect canonical table
column -t -s, reports/benchmark_master/benchmark_results/comparison_summary.csv | head
```

## CCLE subtype annotations (benchmark v2)

Curated luminal / basal labels for **25** breast lines (**8** luminal + **17** basal), HER2-enriched excluded, **CAL120** merged — processed in:

`scripts/benchmark/process_ccle_annotations_v2.R`

See also `reports/benchmark_master/` diagnostics and `configs/tasks/breast_subtype.yaml`.

## Further reading

- [`scripts/preprocessing/README.md`](../scripts/preprocessing/README.md) — preprocessing narrative index  
- [`docs/METHODS.md`](../docs/METHODS.md) — harmonization methods  
- [`docs/REPO_AUDIT.md`](../docs/REPO_AUDIT.md) — machine-generated tree snapshot  
