# `data/scripts/` — exploratory and legacy helpers

**This directory is not the reproducible “manifest → PSM → gene matrix” pipeline.**

That pipeline lives at the **`data/`** root (`run_pipeline_per_manifest.sh`, `pdc_manifest_downloader.py`, `pdc_psm_to_msstatsTMT_protein_matrix.R`). **Start here instead:** [`../../pipeline/psm_to_gene_matrix/README.md`](../../pipeline/psm_to_gene_matrix/README.md).

## What lives here

- **Subtype / Luminal–Basal exploratory DA** (MSstatsTMT on protein summaries, limma on `gene_matrix.csv`, Venns, marker panels, …).
- **Benchmark v1** helpers and diagnostics (`build_benchmark_v1_artifacts.py`, bridge QC, QQ notes, …).
- **CCLE** one-offs (Table S2 paths, MSI vs MSS, reporter QC).

These scripts were developed for **papers, slides, and diagnostics**. They are **tracked** so provenance stays grep-able from `reports/*.md` and `PROJECT_REPORT.md`, but **new** exploratory code should go under **`scripts/exploratory/`** (see [`../../scripts/exploratory/README.md`](../../scripts/exploratory/README.md)) or **`notebooks/exploratory/`** so the main pipeline folder stays easy to explain.

## When you need MSstatsTMT vs limma for *subtype* work

See **`docs/INFERENCE_BASELINES.md`** — benchmark-native vs representation-level paths are separate from most files in this folder.
