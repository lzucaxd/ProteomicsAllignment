# Repository layout plan

**Navigation:** for a contributor-facing map of historical vs canonical paths and **naming rules for new files**, see **[`docs/NAMING_AND_PATHS.md`](docs/NAMING_AND_PATHS.md)**.

## Principles (handoff)

- **Code** (`src/`, `scripts/`) stays separate from **large generated data** (`data/results/`, `data/pdc_psm/`) and **presentation exports** (`presentation_materials/`).
- **Configs** (`configs/`) drive paths where possible; avoid new hardcoded absolute paths.
- **Reports** under `reports/benchmark_master/` remain the benchmark’s “output home”; an index lives under **`reports/final_report/`** without relocating heavy files (avoids breaking R scripts).

## Actual layout (implemented — conservative)

This repository **does not** use a full physical re-rooting (e.g. moving all CPTAC drivers out of `data/`) in this pass, because hundreds of R/shell/Python paths assume the current `data/` layout. Instead:

| Proposed in brief | Actual location | Notes |
|-------------------|-----------------|-------|
| `environment/` | **`environment/README.md`** | Single doc for Python + R setup |
| `results/benchmark/` | **`reports/benchmark_master/benchmark_results/`** | Renaming would break docs and scripts |
| `data/raw/` | **`data/pdc_psm/`**, manifests | Already established names |
| `data/processed/` | **`data/processed/union/`**, `data/processed/methods/` | Overnight v2 contract |
| `notebooks/exploratory/` | **`notebooks/exploratory/`** | New; holds non-pipeline notebooks |

## `src/` vs `scripts/`

| Location | Contents |
|----------|----------|
| **`src/harmonize/`** | Installable Python library: preprocessing, benchmark orchestration helpers, method classes |
| **`scripts/`** | Runnable entrypoints: shell/R/Python drivers, often thin wrappers around `data/*.R` or `src/` |

## `reports/` vs `results/`

Benchmark CSVs and figures are under **`reports/benchmark_master/`** by design. A **`reports/final_report/`** folder provides a **curated index** (`README.md`, `FIGURE_MANIFEST.md`) pointing into `final_tables/` and presentation scripts—**files are not duplicated** unless a future maintainer explicitly copies assets for Zenodo.

## `archive/`

Use for **legacy duplicates** and retired trees. See **`archive/README.md`**. Prefer **move + log** over delete; record moves in **`CLEANUP_LOG.md`**.
