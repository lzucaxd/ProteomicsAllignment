# Benchmark pipeline (CPTAC ↔ CCLE harmonization)

This directory contains the **representation-level benchmark**: same tasks, same genes, same limma protocol per domain — different harmonization methods.

**Orchestrator:** `run_overnight_v2.sh`  
**Canonical results table:** `reports/benchmark_master/benchmark_results/comparison_summary.csv`  
**Disconnect table:** `reports/benchmark_master/benchmark_results/disconnect_scores.csv`

**End-to-end (data → matrices → this script):** [`docs/HOW_TO_RUN_EVERYTHING.md`](../../docs/HOW_TO_RUN_EVERYTHING.md)

---

## Quick start

```bash
# From repo root — full pipeline (permutations + split-half ceilings are SLOW; ~hours)
bash scripts/benchmark/run_overnight_v2.sh
```

Logs: `reports/benchmark_master/logs/overnight_v2_*.log`  
Diagnostics: `reports/benchmark_master/diagnostics/`

**Prerequisites**

- **R** (≥ 4.1 recommended): `limma`, `data.table`, `ggplot2`, MSstatsTMT stack as needed for your local preprocessing (see root `install_r_packages.R`).
- **Python** (≥ 3.10): `pandas`, `numpy`, `pyyaml`, `scikit-learn`; repo uses `PYTHONPATH=src`.
- Input matrices: `data/processed/union/` and `data/results/*/gene_matrix.csv` per `configs/`.

There is **no** `--skip-slow` flag in `run_overnight_v2.sh` today; to iterate quickly, comment out Steps 7–8 locally or run individual R/Python steps below.

---

## What the benchmark measures

Two complementary axes:

1. **Geometry / mixing** — On a **fixed PCA basis** fit on **raw** concatenated samples, do domains overlap after harmonization? (R² on PC1 for domain indicator, silhouette, kNN purity, etc.)
2. **Differential abundance agreement** — For the **same contrast** within each domain (CPTAC vs CCLE separately), does the method preserve **per-gene logFC** relationships across domains? (Pearson correlation, same-direction fraction, marker checks, stratified FC, permutation nulls, concordance ceilings.)

**Key point:** improving (1) does not guarantee improving (2); see `disconnect_scores.csv` and `docs/METHODS.md`.

---

## Pipeline steps (`run_overnight_v2.sh`)

| Step | Script | What it does |
|------|--------|----------------|
| **0** | `process_ccle_annotations_v2.R` | Normalize CCLE labels; CAL120 merge; outputs under `data/processed/` consumed by tasks |
| **1** | `scripts/run_preprocessing.py` | Build per-task **union** matrices + metadata → `data/processed/union/` |
| **1b** | `preflight_diagnostics.R`, `compute_intersection_masks.R` | Coverage audits, intersection gene lists → `reports/benchmark_master/diagnostics/` |
| **2** | `regenerate_methods_union.py` | For each method, write `data/processed/methods/{method}/transformed_{task}.csv` |
| **2b** | `run_structure_batch.py` | Structure metrics on fixed PCA basis → feeds comparison table |
| **3** | `run_all_limma_da.R` | **16** limma runs: 4 methods × 2 tasks × 2 domains → `representation_da/{cptac,ccle}/da_limma_result.csv` |
| **4** | `compute_cross_domain_metrics.R` | FC correlation, same-direction fraction, marker sanity → wide metrics |
| **5** | `compute_stratified_fc.R` | Stratify genes by within-domain significance; per-stratum FC correlation → diagnostics CSVs |
| **6** | `generate_volcanos.R` | Volcano plots (optional QC / slides) |
| **7** | `run_permutation_nulls.R` | Label permutation null for FC correlation (default **1000** permutations) |
| **8** | `run_concordance_ceilings.R` | CPTAC + CCLE split-half concordance ceilings (**200** splits default) |
| **9** | `run_fast_calibration.R` | Calibrated ratios, biology destruction summaries, residual dependence |
| **10** | `compute_disconnect_scores.R` | `disconnect_scores.csv` — geometry vs DA improvement vs biology cost |
| **11** | `assemble_comparison_table.py` | Merges everything → **`comparison_summary.csv`** |
| **12** | `generate_meeting_figures.R` | Exports under `reports/benchmark_master/meeting/` |

Supporting code: `evaluation_helpers.R`, `harmonize_paths.R`, `benchmark_runner.R`, `native_domain_da.R`, etc., are used by individual steps or older entrypoints.

---

## Configuration

| File | Role |
|------|------|
| `configs/preprocessing/default.yaml` | Gene-space construction, study paths |
| `configs/preprocessing/union.yaml` | Union / filtering parameters |
| `configs/tasks/breast_subtype.yaml` | Task A: luminal vs basal |
| `configs/tasks/breast_vs_lung.yaml` | Task B: breast vs lung |
| `configs/methods/*.yaml` | Method-specific options (e.g. Celligner) |

---

## Adding a new harmonization method

1. Implement a transform that reads **`data/processed/union/shared_gene_matrix_{task}.csv`** (and metadata) and writes **`data/processed/methods/{your_method}/transformed_{task}.csv`** with the **same sample columns** and gene identifier column as other methods.
2. Register the method in **`regenerate_methods_union.py`** (and any YAML it reads).
3. Re-run from **Step 2** onward (or the full overnight script).
4. Confirm a new **`method`** row appears in `comparison_summary.csv`.

The evaluation stack is **method-agnostic**: limma, metrics, permutations, and ceilings do not depend on how the representation was produced.

---

## Adding a new benchmark task

1. Add **`configs/tasks/{task}.yaml`** describing CPTAC studies, contrasts, and CCLE column selections.
2. Extend **`scripts/run_preprocessing.py`** / union config so the union builder emits `shared_gene_matrix_{task}.csv` and `sample_meta_{task}.csv`.
3. Re-run **Steps 1–12**.

---

## Key output files

| Path | Description |
|------|-------------|
| `reports/benchmark_master/benchmark_results/comparison_summary.csv` | **Master table** — all methods × tasks × metrics |
| `reports/benchmark_master/benchmark_results/disconnect_scores.csv` | Geometry vs DA disconnect |
| `reports/benchmark_master/benchmark_results/{method}/{task}/representation_da/{cptac,ccle}/da_limma_result.csv` | Per-gene statistics |
| `reports/benchmark_master/diagnostics/` | Stratified FC, QC tables, preflight |
| `reports/benchmark_master/logs/` | Overnight logs (gitignored) |

Bridge offsets used by **bridge_shift** / **bridge_scale** are derived in **`extract_bridge_summaries.R`** + **`bridge_aware_correction.R`** from **`msstats_input.tsv`** Norm rows (not from domain medians of sample abundances).

---

## Presentation / report figures

- **Meeting exports:** Step 12 above.  
- **Paper-style report panels:** `scripts/presentation/generate_report_figures.R` → `presentation_materials/figures/report/` (regenerate locally; directory may be gitignored — see `.gitignore`).

---

## Further reading

- [`docs/BENCHMARK.md`](../../docs/BENCHMARK.md) — short index + deep links  
- [`docs/END_TO_END_TECHNICAL_REPORT.md`](../../docs/END_TO_END_TECHNICAL_REPORT.md) — long-form methods + metrics narrative  
- [`docs/BENCHMARK_V2_AND_PRESENTATION.md`](../../docs/BENCHMARK_V2_AND_PRESENTATION.md) — slide checklist  
