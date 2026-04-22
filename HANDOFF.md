# Lab handoff — Proteomics alignment benchmark

**Audience:** A technically strong lab member who did not work on the repo day-to-day.

## What this project is

A **calibrated benchmark** for **cross-dataset harmonization** of clinical **CPTAC** and preclinical **CCLE** TMT proteomics on a **shared gene space**, with explicit **null** and **ceiling** calibration. Core question: do harmonization methods improve **cross-domain fold-change agreement**, or mainly **PCA geometry**?

## Benchmark tasks (current)

1. **Breast subtype** — Luminal vs Basal (CPTAC + CCLE), curated subtype labels.
2. **Breast vs lung** — Tissue contrast (CPTAC breast vs lung + CCLE arms).

Task definitions: `configs/tasks/breast_subtype.yaml`, `configs/tasks/breast_vs_lung.yaml`.

## Methods compared (overnight / v2)

| ID | Description |
|----|-------------|
| `raw` | No harmonization (baseline) |
| `bridge_shift` | Per-gene offset from **bridge channel** in `msstats_input.tsv` |
| `bridge_scale` | Bridge shift + MAD-based scaling |
| `celligner` | cPCA + MNN (global alignment) |

Definitions: **`docs/METHODS.md`**.

## Key outputs (where to look first)

| Output | Path |
|--------|------|
| Master comparison table | `reports/benchmark_master/benchmark_results/comparison_summary.csv` |
| Disconnect-style metrics | `reports/benchmark_master/benchmark_results/disconnect_scores.csv` |
| Consolidated numbers (Markdown) | `reports/benchmark_master/final_tables/BENCHMARK_NUMBERS_MASTER_FOR_LLM.md` |
| Index to report assets | **`reports/final_report/README.md`** and **`FIGURE_MANIFEST.md`** |

## Main scripts to run (in order)

1. **Environment:** `environment/README.md` → Python venv + `pip install -e .` + `Rscript install_r_packages.R`.
2. **Verify:** `python3 scripts/verify_repro_setup.py` (add `--require-data` when matrices exist).
3. **PSM → gene matrix (CPTAC/CCLE):** from `data/` — `data/PIPELINE_README.md`; typical: `./run_pipeline_per_manifest.sh` or `./run_batch_studies.sh`.
4. **Full benchmark:** **`bash scripts/run_benchmark.sh`** (wrapper for `scripts/benchmark/run_overnight_v2.sh`).
5. **Diagnostics only (subset):** `./scripts/run_diagnostics.sh preflight` or `… structure` or `… all`.
6. **Method registry smoke test (Python):** `./scripts/run_methods.sh` (not the same as overnight union matrices; see script header).

**Single narrative runbook:** **`docs/HOW_TO_RUN_EVERYTHING.md`**.

## What to read first

1. **`README.md`** (repo home)
2. **`docs/HOW_TO_RUN_EVERYTHING.md`**
3. **`data/PIPELINE_README.md`** (manifest → matrix)
4. **`scripts/benchmark/README.md`** (overnight steps)
5. **`docs/METHODS.md`**
6. **`REPO_AUDIT.md`** + **`REPO_LAYOUT_PLAN.md`** (layout and clutter policy)

## Data and reproducibility

- **Raw / large data are not in git.** Manifests expire; place fresh PDC exports under `data/manifests/` per **`data/manifests/README.md`** and **`EXPECTED_INPUTS.md`**.
- **Processed matrices** default under `data/results/` and `data/processed/union/` (see `data/README.md`).
- **Paper freeze:** tag a commit and optionally archive key CSVs externally (Zenodo); see `docs/CLEAN_CLONE_REPRODUCIBILITY.md`.

## Known caveats

- **Historical layout:** heavy CPTAC drivers still live under **`data/`**; `scripts/preprocessing/` is mostly documentation.
- **Two method paths:** Python `harmonize` registry vs overnight **R** union pipeline — both documented in `scripts/benchmark/README.md` and `HANDOFF_SANITY_CHECK.md`.
- **Optional Celligner** requires extra Python deps and a local `models/` checkout.
- **`run_pipeline_per_manifest.sh`** exits non-zero if **no** real manifests are present (only the small `example_pdc_file_manifest.csv`).

## Unfinished / optional next steps

- Refresh **`docs/REPO_AUDIT.md`** tree snapshot after major moves.
- Decide team policy on committing **`reports/benchmark_master/benchmark_results/*.csv`** vs archive-only.
- Consolidate legacy `reports/*.png` exploratory plots into `archive/` if they confuse newcomers (see `CLEANUP_LOG.md`).

## Related checklists

- **`docs/HANDOFF_CHECKLIST.md`** — reproduce / extend tasks.
- **`HANDOFF_SANITY_CHECK.md`** — what was verified in the handoff pass.
