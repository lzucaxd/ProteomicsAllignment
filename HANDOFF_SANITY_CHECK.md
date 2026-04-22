# Handoff sanity check

**Date:** 2026-04-21 (documentation + wrapper pass; no scientific logic changes).

## Verified (static / repo consistency)

- [x] **Wrapper scripts exist** and point to real targets: `scripts/run_benchmark.sh` → `scripts/benchmark/run_overnight_v2.sh`; `scripts/run_methods.sh` → `scripts/run_methods.py`; `scripts/run_diagnostics.sh` calls `preflight_diagnostics.R`, `compute_intersection_masks.R`, `run_structure_batch.py` with **`--meta-dir`** aligned to `data/processed/union` by default (`UNION_DIR` / `META_DIR`).
- [x] **Notebook move** `data/gene_matrix_exploration.ipynb` → `notebooks/exploratory/` — **no** code references found (`grep gene_matrix_exploration`).
- [x] **Manifest policy** documented; dated manifests removed from version control in commit `ea69be5` (see `CLEANUP_LOG.md`).
- [x] **New handoff docs** cross-link to existing canonical guides (`HOW_TO_RUN_EVERYTHING.md`, `PIPELINE_README.md`, `environment/README.md`).

## Needs manual checking (requires your machine / data)

- [ ] **`bash scripts/run_benchmark.sh`** end-to-end after fresh manifests + `gene_matrix.csv` (hours; permutations + ceilings slow).
- [ ] **`./scripts/run_diagnostics.sh all`** after union matrices exist under `data/processed/union/`.
- [ ] **R package versions** vs `install_r_packages.R` on a clean R install.
- [ ] **Presentation figures** regeneration if you rely on `presentation_materials/` (gitignored on purpose).

## Could not test without large external data

- Full CPTAC PSM download + MSstatsTMT protein→gene pipeline.
- Celligner optional path and GPU/CPU runtime.
- Full permutation null and split-half ceiling steps (Step 7–8 of overnight).

## Path / reference audit (sampling)

- Grep-based checks for **`gene_matrix_exploration`** after move: clean.
- **Report-facing paths** in new `reports/final_report/*.md` use relative `../benchmark_master/...` — verify after any future directory rename.

## If something breaks

1. Compare your command to **`docs/HOW_TO_RUN_EVERYTHING.md`**.
2. Read the tail of **`reports/benchmark_master/logs/overnight_v2_*.log`**.
3. Run **`python3 scripts/verify_repro_setup.py --require-data`**.
4. File an issue with OS, Python/R versions, and whether `CPTAC_LOCAL_MIRROR` / Celligner extras are in use (`CONTRIBUTING.md`).
