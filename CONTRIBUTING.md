# Contributing

**New to the repo:** read **[`MASTER.md`](MASTER.md)** first (single narrative). **`START_HERE.md`** is a one-click pointer to the same.

## Repository conventions

- **Harmonization benchmark** entry points live under **`scripts/benchmark/`** (see **[`MASTER.md`](MASTER.md)** §5 and [`docs/HOW_TO_RUN_EVERYTHING.md`](docs/HOW_TO_RUN_EVERYTHING.md) for extra command detail).
- **New exploratory / study-specific** helpers → **`scripts/exploratory/`** (preferred). Old one-offs live in **`archive/data_scripts_legacy/`** (same basenames as former `data/scripts/`). **`data/scripts/`** keeps only the two **mapping** Python scripts — see [`data/scripts/README.md`](data/scripts/README.md). Path history: [`docs/NAMING_AND_PATHS.md`](docs/NAMING_AND_PATHS.md).
- **No machine-specific absolute paths** in code under `scripts/`, `src/`, or portable shell under `data/*.sh`. Use:
  - **`PROTEOMICS_ALIGNMENT_ROOT`** for R scripts (see `scripts/benchmark/harmonize_paths.R`).
  - **`CPTAC_LOCAL_MIRROR`** for CPTAC `*.sample.txt` paths listed in `data/sample_files_msstats_tmt.csv` (see [docs/LAB_ONBOARDING.md](docs/LAB_ONBOARDING.md)).
- Prefer **relative paths** in CSV registries when paths are mirror-relative (`PDC000120/...`).

## Pull requests

1. Branch from `main` with a short descriptive name (e.g. `fix/bridge-paths`, `feat/benchmark-metric`).
2. Keep changes focused; avoid drive-by refactors unrelated to the issue.
3. If you add R scripts under `scripts/benchmark/`, resolve the repo root via `harmonize_paths.R` (copy the preamble from an existing script in that folder).
4. Before opening a PR, from repo root: **`python3 scripts/verify_repro_setup.py`** (use **`--require-data`** if you have matrices). CI runs **`verify_repro_setup.py --skip-r`** plus a scan for **`/Users/`** hardcoded paths under `scripts/`, `src/`, and `data/*.sh`.

## Issues

Include: goal, minimal reproduction command, OS, R and Python versions, and whether `CPTAC_LOCAL_MIRROR` / Celligner extras are in use.
