# Contributing

**New to the repo (clone / GitHub):** read **[`START_HERE.md`](START_HERE.md)** first — data pipeline run order and benchmark entry points.

## Repository conventions

- **Harmonization benchmark** entry points live under **`scripts/benchmark/`** (see [`docs/HOW_TO_RUN_EVERYTHING.md`](docs/HOW_TO_RUN_EVERYTHING.md)).
- **New exploratory / study-specific** helpers → **`scripts/exploratory/`** (preferred). Legacy and slide-linked code remains under **`data/scripts/`** — see [`docs/NAMING_AND_PATHS.md`](docs/NAMING_AND_PATHS.md) for the full map (why paths look scattered and what not to duplicate).
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
