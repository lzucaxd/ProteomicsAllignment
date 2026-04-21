# Contributing

## Repository conventions

- **Exploratory / study-specific** R and Python helpers live under **`data/scripts/`**; the harmonization benchmark entry points live under **`scripts/benchmark/`**.
- **No machine-specific absolute paths** in code under `scripts/`, `src/`, or portable shell under `data/*.sh`. Use:
  - **`PROTEOMICS_ALIGNMENT_ROOT`** for R scripts (see `scripts/benchmark/harmonize_paths.R`).
  - **`CPTAC_LOCAL_MIRROR`** for CPTAC `*.sample.txt` paths listed in `data/sample_files_msstats_tmt.csv` (see [docs/LAB_ONBOARDING.md](docs/LAB_ONBOARDING.md)).
- Prefer **relative paths** in CSV registries when paths are mirror-relative (`PDC000120/...`).

## Pull requests

1. Branch from `main` with a short descriptive name (e.g. `fix/bridge-paths`, `feat/benchmark-metric`).
2. Keep changes focused; avoid drive-by refactors unrelated to the issue.
3. If you add R scripts under `scripts/benchmark/`, resolve the repo root via `harmonize_paths.R` (copy the preamble from an existing script in that folder).
4. CI runs a check for `/Users/` under `scripts/`, `src/`, and `data/` code files.

## Issues

Include: goal, minimal reproduction command, OS, R and Python versions, and whether `CPTAC_LOCAL_MIRROR` / Celligner extras are in use.
