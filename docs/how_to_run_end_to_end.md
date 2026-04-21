# How to run end-to-end

**Canonical benchmark:** from the repo root run **`scripts/benchmark/run_overnight_v2.sh`**. Steps, outputs, and slide assets are documented in **[BENCHMARK_V2_AND_PRESENTATION.md](BENCHMARK_V2_AND_PRESENTATION.md)** and the full narrative in **[END_TO_END_TECHNICAL_REPORT.md](END_TO_END_TECHNICAL_REPORT.md)**.

**PDC → gene matrix:** run from **`data/`** — see **[data/PIPELINE_README.md](../data/PIPELINE_README.md)** and [LAB_ONBOARDING.md](LAB_ONBOARDING.md) (`CPTAC_LOCAL_MIRROR`).

**Legacy (optional):** `scripts/run_all.py` chains older Python entry points (`run_preprocessing.py`, `run_native_baselines.py`, `run_methods.py`, `run_benchmark.py`, `run_meeting_exports.py`). Prefer the overnight shell script unless you are debugging a single Python step.
