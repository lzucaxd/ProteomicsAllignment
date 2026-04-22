# Repository audit (lab handoff)

**Purpose:** High-level inventory for someone cloning the repo. For a machine-generated tree snapshot, see [`docs/REPO_AUDIT.md`](docs/REPO_AUDIT.md) (refresh that file after large moves).

**Last reviewed:** 2026-04-21 (conservative handoff pass; no wholesale directory reshuffle).

---

## Top-level layout (what each area is for)

| Path | Role |
|------|------|
| **`README.md`** | Project entry: goals, layout, links to run guides |
| **`HANDOFF.md`** | Practical handoff for a new lab member |
| **`configs/`** | YAML: preprocessing sources, tasks, methods, benchmark toggles |
| **`src/harmonize/`** | Python package: preprocessing helpers, benchmark task defs, method registry |
| **`data/`** | CPTAC/CCLE pipeline home: PSM download, MSstatsTMT R driver, manifests, results (mostly **gitignored**) |
| **`scripts/benchmark/`** | **Primary benchmark orchestration** (`run_overnight_v2.sh`), metrics, calibration R |
| **`scripts/preprocessing/`** | Documentation index; executables still under `data/` (historical) |
| **`scripts/methods/`**, **`scripts/presentation/`** | Method drivers and figure/slide pipelines |
| **`reports/benchmark_master/`** | Benchmark CSVs, diagnostics, meeting figures, **`final_tables/`** for paper-style bundles |
| **`docs/`** | Design notes, technical reports, reproducibility guides |
| **`archive/`** | Legacy duplicates moved aside (see `archive/README.md`) |
| **`models/`** | Optional Celligner checkout (often local; partial gitignore) |
| **`presentation_materials/`** | Generated slide/report assets (**gitignored**; regenerate via scripts) |

---

## What is safe to **ignore** (git / local)

- **`.venv/`, caches** — in `.gitignore`
- **`data/pdc_psm/`, `data/results/`, large downloads** — gitignored; reproduce via `data/PIPELINE_README.md`
- **Dated `PDC_file_manifest_*.csv`** and **`PDC*_pdc_file_manifest.csv`** — gitignored; URLs expire; see `data/manifests/README.md`
- **`presentation_materials/`** — gitignored; regenerate
- **`reports/benchmark_master/logs/*.log`** — gitignored

## What is safe to **archive** (not delete)

- Duplicate or superseded figure trees → **`archive/`** with a line in `archive/README.md`
- Exploratory notebooks → **`notebooks/exploratory/`** (see `notebooks/exploratory/README.md`)

## What must stay **prominent**

- **`docs/HOW_TO_RUN_EVERYTHING.md`** — single run guide
- **`scripts/benchmark/run_overnight_v2.sh`** — canonical full benchmark
- **`configs/preprocessing/default.yaml`**, **`configs/tasks/*.yaml`**
- **`data/PIPELINE_README.md`**, **`data/manifests/README.md`**, **`data/sample_files_msstats_tmt.csv`** (registry)
- **`reports/benchmark_master/benchmark_results/comparison_summary.csv`** and **`disconnect_scores.csv`** (when committed or archived for a paper tag)
- **`reports/benchmark_master/final_tables/`** — consolidated numbers / LLM-oriented bundle (see `reports/final_report/` index)

---

## Duplicate or overlapping entry points (known)

| Situation | Guidance |
|-----------|------------|
| **`run_overnight_v2.sh`** vs `scripts/run_all.py` | Use **overnight v2** for full runs; `run_all.py` is legacy / debugging |
| **`scripts/run_methods.py`** vs **`regenerate_methods_union.py`** | Overnight benchmark uses **R + union** path (`regenerate_methods_union.py`); `run_methods.py` loads the Python harmonize registry (lighter smoke test) |
| **`configs/benchmark/default.yaml`** | Canonical benchmark YAML; no separate `benchmark.yaml` required |

---

## Obsolete / exploratory (policy)

- **`data/gene_matrix_exploration.ipynb`** — moved to **`notebooks/exploratory/`** (not referenced by pipeline)
- **Notebooks under `models/celligner-master/`** — vendor; do not move

---

## Large / external data

Nothing in **`data/results/`** or **`data/pdc_psm/`** is assumed present on clone. See **`data/manifests/EXPECTED_INPUTS.md`** and **`environment/README.md`**.
