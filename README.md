# Proteomics alignment benchmark

**CPTAC ↔ CCLE TMT proteomics** — a calibrated benchmark on a **shared gene space** with explicit **null** and **ceiling** calibration. The repo evaluates whether harmonization improves **cross-domain fold-change agreement** or mainly **PCA geometry**.

> **Folder name:** the directory is spelled `ProteomicsAllignment` on disk (historical).

---

## Cloning from git — read this first

**Default entry (data pipeline + run order):** **[`START_HERE.md`](START_HERE.md)** — optimized for GitHub / first clone; leads with **how to run the data pipeline** (`PSM` → **`gene_matrix.csv`**) and where the scripts live.

**How sample design becomes MSstatsTMT annotation (CPTAC + CCLE, one page):** **[`docs/ANNOTATION_FROM_SAMPLES.md`](docs/ANNOTATION_FROM_SAMPLES.md)**.

**Single end-to-end runbook:** [`docs/HOW_TO_RUN_EVERYTHING.md`](docs/HOW_TO_RUN_EVERYTHING.md).

**Data pipeline file map (shell + R under `data/`):** [`pipeline/psm_to_gene_matrix/README.md`](pipeline/psm_to_gene_matrix/README.md).

**Keep the doc set small:** `START_HERE` → `ANNOTATION_FROM_SAMPLES` (if needed) → `HOW_TO_RUN_EVERYTHING` → `HANDOFF`. Other files are optional depth.

---

## Quick links

| Goal | Document |
|------|----------|
| **First clone (start here)** | **[`START_HERE.md`](START_HERE.md)** |
| **Sample files → MSstats annotation (CPTAC + CCLE)** | **[`docs/ANNOTATION_FROM_SAMPLES.md`](docs/ANNOTATION_FROM_SAMPLES.md)** |
| **Run everything (clone → data → benchmark)** | [`docs/HOW_TO_RUN_EVERYTHING.md`](docs/HOW_TO_RUN_EVERYTHING.md) |
| **Manifest → PSM → gene matrix (front door)** | [`pipeline/psm_to_gene_matrix/README.md`](pipeline/psm_to_gene_matrix/README.md) |
| **MSstatsTMT (native TMT) vs limma (benchmark)** | [`docs/INFERENCE_BASELINES.md`](docs/INFERENCE_BASELINES.md) |
| **Lab handoff (new teammate)** | [`HANDOFF.md`](HANDOFF.md) |
| **Environment (Python + R)** | [`environment/README.md`](environment/README.md) |
| **Why paths look scattered + naming rules** | [`docs/NAMING_AND_PATHS.md`](docs/NAMING_AND_PATHS.md) |
| **Layout & clutter policy** | [`REPO_AUDIT.md`](REPO_AUDIT.md), [`REPO_LAYOUT_PLAN.md`](REPO_LAYOUT_PLAN.md) |
| **Final tables / figures index** | [`reports/final_report/README.md`](reports/final_report/README.md) |
| **What changed during cleanup** | [`CLEANUP_LOG.md`](CLEANUP_LOG.md) |
| **Subtype / tissue annotations (tracked tables)** | [`data/annotations/README.md`](data/annotations/README.md) |

---

## Repository layout

```
configs/               YAML: preprocessing, tasks, methods, benchmark toggles
src/harmonize/         Python package (preprocessing, benchmark helpers, methods registry)
pipeline/              Doc-only: PSM→matrix front door (executables still under data/)
data/                  Stage-1 CPTAC/CCLE: PSM download, MSstatsTMT R, manifests (see data/README.md)
  scripts/             Legacy exploratory R/py (not the manifest→matrix driver)
scripts/
  benchmark/           Stage-2: run_overnight_v2.sh, limma, metrics, calibration R
  exploratory/         Preferred home for new one-offs
  preprocessing/       Doc index only; drivers live under data/
  run_benchmark.sh     Stable entry → overnight v2
  run_diagnostics.sh   Preflight + structure subset
  run_methods.sh       Python method smoke test
reports/
  benchmark_master/    Numerical outputs live here (name says “reports”; treat as results/)
  final_report/        Curated index + figure manifest
docs/                  Technical reports; see START_HERE.md then HOW_TO_RUN_EVERYTHING.md
notebooks/exploratory/ Non-pipeline notebooks
archive/               Legacy duplicates (see archive/README.md)
environment/           Setup instructions
```

---

## Setup (short)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -U pip && pip install -r requirements.txt && pip install -e .
Rscript install_r_packages.R
python3 scripts/verify_repro_setup.py
```

Details and optional Celligner: **`environment/README.md`**.

---

## Data pipeline and benchmark (how to run)

**Same content, more detail:** **[`START_HERE.md`](START_HERE.md)** and **[`pipeline/psm_to_gene_matrix/README.md`](pipeline/psm_to_gene_matrix/README.md)**.

1. **Inputs:** PDC **manifests**; CPTAC **`data/sample_files_msstats_tmt.csv`** + on-disk **`*.sample.txt`** (design for MSstatsTMT); CCLE either **peptide pipeline** under **`data/ccle_peptide/`** or a **pre-built `gene_matrix.csv`** in **`configs/preprocessing/default.yaml`**. Benchmark **subtype/biospecimen** tables in **`data/annotations/`**, **`data/biospecimen/`**, **`data/ccle/`** are separate (see **`data/manifests/EXPECTED_INPUTS.md`**, **`data/PIPELINE_README.md`**).
2. **Data pipeline — PSM → `gene_matrix.csv`:** run from **`data/`**: `./run_pipeline_per_manifest.sh` or `./run_batch_studies.sh` (often needs **`CPTAC_LOCAL_MIRROR`**; see **`docs/LAB_ONBOARDING.md`**). Driver + input table: **`pipeline/psm_to_gene_matrix/README.md`**.
3. **Benchmark (after matrices exist):** from repo root: **`bash scripts/run_benchmark.sh`** (hours; steps in **`scripts/benchmark/README.md`**).

**Primary outputs:** `reports/benchmark_master/benchmark_results/comparison_summary.csv`, `disconnect_scores.csv`. **Index:** `reports/final_report/README.md`.

---

## Data policy

Large and expiring artifacts are **gitignored** (PSM downloads, `gene_matrix.csv` trees, dated PDC manifests, presentation bundles). The repo stays **small and cloneable**; you attach data locally or via shared storage. See **`.gitignore`** and **`docs/CLEAN_CLONE_REPRODUCIBILITY.md`**.

---

## Methods compared

| ID | Role |
|----|------|
| `raw` | Baseline |
| `bridge_shift` | Bridge-channel offset |
| `bridge_scale` | Bridge offset + MAD scale |
| `celligner` | cPCA + MNN |

**Definitions:** [`docs/METHODS.md`](docs/METHODS.md).

---

## Project status

**Active research codebase** — benchmark v2 shell path is primary; some legacy Python drivers remain for debugging. **Sanity checklist after refactors:** [`HANDOFF_SANITY_CHECK.md`](HANDOFF_SANITY_CHECK.md).

---

## Citation

Zamfira, L.-A. (2025–2026). *A calibrated benchmark for cross-dataset harmonization of clinical and preclinical cancer proteomics.* Vitek Lab, Northeastern University.

## License

[`LICENSE`](LICENSE) (MIT).
