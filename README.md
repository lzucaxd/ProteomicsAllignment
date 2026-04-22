# Proteomics alignment benchmark

**CPTAC ↔ CCLE TMT proteomics** — a calibrated benchmark on a **shared gene space** with explicit **null** and **ceiling** calibration. The repo evaluates whether harmonization improves **cross-domain fold-change agreement** or mainly **PCA geometry**.

> **Folder name:** the directory is spelled `ProteomicsAllignment` on disk (historical).

---

## Quick links

| Goal | Document |
|------|----------|
| **Run everything (clone → data → benchmark)** | [`docs/HOW_TO_RUN_EVERYTHING.md`](docs/HOW_TO_RUN_EVERYTHING.md) |
| **Manifest → PSM → gene matrix (front door)** | [`pipeline/psm_to_gene_matrix/README.md`](pipeline/psm_to_gene_matrix/README.md) |
| **MSstatsTMT (native TMT) vs limma (benchmark)** | [`docs/INFERENCE_BASELINES.md`](docs/INFERENCE_BASELINES.md) |
| **Lab handoff (new teammate)** | [`HANDOFF.md`](HANDOFF.md) |
| **Environment (Python + R)** | [`environment/README.md`](environment/README.md) |
| **Layout & clutter policy** | [`REPO_AUDIT.md`](REPO_AUDIT.md), [`REPO_LAYOUT_PLAN.md`](REPO_LAYOUT_PLAN.md) |
| **Final tables / figures index** | [`reports/final_report/README.md`](reports/final_report/README.md) |
| **What changed during cleanup** | [`CLEANUP_LOG.md`](CLEANUP_LOG.md) |

---

## Repository layout

```
configs/              YAML: preprocessing, tasks, methods, benchmark toggles
src/harmonize/        Python package (preprocessing, benchmark helpers, methods registry)
pipeline/             Doc-only index: PSM→matrix front door (see pipeline/psm_to_gene_matrix/)
data/                 CPTAC/CCLE drivers: PSM download, MSstatsTMT R, manifests (see data/README.md)
scripts/
  benchmark/          run_overnight_v2.sh (implementation), metrics, calibration R
  exploratory/        Preferred home for new one-offs (not data/scripts/)
  preprocessing/    Documentation index (executables still under data/)
  run_benchmark.sh    Stable entry → overnight v2
  run_diagnostics.sh  Preflight + structure subset
  run_methods.sh      Python method smoke test
reports/
  benchmark_master/   CSV outputs, diagnostics, final_tables/
  final_report/       Index README + figure manifest (no duplicate large assets)
docs/                 Technical reports, methods, reproducibility guides
notebooks/exploratory/ Non-pipeline notebooks
archive/              Legacy duplicates (see archive/README.md)
environment/          Setup instructions
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

## Run the main pipeline

1. **Inputs:** PDC manifests + CPTAC `.sample.txt` paths + CCLE matrices — **`data/manifests/README.md`**, **`data/manifests/EXPECTED_INPUTS.md`**, **`data/PIPELINE_README.md`**.
2. **PSM → `gene_matrix.csv`:** from `data/`: `./run_pipeline_per_manifest.sh` (or `./run_batch_studies.sh` with `CPTAC_LOCAL_MIRROR`).
3. **Benchmark:** from repo root: **`bash scripts/run_benchmark.sh`** (hours; slow steps documented in `scripts/benchmark/README.md`).

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
