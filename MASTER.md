# Proteomics alignment benchmark — master documentation

**Single entry document.** Read this first; everything else is optional depth. Repo folder on disk: `ProteomicsAllignment` (historical spelling).

---

## Table of contents

1. [What this is](#1-what-this-is)  
2. [Install and verify](#2-install-and-verify)  
3. [Data pipeline (PSM → gene matrix)](#3-data-pipeline-psm--gene_matrixcsv)  
4. [Sample files and MSstats annotation](#4-sample-files-and-msstats-annotation)  
5. [Harmonization benchmark](#5-harmonization-benchmark)  
6. [Configs, paths, and outputs](#6-configs-paths-and-outputs)  
7. [Deep dives (only when you need them)](#7-deep-dives-only-when-you-need-them)  
8. [Repository layout and tidying policy](#8-repository-layout-and-tidying-policy)  
9. [Contributing](#9-contributing)

---

## 1. What this is

A **calibrated benchmark** for **CPTAC (clinical)** vs **CCLE (preclinical)** TMT proteomics on a **shared gene space**, with null/ceiling calibration. Core question: do harmonization methods improve **cross-domain fold-change agreement**, or mainly **PCA geometry**?

**Two statistical layers (do not mix them up):** CPTAC matrices are built with **MSstatsTMT** (TMT design). The overnight benchmark compares methods using **limma** on harmonized **gene** matrices. One page: **`docs/INFERENCE_BASELINES.md`**.

---

## 2. Install and verify

From **repository root**:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -U pip && pip install -r requirements.txt && pip install -e .
Rscript install_r_packages.R
python3 scripts/verify_repro_setup.py
```

After large matrices exist: `python3 scripts/verify_repro_setup.py --require-data`.

Details: **`environment/README.md`**. Mirror variables: **`docs/LAB_ONBOARDING.md`**.

---

## 3. Data pipeline (PSM → gene matrix)

**Goal:** PDC manifest → downloaded `.psm` → **MSstatsTMT** → **`gene_matrix.csv`** per study (under `data/results/{study_id}/`, gitignored when large).

| Step | Where |
|------|--------|
| **Script / shell map** | **`pipeline/psm_to_gene_matrix/README.md`** |
| **Long technical narrative** | **`data/PIPELINE_README.md`** |
| **Run (typical)** | `cd data` then `./run_pipeline_per_manifest.sh` or `./run_batch_studies.sh` (batch often needs `CPTAC_LOCAL_MIRROR`) |
| **Manifest + disk checklist** | **`data/manifests/README.md`**, **`data/manifests/EXPECTED_INPUTS.md`** |

**Not** this pipeline: exploratory code formerly under `data/scripts/` — most one-offs live in **`archive/data_scripts_legacy/`** (see below). The **reproducible** matrix drivers stay under **`data/`** root (`pdc_manifest_downloader.py`, `pdc_psm_to_msstatsTMT_protein_matrix.R`, shell runners).

---

## 4. Sample files and MSstats annotation

**CPTAC:** PDC manifests; **`data/sample_files_msstats_tmt.csv`** lists each study’s **`*.sample.txt`** path; files must exist on disk (mirror). **CCLE:** peptide TSV + sample sheet → Python converter, then the **same** R script with **`--msstats_input_dir`**, *or* use a pre-built **`gene_matrix.csv`** only (benchmark path in YAML).

**How sample design becomes the MSstats annotation table (CPTAC + CCLE):** **`docs/ANNOTATION_FROM_SAMPLES.md`**.

**Benchmark subtype / biospecimen labels (small tables in git):** **`data/annotations/README.md`**, **`data/biospecimen/README.md`**, **`data/ccle/README.md`**. These label samples for **tasks** after matrices exist; they do **not** replace CPTAC **`.sample.txt`** for MSstatsTMT.

**Regenerating CPTAC column ↔ subtype mapping:** `data/scripts/build_PDC000120_subtype_mapping.py` (kept on purpose). **Mixture subset rules (optional science step):** `data/scripts/build_cptac_basal_luminal_mixture_subset.py`.

---

## 5. Harmonization benchmark

**After** CPTAC + CCLE `gene_matrix.csv` exist at the paths in **`configs/preprocessing/default.yaml`**:

```bash
bash scripts/run_benchmark.sh
```

Orchestration and steps: **`scripts/benchmark/README.md`**. Primary tables (when committed): `reports/benchmark_master/benchmark_results/comparison_summary.csv`, `disconnect_scores.csv`. Index: **`reports/final_report/README.md`**.

---

## 6. Configs, paths, and outputs

| Area | Path |
|------|------|
| Preprocessing + study paths | **`configs/preprocessing/default.yaml`**, **`configs/preprocessing/union.yaml`** |
| Tasks | **`configs/tasks/breast_subtype.yaml`**, **`breast_vs_lung.yaml`** |
| Methods | **`configs/methods/`** |
| Union / method matrices | **`data/processed/union/`**, **`data/processed/methods/`** (see overnight README) |

Large downloads and `data/results/` are **gitignored**; see **`.gitignore`** and **`docs/CLEAN_CLONE_REPRODUCIBILITY.md`** for paper / clone policy.

---

## 7. Deep dives (only when you need them)

**Doc policy:** Add new **onboarding** text to **this file** first. Add new **technical depth** only by extending an existing row below or a file already listed—avoid new top-level “start here” guides. Optional notes belong in **`docs/archive/topic_notes/`** or **`archive/`**.

| Topic | File |
|-------|------|
| Sample → MSstats annotation (CPTAC + CCLE) | **`docs/ANNOTATION_FROM_SAMPLES.md`** |
| MSstatsTMT vs limma | **`docs/INFERENCE_BASELINES.md`** |
| Harmonization method definitions | **`docs/METHODS.md`** |
| Paper-length narrative | **`docs/END_TO_END_TECHNICAL_REPORT.md`** |
| Overnight shell steps | **`scripts/benchmark/README.md`** |
| YAML / new tasks / extending | **`docs/HANDOFF_CHECKLIST.md`**, **`docs/archive/topic_notes/config_system_overview.md`** |
| Why folders look scattered | **`docs/NAMING_AND_PATHS.md`** |
| Slides + v2 checklist | **`docs/BENCHMARK_V2_AND_PRESENTATION.md`** |
| Lab handoff + caveats (short) | **`HANDOFF.md`** |
| Inventory / long prose | **`PROJECT_REPORT.md`** |
| What moved during cleanup | **`CLEANUP_LOG.md`** |

**Archived topic notes** (benchmark philosophy, structure metrics, etc.): **`docs/archive/topic_notes/README.md`**.

---

## 8. Repository layout and tidying policy

- **`src/harmonize/`** — Python package (preprocessing, methods registry).  
- **`data/`** — PSM download, MSstatsTMT R, manifests; **two** Python scripts remain under **`data/scripts/`** (subtype mapping builders); older exploratory R/Python moved to **`archive/data_scripts_legacy/`** (paths in old slide markdown may still say `data/scripts/` — use **`archive/data_scripts_legacy/`** + same filename).  
- **`scripts/benchmark/`** — overnight benchmark.  
- **`archive/`** — legacy outputs and migrated exploratory scripts.  
- **`reports/benchmark_master/`** — numerical outputs (name says “reports”; treat as results).  

---

## 9. Contributing

See **`CONTRIBUTING.md`**. New one-off scripts: **`scripts/exploratory/`**, not `data/scripts/`.

**CI / checks:** `python3 scripts/verify_repro_setup.py` before PRs.

---

## Legacy entry files (redirects)

**`START_HERE.md`** and the top of **`README.md`** point here so GitHub clones still work. Prefer opening **this file** (`MASTER.md`) for the full picture in one place.
