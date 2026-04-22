# PDC manifest → PSM → gene matrix (reproducible CPTAC/CCLE driver)

This is the **lab handoff “front door”** for the heavy preprocessing stage. Implementations remain under **`data/`** (historical paths used by shell, R, and Python across the repo).

## Sample files → MSstatsTMT annotation (read this first)

**One canonical explanation (CPTAC `*.sample.txt` + CCLE Sheet2 / converter → same R script):** **[`../../docs/ANNOTATION_FROM_SAMPLES.md`](../../docs/ANNOTATION_FROM_SAMPLES.md)**.

**Checklist of paths on disk:** [`../../data/manifests/EXPECTED_INPUTS.md`](../../data/manifests/EXPECTED_INPUTS.md). **Mirror / env:** [`../../docs/LAB_ONBOARDING.md`](../../docs/LAB_ONBOARDING.md).

Do not confuse **MSstats design** (above) with **benchmark subtype tables** in **`data/annotations/`** / **`data/biospecimen/`** — those label samples for tasks **after** `gene_matrix.csv` exist.

---

## What to run (orchestration)

| Step | Script / entry | Directory |
|------|----------------|-----------|
| All manifests in `data/manifests/` | **`data/run_pipeline_per_manifest.sh`** | run from **`data/`** |
| Batch by study + `CPTAC_LOCAL_MIRROR` | **`data/run_batch_studies.sh`** | `data/` |
| Single study shortcuts | **`data/run_study_PDC000127.sh`**, **`data/run_study_PDC000153.sh`** | `data/` |

## Core programs (same folder: `data/`)

| Role | File |
|------|------|
| Download `.psm` from PDC manifest CSV | **`data/pdc_manifest_downloader.py`** |
| PSM parsing → MSstatsTMT → protein → **gene** matrix | **`data/pdc_psm_to_msstatsTMT_protein_matrix.R`** |
| Manifest + export policy | **`data/manifests/README.md`**, **`data/manifests/EXPECTED_INPUTS.md`** |
| Full narrative + data flow diagram | **`data/PIPELINE_README.md`** |

## Outputs (typical; gitignored when large)

- `data/pdc_psm/{study_id}/…/*.psm`
- `data/results/{study_id}/gene_matrix.csv`, `msstats_input.tsv`, etc.

## What this is **not**

- **`data/scripts/`** — exploratory subtype / v1 / one-off DA and QC (see **`data/scripts/README.md`**). Those files are **not** required to reproduce the main manifest → matrix → benchmark path.
- **`scripts/benchmark/`** — starts **after** `gene_matrix.csv` exist (union, methods, limma evaluation).

## MSstatsTMT vs limma (where inference happens)

- **Inside this stage:** CPTAC matrices come from **MSstatsTMT** summarization in **`pdc_psm_to_msstatsTMT_protein_matrix.R`** (TMT-aware).
- **Downstream benchmark:** per-domain **limma** on harmonized (or raw) **gene** matrices — see **`docs/INFERENCE_BASELINES.md`**.
