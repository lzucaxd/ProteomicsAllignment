# PDC manifest → PSM → gene matrix (reproducible CPTAC/CCLE driver)

This is the **lab handoff “front door”** for the heavy preprocessing stage. Implementations remain under **`data/`** (historical paths used by shell, R, and Python across the repo).

## Required inputs: CPTAC sample files, annotations, CCLE

Do not confuse **(A) MSstatsTMT design files** (needed to run the matrix pipeline) with **(B) benchmark subtype / biospecimen tables** (small CSVs in git for Luminal–Basal tasks **after** `gene_matrix.csv` exist).

### CPTAC (clinical TMT)

| Input | Where it lives | Role |
|--------|----------------|------|
| **PDC manifest CSV** | `data/manifests/` (local export; URLs expire) | Lists `.psm` URLs for **`pdc_manifest_downloader.py`**. |
| **Study registry** | **`data/sample_files_msstats_tmt.csv`** | One row per `study_id`: mirror-relative **`path`** to that study’s **`*.sample.txt`**, TMT format, `reference_channel`, optional **`annotation_path`**, `use_for_msstats_tmt`. |
| **`*.sample.txt` on disk** | Resolved via **`path`** + often **`CPTAC_LOCAL_MIRROR`** (parent of `PDC000120/`, …) | **MSstatsTMT channel / mixture / bridge design**; R uses `--sample_txt` to build or correct **`annotation_filled_corrected.csv`** under `data/results/{study_id}/`. |
| **Per-study annotation CSVs** | Usually created under **`data/results/{study_id}/`** (`annotation_filled*.csv`) | Produced or updated by the R driver from the sample file; not the same as `data/annotations/` (benchmark mapping). |

Details: **`data/manifests/EXPECTED_INPUTS.md`** §1–2, **`data/PIPELINE_README.md`** (registry + `--sample_txt`), **`docs/LAB_ONBOARDING.md`** (mirror env vars).

### CCLE (preclinical path through the same R driver)

| Input | Where | Role |
|-------|--------|------|
| **Peptide TSV + sample table** | Typically **`data/ccle_peptide/`** (see READMEs there) | Convert to **`msstats_input.tsv`** + annotation via **`ccle_to_msstats_input.py`**. |
| **Run R with `--msstats_input_dir`** | e.g. `results/CCLE` | Same **`pdc_psm_to_msstatsTMT_protein_matrix.R`** path as CPTAC, skipping PSM parse when input is already MSstats-shaped. |

Benchmark-only CCLE **gene matrix** path (often already summarized): **`configs/preprocessing/default.yaml`** → e.g. `data/results/CCLE_corrected/gene_matrix.csv` — that is **downstream** of peptide conversion unless you only ingest a pre-built matrix.

### Benchmark tables in git (not a substitute for `.sample.txt`)

**`data/annotations/`**, **`data/biospecimen/`**, **`data/ccle/ccle_breast_subtype_annotations_v2.csv`** — used for **task labels** (subtype, tissue contrasts) and preprocessing metadata **after** matrices exist. They do **not** replace CPTAC **`sample_files_msstats_tmt.csv`** + **`.sample.txt`** for running MSstatsTMT on PSMs.

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
