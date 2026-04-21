# Preprocessing pipeline: PSM → gene matrices

## Overview

This stage converts **raw TMT proteomics** (PSM- or peptide-level reporter-ion data) into **gene × sample** abundance matrices in log₂ space. It is intended to run **once per study** (or batch of studies). The **harmonization benchmark** (`scripts/benchmark/run_overnight_v2.sh`) consumes the resulting `gene_matrix.csv` files under `data/results/`.

**Where the code lives today:** most orchestration and the heavy **MSstatsTMT** driver still live under `data/` for historical reasons:

| Artifact / role | Location in this repo |
|-------------------|------------------------|
| Deep pipeline reference (channels, MSstatsTMT stages, QC) | [`data/PIPELINE_README.md`](../../data/PIPELINE_README.md) |
| Main R driver (PSM → protein → gene, per study) | [`data/pdc_psm_to_msstatsTMT_protein_matrix.R`](../../data/pdc_psm_to_msstatsTMT_protein_matrix.R) |
| Batch runner (download + R, space-efficient) | [`data/run_batch_studies.sh`](../../data/run_batch_studies.sh) |
| Per-manifest orchestration | `data/run_pipeline_per_manifest.sh` (see PIPELINE_README) |
| PDC manifest CSVs | `data/manifests/` |
| Downloaded PSMs (large, gitignored) | `data/pdc_psm/{study_id}/` |
| Outputs | `data/results/{study_id}/` (e.g. `msstats_input.tsv`, `ProteinLevelData.csv`, `gene_matrix.csv`) |

This folder (`scripts/preprocessing/`) is the **documentation home** for that pipeline; executable paths above remain under `data/` until a future refactor moves them here without breaking lab workflows.

---

## Input data

### CPTAC (clinical tumors)

**Source:** Proteomics Data Commons (PDC). Benchmark studies used here include **PDC000120** (breast) and **PDC000153** (lung).

Each study provides **PSM-level** tables: one row per spectrum matched to a peptide, with reporter-ion intensities for each TMT channel in the plex.

**Typical workflow**

1. Obtain a **PDC file manifest** CSV for the study and place it under `data/manifests/` (see `data/sample_files_msstats_tmt.csv` for how runs map to files).
2. From the `data/` directory, run the pipeline scripts described in [`data/PIPELINE_README.md`](../../data/PIPELINE_README.md) or use `data/run_batch_studies.sh` if you maintain a local CPTAC mirror (`CPTAC_LOCAL_MIRROR`).

There is **no** committed `download_pdc_manifests.sh` under `scripts/preprocessing/` yet; manifests are usually downloaded from the PDC web UI or API and checked into `data/manifests/` locally.

### CCLE (cell lines)

**Source:** DepMap / Broad peptide-level reporter-ion summaries used by this project live under **`data/ccle_peptide/`** (see that directory and `data/ccle_peptide/ccle_to_msstats_input.py`). CCLE gene matrices for the benchmark are expected at **`data/results/CCLE_corrected/gene_matrix.csv`** after running the same MSstatsTMT-style pipeline with CCLE-specific annotation.

---

## Pipeline stages (conceptual)

### Stage 1: TMT channel annotation

For each study, build the MSstatsTMT annotation table with columns such as **Run**, **Channel**, **Condition**, **Mixture**, **BioReplicate**, **TechRepMixture**. The **Norm** condition marks the **bridge / reference** channel per mixture. MSstatsTMT uses the bridge during summarization; bridge channels are **not** carried as sample columns in the final gene matrix.

### Stage 2: MSstatsTMT protein summarization

`MSstatsTMT::dataProcess` (and related steps in `pdc_psm_to_msstatsTMT_protein_matrix.R`) perform normalization and **protein-level summarization** (e.g. Tukey median polish / model-based summarization depending on configuration). Output includes long **ProteinLevelData** tables.

### Stage 3: Gene-level abundance

Protein summaries are mapped to **HGNC gene symbols** (via `org.Hs.eg.db` in the R pipeline) and collapsed to **`gene_matrix.csv`** (genes × samples, log₂ scale).

### Stage 4: CCLE-specific notes

CAL120 and related duplicate column handling are documented in the benchmark and CCLE utility scripts; see **`scripts/benchmark/process_ccle_annotations_v2.R`** for subtype labels used in the benchmark.

---

## Outputs (per study directory under `data/results/{study_id}/`)

Typical files:

- `msstats_input.tsv` — formatted MSstatsTMT input (includes **Norm** rows used later for **bridge-aware** harmonization).
- `ProteinLevelData.csv` (or `.rds` / variants depending on run) — protein summaries.
- `gene_matrix.csv` — final **gene × sample** matrix consumed by union construction.

---

## Relation to the benchmark

`scripts/run_preprocessing.py` reads **`gene_matrix.csv`** paths declared in `configs/preprocessing/default.yaml` and `configs/tasks/*.yaml`, applies prevalence / variance filters, and writes **`data/processed/union/`** matrices. Bridge statistics for harmonization are extracted from **`msstats_input.tsv`** in **`scripts/benchmark/extract_bridge_summaries.R`** (not from domain medians of `gene_matrix.csv`).

---

## Further reading

- [`data/PIPELINE_README.md`](../../data/PIPELINE_README.md) — authoritative technical detail for PDC → MSstatsTMT.
- [`docs/PREPROCESSING.md`](../../docs/PREPROCESSING.md) — short index + links.
- [`data/README.md`](../../data/README.md) — provenance table and regeneration commands.
