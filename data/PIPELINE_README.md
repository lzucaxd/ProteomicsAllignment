# PDC TMT → MSstatsTMT → Gene Matrix Pipeline

End-to-end documentation for going from PDC file manifests to sample × gene matrices and QC plots.

The pipeline supports **arbitrary TMT plex sizes** (TMT6, TMT10, TMT11, TMT16/TMTpro, TMT18) for **both** CPTAC and CCLE (or other pre-built input):

- **CPTAC path**: TMT channel columns in PSM files are auto-detected from the header (e.g. `TMT10-126`, `TMT11-131C`); `*.sample.txt` channel columns are detected by pattern (e.g. `126`, `127N`, `131C`, `134N`). No hardcoded channel list.
- **CCLE / pre-built path**: Reporter-ion columns in the peptide TSV are auto-detected (`rq_*_sn`); sample table channels are normalized to the same canonical labels.

The annotation schema is fixed (Run, Channel, Condition, BioReplicate, Mixture, TechRepMixture); **exactly one bridge (Condition = Norm) per mixture** is required.

---

## Purpose of this pipeline

This pipeline converts CPTAC PSM-level proteomics data from the Proteomic Data Commons (PDC) into **sample × gene matrices** using MSstatsTMT.

The goal is to produce **clean, standardized gene abundance matrices** that can be used for downstream machine learning models and cross-study analysis.

**This pipeline does not perform machine learning preprocessing.** Filtering, imputation, and ML-specific preprocessing are not performed here and will happen later during model training.

---

## Data flow

The pipeline transforms data through the following stages and **stops at the gene matrix stage**:

```
PDC Manifest
   ↓
PSM Files (.psm)
   ↓
PSM Parsing
   ↓
MSstatsTMT Input Formatting
   ↓
MSstatsTMT Protein Summarization
   ↓
Protein Abundance Matrix
   ↓
Protein → Gene Mapping
   ↓
Sample × Gene Matrix
```

---

## What this pipeline does

The pipeline performs the following processing steps:

- **Download of PSM files** from PDC manifests into `pdc_psm/{study_id}/`.
- **Parsing of PSM files** into long-format peptide data (Run × Channel × PSM).
- **Construction of MSstatsTMT input tables** with Run, Channel, Condition, BioReplicate, Mixture, Fraction, and intensity columns.
- **Reference channel normalization** for TMT experiments (bridge channel used as reference).
- **Peptide-to-protein summarization** using MSstatsTMT (median polish / model-based summarization to estimate protein abundances).
- **Removal of bridge/reference channels** from the final matrix.
- **Mapping of proteins to gene symbols** using `org.Hs.eg.db` (or `org.Mm.eg.db` for mouse).
- **Generation of a sample × gene abundance matrix** per study.

Protein abundances are estimated using MSstatsTMT summarization methods (e.g. median polish, model-based summarization). This pipeline performs **minimal processing only** and stops at the gene matrix stage.

---

## What this pipeline does NOT do

The pipeline **intentionally does not perform** additional preprocessing steps. Specifically, it does **not** perform:

- Missing value imputation
- Protein or gene filtering based on missingness
- Feature selection
- Cross-study normalization
- Batch correction
- Dimensionality reduction
- Scaling or standardization

These steps will be performed later during model training and downstream analysis. They are kept separate from this pipeline so that raw gene matrices remain consistent and reusable.

---

## MSstatsTMT normalization (clarification)

MSstatsTMT performs:

- **TMT reference channel normalization** (using the bridge/reference channel to normalize within-plex intensities).
- **Peptide-level summarization to proteins** (aggregating peptide abundances to protein level).

This is **required for TMT proteomics data** to obtain comparable protein abundances across channels and runs. It is **not** equivalent to dataset-wide normalization for machine learning (e.g. cross-study scaling or batch correction), which is done later in downstream pipelines.

---

## Gene mapping

Protein accessions are mapped to gene symbols using the Bioconductor package **org.Hs.eg.db** (human) or **org.Mm.eg.db** (mouse).

If **multiple proteins map to the same gene symbol**, they are aggregated to generate gene-level abundance values in the final sample × gene matrix.

---

## Overview

| Step | What it does |
|------|----------------|
| 1. Manifests | One CSV per study (from PDC portal). Duplicate manifests removed. |
| 2. Download | PSM files (`.psm`) downloaded per manifest into `pdc_psm/{study_id}/`. |
| 3. Sample-file check | Each study must have a row in `sample_files_msstats_tmt.csv` (path to `.sample.txt`). |
| 4. R pipeline | PSM → MSstatsTMT (annotation, summarization) → gene matrix and per-study annotation. |
| 5. QC plots | Optional: publication-style QC and profile plots per study. |

Each study gets its **own** outputs under `results/{study_id}/` and its **own** annotation files there.

---

## Directory layout

```
.
├── manifests/                          # PDC file manifest CSVs (one per study)
│   └── PDC_file_manifest_*.csv
├── sample_files_msstats_tmt.csv        # Registry: study_id → path to .sample.txt, format, reference_channel
├── pdc_psm/                            # Downloaded PSM files (created by pipeline)
│   └── {study_id}/
│       └── Peptide_Spectral_Matches/
│           └── *.psm
├── results/                            # Pipeline outputs (created by pipeline)
│   └── {study_id}/
│       ├── annotation_filled.csv
│       ├── annotation_filled_corrected.csv   # From sample.txt when using --sample_txt
│       ├── annotation_audit.txt
│       ├── msstats_input.tsv
│       ├── protein_summary.tsv
│       ├── gene_matrix.csv
│       ├── qc_summary.txt
│       └── plots/                     # MSstatsTMT QCPlot (and optional ProfilePlots) for lab debugging
├── run_pipeline_per_manifest.sh        # Main runner: download + R per manifest
├── pdc_manifest_downloader.py          # Download files from manifest
├── pdc_psm_to_msstatsTMT_protein_matrix.R   # PSM → gene matrix (MSstatsTMT)
├── msstatsTMT_qc_plots.R               # Standalone: MSstatsTMT built-in QC/Profile plots only
├── check_manifests.py                  # List studies, find duplicate manifests
└── check_studies_sample_file.py        # Check studies have entry in sample_files_msstats_tmt.csv
```

---

## Prerequisites

- **Python 3** with `requests` (e.g. `pip install requests` or use `.venv`).
- **R** with:
  - `MSstatsTMT` (BiocManager)
  - `org.Hs.eg.db` (human) or `org.Mm.eg.db` (mouse) for gene mapping
  - `data.table` (for optional standalone QC script)

---

## Inputs

### 1. Manifests (`manifests/`)

- PDC file manifest CSVs downloaded from the [PDC portal](https://pdc.cancer.gov/).
- Each file should contain **one study** (one “PDC Study ID” value).
- Required columns include: **File Name**, **PDC Study ID**, **Data Category**, **File Download Link**.

### 2. Sample file registry (`sample_files_msstats_tmt.csv`)

- CSV with columns: `study_id`, `path`, `file_name`, `format`, `reference_channel`, `n_channels`, `delimiter`, `use_for_msstats_tmt`, `annotation_path`, `notes`.
- **path**: CPTAC `.sample.txt` for that study. In git this is usually **mirror-relative**, e.g. `PDC000120/CPTAC2_Breast_....sample.txt` (study folder + file name). Resolution order:
  1. Absolute path, if the file exists.
  2. Relative to the `data/` directory (`data/<path>`).
  3. If **`CPTAC_LOCAL_MIRROR`** is set: `<CPTAC_LOCAL_MIRROR>/<path>` (typical: export mirror to the parent of `PDC000120/`, `PDC000153/`, …).
- **format**: e.g. `TMT10`, `TMT11`; **reference_channel**: e.g. `131`, `126C` (bridge channel).
- Every study you run **must** have a row here (or the pipeline will warn and may need `--reference_channel` / manual annotation). Lab setup: [docs/LAB_ONBOARDING.md](../docs/LAB_ONBOARDING.md).

---

## Running the pipeline

### Option A: One command (all manifests)

From the pipeline root (where `run_pipeline_per_manifest.sh` lives):

```bash
./run_pipeline_per_manifest.sh
```

This script:

1. Loops over every `manifests/PDC_file_manifest_*.csv`.
2. Extracts **PDC Study ID** from each.
3. **Downloads** PSM files (`.psm` only) into `pdc_psm/{study_id}/`.
4. **Checks** that the study appears in `sample_files_msstats_tmt.csv` (warns if not).
5. Runs **R pipeline** with `--psm_dir pdc_psm/{study_id}` and `--outdir results/{study_id}`.

Each study’s results and annotation files are written under `results/{study_id}/`.

#### Efficient run (low disk)

To avoid filling your disk, run with **cleanup after each study**. Peak usage stays at roughly one study's size (~5–10 GB) instead of all 12.

```bash
# After each study: keep gene_matrix + annotation, remove heavy intermediates (parsed_psm_long, msstats_input)
./run_pipeline_per_manifest.sh --cleanup-after

# Also remove raw PSM download for each study after its pipeline succeeds (maximum space saving)
./run_pipeline_per_manifest.sh --cleanup-after --delete-psm
```

With `--cleanup-after`, once a study's `gene_matrix.csv` is written, the script deletes `parsed_psm_long.tsv` and `msstats_input.tsv` for that study. With `--delete-psm` it also removes `pdc_psm/{study_id}/`. You keep `gene_matrix.csv`, annotation, `protein_summary.tsv`, and QC files. Re-running a study later requires re-downloading and re-running the R pipeline if you used `--delete-psm`. **Note:** MSstatsTMT QC plots are written to `results/{study_id}/plots/` during the pipeline. The standalone `msstatsTMT_qc_plots.R` needs `msstats_input.tsv` and `protein_summary.tsv`; run it before using `--cleanup-after` if you want to regenerate those plots later.

### Option B: Single study (manual)

```bash
# 1. Download (one manifest)
python3 pdc_manifest_downloader.py \
  --manifest manifests/PDC_file_manifest_03142026_145054.csv \
  --outdir pdc_psm \
  --include-category "Peptide Spectral Matches" \
  --ext .psm

# 2. R pipeline (use --sample_txt if you have the .sample.txt path from sample_files_msstats_tmt.csv)
Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R \
  --psm_dir pdc_psm/PDC000614 \
  --outdir results/PDC000614 \
  --sample_txt /path/to/PDC000614/CPTAC4_Gastric_Cancer_JHU_Proteome.sample.txt
```

Without `--sample_txt`, the R script uses existing annotation in `outdir` (e.g. `annotation_filled.csv`) or writes a template and exits; with `--sample_txt` it builds or corrects annotation from the CPTAC sample file.

### Option C: MSstatsTMT QC plots (for lab debugging)

The pipeline saves **MSstatsTMT’s built-in QC plot** (box plots of log intensities across channels and MS runs) to `results/{study_id}/plots/` by default. To disable: `--msstats_qc_plots FALSE`. To also save ProfilePlots for a few proteins: `--n_profile_proteins 5`.

To regenerate only MSstatsTMT’s native plots after the pipeline (no custom plots):

```bash
Rscript --no-init-file msstatsTMT_qc_plots.R --outdir results/PDC000120
Rscript --no-init-file msstatsTMT_qc_plots.R --outdir results/PDC000120 --n_profile_proteins 5
```

Outputs go to `results/{study_id}/plots/` (QCplot.pdf and optionally ProfilePlot per protein).

### Option D: CCLE (or other pre-built MSstats input)

For **CCLE** peptide-level TSV + sample table (Sheet2), use the converter in `data/ccle_peptide/` to produce `msstats_input.tsv` and `annotation_filled.csv`, then run the **same** R pipeline with `--msstats_input_dir`. CPTAC flow is unchanged (no `--msstats_input_dir` = normal PSM + annotation).

1. Export sample sheet: `python3 ccle_peptide/export_sample_sheet2_csv.py`
2. Convert: `python3 ccle_peptide/ccle_to_msstats_input.py --tsv ... --sample_csv ccle_peptide/sample_info_ccle.csv --outdir results/CCLE`
3. Run R: `Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R --msstats_input_dir results/CCLE --outdir results/CCLE`

See **data/ccle_peptide/README_CCLE.md** for details.

---

## R pipeline: main options

| Argument | Description |
|----------|-------------|
| `--psm_dir` | Directory of `.psm` files (e.g. `pdc_psm/PDC000120`). Required unless `--msstats_input_dir` is set. |
| `--msstats_input_dir` | Use pre-built `msstats_input.tsv` (e.g. from CCLE converter). Skips PSM parse and annotation. |
| `--outdir` | Output directory (e.g. `results/PDC000120`). All outputs and **per-study annotation** go here. |
| `--sample_txt` | Path to CPTAC `*.sample.txt`. If set, annotation is built/audited/rebuilt from it. |
| `--annotation` | Path to existing annotation CSV (Run, Channel, Condition, BioReplicate, Mixture, Fraction, TechRepMixture). |
| `--reference_channel` | Bridge channel (e.g. `131`, `126C`). Used when no annotation/sample.txt to auto-fill Norm. |
| `--species` | `Hs` or `Mm` for gene mapping. Default: Hs. |
| `--force_parse` | Re-parse PSM even if `parsed_psm_long.tsv` exists. |
| `--force_summarize` | Re-run protein summarization even if `protein_summary.tsv` exists. |

---

## Outputs per study (`results/{study_id}/`)

| File | Description |
|------|-------------|
| **annotation_filled.csv** | Filled annotation used for this run (or from `--annotation`). |
| **annotation_filled_corrected.csv** | Annotation built/corrected from `*.sample.txt` when using `--sample_txt`. |
| **annotation_template.csv** | Written only if no annotation and no `--reference_channel` (then script exits). |
| **annotation_audit.txt** | Audit of annotation vs sample.txt (when `--sample_txt` used). |
| **normalization_audit.txt** | Norm channel / POOL validation (when `--sample_txt` used). |
| **parsed_psm_long.tsv** | Long-format PSM (Run × Channel × PSM). |
| **msstats_input.tsv** | Input to MSstatsTMT `proteinSummarization`. |
| **protein_summary.tsv** | Protein-level abundances from MSstatsTMT. |
| **gene_matrix.csv** | Sample × gene (gene symbols), Norm channel removed. |
| **qc_summary.txt** | Counts and intensity summary. |
| **plots/** | MSstatsTMT QCPlot (and optional ProfilePlots) for lab debugging. |

So **each study has its own unique annotation files** in its own `results/{study_id}/` folder.

### Expected output dimensions

Typical outputs per CPTAC study:

- **Samples:** ~80–150 (exact number depends on study design and TMT plex size).
- **Genes:** ~4000–6000 (depends on study coverage and missingness).

Exact values vary by study. The pipeline does not filter genes by missingness; downstream preprocessing may apply such filters.

---

## Disk space

Pipeline outputs can use **a lot of disk space** on a laptop. Approximate sizes per CPTAC study:

| What | Typical size per study |
|------|------------------------|
| **Raw PSM downloads** (`pdc_psm/{study_id}/`) | ~2–4 GB |
| **parsed_psm_long.tsv** | ~2–10 GB |
| **msstats_input.tsv** | ~2–3 GB |
| **protein_summary.tsv** | ~2–4 GB |
| **gene_matrix.csv** | ~50–150 MB |
| Annotation / QC files | under 50 MB |

For **12 studies**, keeping everything can mean **~50–100+ GB** total. To save space:

- **Minimal retention:** After the pipeline finishes for a study, you can keep only what you need for downstream ML and delete the rest. The only file strictly required for modeling is **`gene_matrix.csv`**. You may also want to keep **annotation_filled.csv** (or **annotation_filled_corrected.csv**) and optionally **protein_summary.tsv** for QC.
- **Delete heavy intermediates:** Removing `parsed_psm_long.tsv` and `msstats_input.tsv` from `results/{study_id}/` frees several GB per study. You cannot re-run summarization without re-running the full R pipeline (and re-parsing PSM).
- **Delete PSM after success:** Once `gene_matrix.csv` exists for a study, you can delete `pdc_psm/{study_id}/` to free ~2–4 GB per study. Re-running that study would require re-downloading from PDC.

Use the optional script **`cleanup_study_disk.sh`** to remove heavy intermediates (and optionally PSM) for one or all studies after the pipeline has produced `gene_matrix.csv`. See the script header for usage.

---

## Downstream preprocessing (not part of this pipeline)

After generating gene matrices with this pipeline, **downstream analysis pipelines** may perform:

- Missing value imputation
- Feature filtering (e.g. by missingness or variance)
- Cross-study normalization
- Batch correction
- Dimensionality reduction
- Machine learning model training

These steps are **intentionally separated** from this pipeline. This pipeline only produces the sample × gene matrices; all ML preprocessing and modeling happen in later workflows.

---

## Helper scripts

### Check manifests (studies and duplicates)

```bash
python3 check_manifests.py
```

- Lists which study is in each manifest file.
- Lists all unique studies.
- Reports duplicate **File ID** / **File name** across manifests and which studies appear in more than one manifest file.

### Check studies vs sample file registry

```bash
python3 check_studies_sample_file.py --manifests manifests
# or after downloads, from pdc_psm:
python3 check_studies_sample_file.py --psm_dir pdc_psm
```

- Ensures every study (from manifests or `pdc_psm/`) has a row in `sample_files_msstats_tmt.csv`.
- Reports `path_exists=yes/no` for each study. Exit code 1 if any are missing.

### Single-study check

```bash
python3 check_studies_sample_file.py --study PDC000614
```

Exit 0 if that study has an entry in the CSV, else 1.

---

## Typical workflow

1. **Export manifests** from PDC (one manifest per study or one per download batch).
2. Put them in `manifests/` and **remove duplicates** (keep one manifest per study); use `check_manifests.py` to see studies and duplicates.
3. Ensure **every study** has a row in `sample_files_msstats_tmt.csv` with the correct path to its `.sample.txt`; use `check_studies_sample_file.py`.
4. Run **`./run_pipeline_per_manifest.sh`** to download and process each manifest.
5. QC plots: the pipeline writes MSstatsTMT’s built-in QC plot to `results/{study_id}/plots/`. Optionally run **`msstatsTMT_qc_plots.R`** to regenerate those plots only.

To use **sample.txt-driven annotation** for each study, extend the runner to read `sample_files_msstats_tmt.csv`, look up `path` for the current `study_id`, and pass `--sample_txt "$path"` to the R script. Until then, the R pipeline uses whatever annotation already exists in `results/{study_id}/` or a template/reference_channel fallback.

---

## References

- [MSstatsTMT](https://bioconductor.org/packages/MSstatsTMT/) – protein summarization and normalization for TMT.
- [PDC](https://pdc.cancer.gov/) – Proteomic Data Commons.
- CPTAC `*.sample.txt`: study design (FileNameRegEx, AnalyticalSample, TMT channel columns, POOL = bridge channel).
