# Expected inputs (manifest → PSM → gene matrix)

Large files are **not** committed. After a fresh clone, obtain the following **locally** (or on shared lab storage).

## 1. PDC file manifest CSVs

- **What:** CSV export from the PDC portal listing **Peptide Spectral Matches** as **Text** / **`.psm`** rows (not mzIdentML-only exports).
- **Where to save:** `data/manifests/` — see **`README.md`** in this folder for naming (`PDC*_pdc_file_manifest.csv` or portal default `PDC_file_manifest_*.csv`).
- **Produced by:** you, from [PDC](https://pdc.cancer.gov/) or [Data Commons](https://proteomic.datacommons.cancer.gov/).
- **Consumed by:** `data/pdc_manifest_downloader.py`, `data/run_pipeline_per_manifest.sh`, per-study `run_study_*.sh`.

## 2. CPTAC MSstatsTMT design (`*.sample.txt`) and per-run annotation

- **What:** Per-study TMT **channel / mixture / bridge** design (CPTAC `*.sample.txt` from the study bundle or PDC).
- **Registry (required):** **`data/sample_files_msstats_tmt.csv`** — for each `study_id`, the **`path`** column must point to that study’s `.sample.txt` (usually **mirror-relative**, e.g. `PDC000120/CPTAC2_Breast_....sample.txt`). Optional columns include **`annotation_path`** (template or starter annotation) and **`reference_channel`**.
- **On disk:** the `.sample.txt` files themselves must exist where `path` resolves (lab mirror + **`CPTAC_LOCAL_MIRROR`** when paths are not absolute). Some studies use paths under **`data/cptac_samples/<PDC_ID>/`**; see the `notes` column and per-study READMEs there.
- **Outputs (local, under `data/results/{study_id}/`):** the R pipeline writes or updates **MSstatsTMT annotation CSVs** such as `annotation_filled_corrected.csv` when driven by `--sample_txt`. Those are **run artifacts**, not the small benchmark tables under **`data/annotations/`** (subtype mapping for tasks).
- **Consumed by:** `data/pdc_psm_to_msstatsTMT_protein_matrix.R` (via runners that pass `--sample_txt`).

## 3. CCLE inputs (if running the peptide → matrix arm)

- **Peptide-level path:** peptide TSV + sample sheet under **`data/ccle_peptide/`** → convert to `msstats_input.tsv` + annotation, then run the **same** R script with **`--msstats_input_dir`** (see **`../PIPELINE_README.md`** “Option D: CCLE”).
- **Pre-built gene matrix only:** if you already have **`gene_matrix.csv`**, the **benchmark** reads the path in **`configs/preprocessing/default.yaml`** (default example: `data/results/CCLE_corrected/gene_matrix.csv`). That skips PSM/peptide ingestion for CCLE but is still a required **file** for the harmonization stage.
- **Benchmark subtype labels (small, in git):** Luminal/Basal line lists live under **`data/ccle/`** and `data/processed/` processed copies — used when building **task** metadata, not as a substitute for MSstatsTMT **`.sample.txt`** on the CPTAC side.
- **Docs:** `data/ccle_peptide/` READMEs if present; **`data/ccle/README.md`** for the curated subtype CSV.

## 4. Outputs this stage produces (gitignored by default)

| Path | Produced by |
|------|-------------|
| `data/pdc_psm/{study}/` | `pdc_manifest_downloader.py` |
| `data/results/{study}/gene_matrix.csv` | `pdc_psm_to_msstatsTMT_protein_matrix.R` |

Full pipeline narrative: **`../PIPELINE_README.md`**.

## 5. Downstream (benchmark) inputs

After `gene_matrix.csv` exist for CPTAC + CCLE, the **overnight benchmark** reads paths from **`configs/preprocessing/default.yaml`** and writes union matrices under `data/processed/union/` (see `scripts/benchmark/run_overnight_v2.sh`).
