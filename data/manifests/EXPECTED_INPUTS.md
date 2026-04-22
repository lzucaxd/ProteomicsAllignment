# Expected inputs (manifest → PSM → gene matrix)

Large files are **not** committed. After a fresh clone, obtain the following **locally** (or on shared lab storage).

## 1. PDC file manifest CSVs

- **What:** CSV export from the PDC portal listing **Peptide Spectral Matches** as **Text** / **`.psm`** rows (not mzIdentML-only exports).
- **Where to save:** `data/manifests/` — see **`README.md`** in this folder for naming (`PDC*_pdc_file_manifest.csv` or portal default `PDC_file_manifest_*.csv`).
- **Produced by:** you, from [PDC](https://pdc.cancer.gov/) or [Data Commons](https://proteomic.datacommons.cancer.gov/).
- **Consumed by:** `data/pdc_manifest_downloader.py`, `data/run_pipeline_per_manifest.sh`, per-study `run_study_*.sh`.

## 2. CPTAC MSstatsTMT design (`*.sample.txt`)

- **What:** Per-study TMT channel / mixture / bridge annotation.
- **Registry:** `data/sample_files_msstats_tmt.csv` (paths often mirror-relative).
- **Consumed by:** `data/pdc_psm_to_msstatsTMT_protein_matrix.R`.

## 3. CCLE inputs (if running CCLE arm)

- **What:** Peptide-level tables and processed gene matrix per your lab pipeline.
- **Benchmark default path:** `data/results/CCLE_corrected/gene_matrix.csv` (override in `configs/preprocessing/default.yaml`).
- **Docs:** `data/ccle_peptide/` READMEs if present.

## 4. Outputs this stage produces (gitignored by default)

| Path | Produced by |
|------|-------------|
| `data/pdc_psm/{study}/` | `pdc_manifest_downloader.py` |
| `data/results/{study}/gene_matrix.csv` | `pdc_psm_to_msstatsTMT_protein_matrix.R` |

Full pipeline narrative: **`../PIPELINE_README.md`**.

## 5. Downstream (benchmark) inputs

After `gene_matrix.csv` exist for CPTAC + CCLE, the **overnight benchmark** reads paths from **`configs/preprocessing/default.yaml`** and writes union matrices under `data/processed/union/` (see `scripts/benchmark/run_overnight_v2.sh`).
