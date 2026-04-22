# PDC000127 — CPTAC3 ccRCC (Clear Cell Renal Cell Carcinoma) discovery proteome

MSstatsTMT/TMT10 pipeline expects this file **in this folder**:

`CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.sample.txt`

## How to obtain it

1. **From your existing CPTAC mirror** (if you have one), copy or symlink the study’s `.sample.txt` into this directory under the exact name above.
2. **From NCI PDC:** open study **PDC000127** (CPTAC CCRCC Discovery — Proteome), download the study **sample annotation** / design file that lists TMT plexes, channels, and `POOL`/bridge rows (same style as other CPTAC3 `.sample.txt` bundles).

`data/sample_files_msstats_tmt.csv` points at `cptac_samples/PDC000127/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.sample.txt` (path is **relative to `data/`**).

## Run the pipeline

From **`data/`**:

```bash
./run_study_PDC000127.sh
```

Optional:

- `PDC000127_SAMPLE_TXT=/path/to/…sample.txt` — override sample file location.
- `PDC000127_MANIFEST=/path/to/your_export.csv` — manifest containing **only** PDC000127 PSM rows (default: `manifests/PDC000127_pdc_file_manifest.csv` after you save a PDC export there; see `data/manifests/README.md`).
- `MAX_PSM_FILES=20` — download only the first *N* `.psm` files (quick sanity check; not for production).

Outputs: `results/PDC000127/` — `gene_matrix.csv`, `protein_summary.tsv`, `msstats_input.tsv`, etc.  
PSM downloads: `pdc_psm/PDC000127/Peptide_Spectral_Matches/*.psm` (~575 `.psm` rows for the full manifest filter).

**Requirements:** Python `requests`, R + **MSstatsTMT**, network for PDC downloads.

## Manifest must list `.psm` files (not only `.mzid.gz`)

PDC often offers **Peptide Spectral Matches** as:

- **mzIdentML** (`.mzid.gz`) — *not* used by `pdc_psm_to_msstatsTMT_protein_matrix.R` in this repo.
- **Text / TSV** (`.psm`) — **required** for the parser.

If your exported manifest has **575 rows all ending in `.mzid.gz`**, the downloader will match **0** `.psm` files and `run_study_PDC000127.sh` will **exit with an error** after preflight.

**Fix:** In [PDC](https://proteomic.datacommons.cancer.gov/), open **PDC000127** → **Files**, filter **Data Category** = Peptide Spectral Matches and **File Type** = **Text** (or file name `*.psm`), then **export manifest** again. You should see **~575** rows with names like `15CPTAC_CCRCC_..._LUMOS_f01.psm`.

## HTTP 403 on download (expired manifest URLs)

PDC file manifests embed **time-limited** CloudFront download links. If `pdc_manifest_downloader.py` prints `HTTP 403` for every file, the manifest CSV in `data/manifests/` is **stale**.

**Fix:** In the [NCI Proteomic Data Commons](https://proteomic.datacommons.cancer.gov/) portal, open study **PDC000127**, go to **Files**, filter to **Peptide Spectral Matches** / `.psm` as needed, and **export a new file manifest** (CSV). Save it under `data/manifests/` (e.g. `PDC_file_manifest_PDC000127_fresh.csv`) and run:

```bash
PDC000127_MANIFEST="$PWD/manifests/PDC_file_manifest_PDC000127_fresh.csv" ./run_study_PDC000127.sh
```

Alternatively, place already-downloaded `.psm` files under `data/pdc_psm/PDC000127/` (recursive subfolders are fine) and run only the R step:

```bash
Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R \
  --psm_dir pdc_psm/PDC000127 \
  --outdir results/PDC000127 \
  --sample_txt cptac_samples/PDC000127/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.sample.txt \
  --replace_annotation
```
