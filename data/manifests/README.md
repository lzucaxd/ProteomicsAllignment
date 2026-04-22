# PDC file manifests (local only)

**Committed manifests are not stored in git.** PDC “File Download Link” URLs **expire** (often within days), so exports are **not reproducible** as long-lived repo artifacts.

## What to put here

1. In [NCI Proteomic Data Commons](https://pdc.cancer.gov/) (or [Data Commons](https://proteomic.datacommons.cancer.gov/)), open your study (e.g. **PDC000120**, **PDC000153**).
2. Go to **Files**.
3. Filter **Data Category** = **Peptide Spectral Matches** and **File Type** = **Text** (files ending in **`.psm`**). Do **not** export a manifest that lists only **`.mzid.gz`** — the downloader will match zero PSM tables.
4. **Export** the file manifest as **CSV** and save it under this folder.

Naming (recommended so per-study scripts find your file without environment variables):

| Study   | Suggested filename                 |
|---------|------------------------------------|
| PDC000127 | `PDC000127_pdc_file_manifest.csv` |
| PDC000153 | `PDC000153_pdc_file_manifest.csv` |

You can also keep PDC’s default export name (`PDC_file_manifest_YYYYMMDD_*.csv`); **`run_pipeline_per_manifest.sh`** picks up any `PDC_file_manifest_*.csv` or `PDC*_pdc_file_manifest.csv` in this directory.

## Example layout (committed)

- **`example_pdc_file_manifest.csv`** — same **columns** as a real PDC export, with **placeholder** download URLs. Use it only to see the CSV shape; **replace** with your own export before downloading.

## Checks

```bash
cd data
python3 check_manifests.py
python3 check_studies_sample_file.py
```

Full pipeline: **[`../PIPELINE_README.md`](../PIPELINE_README.md)**.
