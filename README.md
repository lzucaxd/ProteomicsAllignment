# ProteomicsAllignment

Pipeline to convert PDC (Proteomic Data Commons) PSM-level proteomics data into sample × gene matrices using **MSstatsTMT** (R) and Python for download and orchestration. Supports TMT6, TMT10, TMT11, TMT16/TMTpro, and TMT18; channels are detected from the data.

## Requirements

- **Python 3.8+** — downloads (PDC manifests, PSM files)
- **R 4.0+** — TMT normalization, protein summarization, gene mapping, and optional QC plots

## Setup (one-time)

Run from the **repo root**:

```bash
# 1. Python dependencies
pip install -r requirements.txt
# Or with a virtualenv:
#   python3 -m venv .venv && source .venv/bin/activate  # Windows: .venv\Scripts\activate
#   pip install -r requirements.txt

# 2. R dependencies (installs MSstatsTMT, org.Hs.eg.db, data.table, ggplot2, tidyr)
Rscript install_r_packages.R
```

R packages installed by the script:

| Source         | Packages |
|----------------|----------|
| **CRAN**       | `data.table`, `ggplot2`, `tidyr` |
| **Bioconductor** | `MSstatsTMT`, `org.Hs.eg.db` |

For **mouse** gene mapping, uncomment the `org.Mm.eg.db` line in `install_r_packages.R` and run it again. Package list (for reference): [requirements-R.txt](requirements-R.txt).

## Run the pipeline

All pipeline steps are run from the **`data/`** directory:

```bash
cd data
./run_pipeline_per_manifest.sh
```

This downloads PSM files from each manifest in `manifests/`, then runs the R pipeline (MSstatsTMT → gene matrix) per study. Outputs go to `data/results/{study_id}/`.

**First time?** Put your PDC manifest CSVs in `data/manifests/` and ensure each study has a row in `data/sample_files_msstats_tmt.csv` (see [data/PIPELINE_README.md](data/PIPELINE_README.md)).

## Documentation

Full pipeline documentation, options, and usage: **[data/PIPELINE_README.md](data/PIPELINE_README.md)**

## Layout

| Path | Purpose |
|------|---------|
| `requirements.txt` | Python deps (`requests`) |
| `requirements-R.txt` | List of R packages (reference) |
| `install_r_packages.R` | One-command R setup |
| `data/` | Pipeline scripts; **run from here** |
| `data/manifests/` | PDC file manifest CSVs |
| `data/pdc_psm/` | Downloaded PSM files (generated, git-ignored) |
| `data/results/` | Gene matrices and outputs (generated, git-ignored) |
