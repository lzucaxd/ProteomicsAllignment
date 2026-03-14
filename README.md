# ProteomicsAllignment

Pipeline to convert PDC (Proteomic Data Commons) PSM-level proteomics data into sample × gene matrices using MSstatsTMT.

## Quick start

- **Python:** `pip install -r requirements.txt` (or use a venv: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`)
- **R:** MSstatsTMT, org.Hs.eg.db (or org.Mm.eg.db), data.table, ggplot2, tidyr (see [data/PIPELINE_README.md](data/PIPELINE_README.md))

## Documentation

Full pipeline documentation, usage, and options: **[data/PIPELINE_README.md](data/PIPELINE_README.md)**

## Layout

- `data/` — Pipeline scripts, manifests, and config; run from here.
  - Download PSM files from PDC manifests → MSstatsTMT → gene matrices.
  - `pdc_psm/` and `results/` are generated (ignored by git; see `.gitignore`).
