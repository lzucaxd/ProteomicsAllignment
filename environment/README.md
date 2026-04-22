# Environment setup

This project uses **Python ≥ 3.10** and **R ≥ 4.1** (Bioconductor + CRAN).

## Python

From the **repository root**:

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -U pip
pip install -r requirements.txt
pip install -e .
# Optional Celligner extras:
# pip install -e ".[celligner]"
```

**Declared dependencies:** `pyproject.toml` (package `harmonize`). **`requirements.txt`** is kept aligned for `pip install -r` workflows.

**Verify:**

```bash
python3 scripts/verify_repro_setup.py
python3 scripts/verify_repro_setup.py --require-data   # after gene_matrix.csv exist
```

## R

From the **repository root**:

```bash
Rscript install_r_packages.R
```

Installs CRAN packages (`data.table`, `ggplot2`, `limma`, …) and Bioconductor **`MSstatsTMT`**, **`org.Hs.eg.db`**, etc. See the script header for details.

## External requirements

- **Disk:** CPTAC PSM downloads and `gene_matrix.csv` outputs are large (tens–hundreds of GB possible).
- **Network:** PDC downloads; manifest URLs **expire** — export fresh CSVs (`data/manifests/README.md`).
- **`CPTAC_LOCAL_MIRROR`:** optional but typical for `data/run_batch_studies.sh` (see `docs/LAB_ONBOARDING.md`).

## Optional / heavy

- **Celligner:** `pip install -e ".[celligner]"` plus vendored tree under `models/` per `docs/METHODS.md` and `scripts/benchmark/README.md`.
- **Conda:** not required; a venv is enough. If your lab standardizes on conda, mirror the same package list from `pyproject.toml`.
