# Clean clone and reproducibility

**Primary run guide (all steps in one place):** **[HOW_TO_RUN_EVERYTHING.md](HOW_TO_RUN_EVERYTHING.md)**.

This page focuses on **clone → install → verify → commit policy** after `git clone`. Large data are **not** in the repository; reproducibility means **same code + same inputs → same outputs** once those inputs are placed as documented.

## 1. Clone and enter the repo

```bash
git clone <repo-url> ProteomicsAllignment
cd ProteomicsAllignment
```

The directory name `ProteomicsAllignment` is historical; keep it or adjust paths in your own notes.

## 2. Python environment (required)

**Supported:** Python **≥ 3.10** (see `pyproject.toml` → `harmonize`).

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -U pip
pip install -r requirements.txt
```

For an editable install of the `harmonize` package (same dependency set as declared in `pyproject.toml`):

```bash
pip install -e .
```

Optional Celligner extras:

```bash
pip install -e ".[celligner]"
```

## 3. R packages (required for overnight benchmark)

```bash
Rscript install_r_packages.R
```

Uses **R ≥ 4.1** in practice (script header mentions R ≥ 4.0). CRAN + Bioconductor: `limma`, `MSstatsTMT`, `data.table`, `ggplot2`, etc.

## 4. Verify the machine (no data required)

From **repository root**:

```bash
python3 scripts/verify_repro_setup.py
```

With CPTAC/CCLE `gene_matrix.csv` files in place (paths from `configs/preprocessing/default.yaml`):

```bash
python3 scripts/verify_repro_setup.py --require-data
```

**Continuous integration:** `.github/workflows/repro_check.yml` runs `verify_repro_setup.py --skip-r` on each push/PR (Python stack only; install R locally for the full pipeline).

## 5. Data you must supply (not in git)

See root **`README.md`** §0 and **`data/README.md`**. Minimum for the **overnight benchmark** after matrices exist:

| Path | Role |
|------|------|
| `data/results/PDC000120/gene_matrix.csv` | CPTAC breast |
| `data/results/PDC000153/gene_matrix.csv` | CPTAC lung |
| `data/results/CCLE_corrected/gene_matrix.csv` | CCLE |

Paths are configurable in **`configs/preprocessing/default.yaml`**. PDC manifests expire; refresh them as documented in the README.

## 6. Run the full benchmark pipeline

From repo root (overnight script sets `PYTHONPATH` and prefers `.venv/bin/python3` if present):

```bash
bash scripts/benchmark/run_overnight_v2.sh
```

Primary outputs:

- `reports/benchmark_master/benchmark_results/comparison_summary.csv`
- `reports/benchmark_master/benchmark_results/disconnect_scores.csv`
- `reports/benchmark_master/logs/overnight_v2_*.log`

## 7. What to commit for “reproducible paper numbers”

- **Always commit:** code, `configs/`, `install_r_packages.R`, `requirements.txt`, `pyproject.toml`, scripts, and documentation.
- **Policy choice:** whether to commit frozen **`reports/benchmark_master/benchmark_results/*.csv`** for a given publication tag. Large generated trees are often gitignored; if so, archive those CSVs (e.g. Zenodo) next to a **git tag** and record the tag in the supplement.

## 8. Extending tasks or methods

YAML alone is not always sufficient. See **`docs/HANDOFF_CHECKLIST.md`** and **`docs/config_system_overview.md`**. New benchmark tasks require Python changes in **`scripts/run_preprocessing.py`** and **`src/harmonize/preprocessing/metadata.py`** (`build_sample_meta`), plus updates to **`scripts/benchmark/regenerate_methods_union.py`** for the overnight path.

## See also

| Doc | Purpose |
|-----|---------|
| [`HANDOFF_CHECKLIST.md`](HANDOFF_CHECKLIST.md) | Checklists: reproduce, extend |
| [`LAB_ONBOARDING.md`](LAB_ONBOARDING.md) | `CPTAC_LOCAL_MIRROR`, layout |
| [`config_system_overview.md`](config_system_overview.md) | YAML layout |
| [`HOW_TO_RUN_EVERYTHING.md`](HOW_TO_RUN_EVERYTHING.md) | Central end-to-end run guide |
| [`../README.md`](../README.md) | Repo home, findings, doc map |
