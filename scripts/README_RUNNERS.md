# Runnable entry points (quick reference)

| Script | What it runs |
|--------|----------------|
| **`run_benchmark.sh`** | Full overnight benchmark (`run_overnight_v2.sh`): preprocessing union, methods, limma, nulls, ceilings, tables, meeting figures. |
| **`run_diagnostics.sh`** | Subset: `preflight` (R), `structure` (Python batch), or `all`. Env: `PROCESSED_UNION_DIR`, `META_DIR`, `INTERSECTION_DIR` override defaults. |
| **`run_methods.py`** / **`run_methods.sh`** | Python harmonize method registry smoke test (not the overnight R union matrix rebuild). For union matrices: `python scripts/benchmark/regenerate_methods_union.py`. |
| **`run_benchmark.py`** | Python benchmark driver (legacy / debugging); prefer **`run_benchmark.sh`** for parity with lab paper runs. |
| **`run_preprocessing.py`** | Union / intersection shared matrices from existing `gene_matrix.csv` files. |

Orchestration detail: **`benchmark/README.md`**.
