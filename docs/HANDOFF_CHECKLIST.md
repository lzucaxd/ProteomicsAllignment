# Handoff checklist (Vitek Lab)

Use this when **cloning fresh**, **auditing reproducibility**, or **extending** the benchmark.

## Reproduce benchmark outputs (given gene matrices exist)

- [ ] Clone repository (`ProteomicsAllignment` directory name on disk).
- [ ] Python **≥ 3.10**: `pip install -r requirements.txt` (optional `.venv/`).
- [ ] R **≥ 4.1**: `Rscript install_r_packages.R`.
- [ ] Confirm **`data/results/PDC000120/gene_matrix.csv`**, **`PDC000153`**, **`CCLE_corrected`** exist (large; see [`data/README.md`](../data/README.md) if missing).
- [ ] `export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"`
- [ ] Run `bash scripts/benchmark/run_overnight_v2.sh`
- [ ] Inspect `reports/benchmark_master/benchmark_results/comparison_summary.csv`
- [ ] Inspect `reports/benchmark_master/benchmark_results/disconnect_scores.csv`
- [ ] Check log: `ls -t reports/benchmark_master/logs/ | head`

## Reproduce preprocessing (PDC PSM → gene matrices)

- [ ] Read [`scripts/preprocessing/README.md`](../scripts/preprocessing/README.md) and [`data/PIPELINE_README.md`](../data/PIPELINE_README.md).
- [ ] Place PDC manifest CSVs under `data/manifests/` and configure `data/sample_files_msstats_tmt.csv`.
- [ ] From `data/`, run `./run_pipeline_per_manifest.sh` **or** `./run_batch_studies.sh` (requires `CPTAC_LOCAL_MIRROR` for the batch path).
- [ ] Confirm each study folder contains `gene_matrix.csv`.

## Add a new harmonization method

- [ ] Read **“Adding a new harmonization method”** in [`scripts/benchmark/README.md`](../scripts/benchmark/README.md).
- [ ] Emit `data/processed/methods/{name}/transformed_{task}.csv` compatible with union sample ordering.
- [ ] Wire method into `scripts/benchmark/regenerate_methods_union.py` (+ configs if needed).
- [ ] Re-run overnight from **Step 2** (or full script).
- [ ] Verify new rows in `comparison_summary.csv`.

## Add a new benchmark task

- [ ] Copy and edit `configs/tasks/*.yaml` for the new task.
- [ ] Extend `scripts/run_preprocessing.py` / preprocessing YAML to emit union outputs.
- [ ] Re-run overnight from **Step 1**.

## First files to open

| File | Why |
|------|-----|
| [`README.md`](../README.md) | Entry point |
| [`comparison_summary.csv`](../reports/benchmark_master/benchmark_results/comparison_summary.csv) | All numeric results |
| [`run_overnight_v2.sh`](../scripts/benchmark/run_overnight_v2.sh) | Ordered steps |
| [`extract_bridge_summaries.R`](../scripts/benchmark/extract_bridge_summaries.R) | Bridge channel recovery |

## Path / doc verification (optional)

```bash
cd /path/to/ProteomicsAllignment
grep -roh "scripts/[^ ]*\\.\\(R\\|py\\|sh\\)" docs/ scripts/*/README.md README.md 2>/dev/null | sort -u | while read -r f; do
  test -f "$f" && echo "OK: $f" || echo "MISSING: $f"
done
for f in configs/preprocessing/default.yaml configs/tasks/breast_subtype.yaml configs/tasks/breast_vs_lung.yaml; do
  test -f "$f" && echo "OK: $f" || echo "MISSING: $f"
done
for f in reports/benchmark_master/benchmark_results/comparison_summary.csv reports/benchmark_master/benchmark_results/disconnect_scores.csv; do
  test -f "$f" && echo "OK: $f" || echo "MISSING: $f"
done
```
