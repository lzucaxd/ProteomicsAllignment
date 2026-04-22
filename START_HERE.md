# Start here (after `git clone`)

This file is the **default entry** for anyone who landed in the repo from **git** or GitHub and wants to run something without reading the whole doc tree.

---

## 1) Data pipeline — how to run (PSM → gene matrix)

This stage is **CPTAC / CCLE acquisition + MSstatsTMT → `gene_matrix.csv`**. It is **not** under `data/scripts/` (those are exploratory helpers).

| Step | What to open |
|------|----------------|
| **Map of every shell / R driver** | **[`pipeline/psm_to_gene_matrix/README.md`](pipeline/psm_to_gene_matrix/README.md)** ← read this first |
| **Full narrative (channels, QC, data flow)** | [`data/PIPELINE_README.md`](data/PIPELINE_README.md) |
| **What files you must place locally** | [`data/manifests/README.md`](data/manifests/README.md), [`data/manifests/EXPECTED_INPUTS.md`](data/manifests/EXPECTED_INPUTS.md) |
| **Curated subtype / biospecimen tables (in git)** | [`data/annotations/README.md`](data/annotations/README.md) |

**Run** (typical; working directory **`data/`**):

```bash
cd data
./run_pipeline_per_manifest.sh
# batch / mirror layout:
# export CPTAC_LOCAL_MIRROR=/path/to/parent/of/PDC000120
# ./run_batch_studies.sh
```

Mirror and env vars: [`docs/LAB_ONBOARDING.md`](docs/LAB_ONBOARDING.md).

---

## 2) Harmonization benchmark (after matrices exist)

From **repository root**:

```bash
bash scripts/run_benchmark.sh
```

Details: [`scripts/benchmark/README.md`](scripts/benchmark/README.md).

---

## 3) Single runbook (clone → install → data → benchmark)

**[`docs/HOW_TO_RUN_EVERYTHING.md`](docs/HOW_TO_RUN_EVERYTHING.md)** — one file, end to end.

---

## 4) Install once

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -U pip && pip install -r requirements.txt && pip install -e .
Rscript install_r_packages.R
python3 scripts/verify_repro_setup.py
```

More: [`environment/README.md`](environment/README.md).

---

## 5) “Too many docs?” — smallest reading order

**[`HANDOFF.md`](HANDOFF.md#if-this-feels-like-too-many-docs-at-once)** — what to read first vs later.

**MSstatsTMT (matrix build) vs limma (benchmark):** [`docs/INFERENCE_BASELINES.md`](docs/INFERENCE_BASELINES.md).

**Why folders look scattered:** [`docs/NAMING_AND_PATHS.md`](docs/NAMING_AND_PATHS.md).

---

## 6) Lab context and outputs index

- Project + caveats: [`HANDOFF.md`](HANDOFF.md)
- Tables / figures index: [`reports/final_report/README.md`](reports/final_report/README.md)
