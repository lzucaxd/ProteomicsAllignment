# How to run everything

**This is the single entry-point guide** for going from a **fresh clone** to **benchmark outputs**, including **where to point the pipeline at your files**. Deep design and methods live elsewhere (see [Documentation map](#documentation-map) at the end).

| If you need… | Jump to |
|--------------|---------|
| Clone, install, verify, what to commit for papers | [Clean clone and checks](#1-clone-install-and-verify) · [CLEAN_CLONE_REPRODUCIBILITY.md](CLEAN_CLONE_REPRODUCIBILITY.md) |
| PDC manifests → PSM → `gene_matrix.csv` (front-door table) | [`../pipeline/psm_to_gene_matrix/README.md`](../pipeline/psm_to_gene_matrix/README.md) · [PSM / MSstatsTMT preprocessing](#3-preprocessing--psm--gene_matrixcsv) |
| MSstatsTMT vs limma (what runs where) | [`INFERENCE_BASELINES.md`](INFERENCE_BASELINES.md) |
| Full benchmark after matrices exist | [Overnight benchmark](#4-harmonization-benchmark) |
| Your paths / other studies (same matrix shape) | [Running on your own data](#6-running-on-your-own-data-paths) |
| New task or new harmonization method | [HANDOFF_CHECKLIST.md](HANDOFF_CHECKLIST.md) · [config_system_overview.md](config_system_overview.md) |

---

## 0) What you are running (two layers)

1. **Heavy CPTAC/CCLE acquisition and PSM → protein/gene matrices** — mostly under **`data/`** (see [`data/PIPELINE_README.md`](../data/PIPELINE_README.md)). **Front-door file map:** [`../pipeline/psm_to_gene_matrix/README.md`](../pipeline/psm_to_gene_matrix/README.md) — this is **not** the exploratory code under **`data/scripts/`** (see [`../data/scripts/README.md`](../data/scripts/README.md)).
2. **Harmonization benchmark** (shared gene space, methods, limma, nulls, ceilings, tables) — orchestrated by **`scripts/benchmark/run_overnight_v2.sh`** from the **repo root**.

**Inference map:** MSstatsTMT builds CPTAC matrices; limma evaluates harmonized representations — **[`INFERENCE_BASELINES.md`](INFERENCE_BASELINES.md)**.

Large files are **not** in git; you supply data and manifests as below.

---

## 1) Clone, install, and verify

```bash
git clone <repo-url> ProteomicsAllignment
cd ProteomicsAllignment

python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -U pip
pip install -r requirements.txt
pip install -e .

# Optional: Celligner extras
# pip install -e ".[celligner]"

Rscript install_r_packages.R
```

**Sanity check (no large data):**

```bash
python3 scripts/verify_repro_setup.py
```

After `gene_matrix.csv` files exist at the paths in `configs/preprocessing/default.yaml`:

```bash
python3 scripts/verify_repro_setup.py --require-data
```

**Environment variables** (typical; see [`LAB_ONBOARDING.md`](LAB_ONBOARDING.md) for detail):

| Variable | Purpose |
|----------|---------|
| `PYTHONPATH` | Include `src/` for ad-hoc Python. **`run_overnight_v2.sh` sets this.** |
| `CPTAC_LOCAL_MIRROR` | Parent of your local `PDC000120/`, … tree for `.sample.txt` paths in `data/sample_files_msstats_tmt.csv`. |
| `PROTEOMICS_ALIGNMENT_ROOT` | Rare: absolute repo root for R if auto-detection fails. |
| `R_PROFILE_USER=/dev/null` | If a broken parent **renv** breaks `Rscript`. |

Reproducibility / commit policy for paper runs: **[CLEAN_CLONE_REPRODUCIBILITY.md](CLEAN_CLONE_REPRODUCIBILITY.md)**.

---

## 2) Data you must supply (nothing auto-downloads on clone)

| What | Where it comes from | Where it goes in this repo |
|------|---------------------|----------------------------|
| **PDC file manifest CSV** per study | [NCI PDC](https://pdc.cancer.gov/) → study → **Files** → filter to **PSM** / **`.psm`** → **Export** manifest | `data/manifests/` — **how to name and refresh:** [`data/manifests/README.md`](../data/manifests/README.md) |
| **CPTAC MSstatsTMT design** (`.sample.txt`) | Study bundle or lab mirror | Paths in **`data/sample_files_msstats_tmt.csv`** (column `path`); resolve with **`CPTAC_LOCAL_MIRROR`** if relative |
| **CCLE** inputs + processed matrix | Your lab pipeline | Benchmark expects **`data/results/CCLE_corrected/gene_matrix.csv`** (and inputs under `data/ccle_peptide/` as documented) |

**Manifest URLs expire** (~often a week); refresh the manifest from PDC if downloads return **HTTP 403**. The manifest must list **`.psm`** files, not only `.mzid.gz`. See root **`README.md`** (repo home) for the same table with footnotes, and **`data/cptac_samples/`** READMEs for study-specific notes.

**Checks before long CPTAC runs:**

```bash
cd data
python3 check_manifests.py
python3 check_studies_sample_file.py
```

---

## 3) Preprocessing — PSM → `gene_matrix.csv`

From **`data/`** after manifests + `sample_files_msstats_tmt.csv` + mirror are correct:

```bash
cd data

# All manifests in data/manifests/
./run_pipeline_per_manifest.sh
# Optional space-saving: ./run_pipeline_per_manifest.sh --cleanup-after

# Or batch by study with mirror set:
# export CPTAC_LOCAL_MIRROR=/path/to/parent/of/PDC000120
./run_batch_studies.sh
```

Authoritative detail: **[`data/PIPELINE_README.md`](../data/PIPELINE_README.md)**. Script map: **[`scripts/preprocessing/README.md`](../scripts/preprocessing/README.md)**.

**Minimum inputs for the overnight benchmark** (default layout):

- `data/results/PDC000120/gene_matrix.csv`
- `data/results/PDC000153/gene_matrix.csv`
- `data/results/CCLE_corrected/gene_matrix.csv`

You may **change those paths** in **`configs/preprocessing/default.yaml`** (see [§6](#6-running-on-your-own-data-paths)).

---

## 4) Harmonization benchmark

From **repository root** (after §3 matrices exist):

```bash
export PYTHONPATH="${PWD}/src${PYTHONPATH:+:$PYTHONPATH}"   # optional if using the shell script alone
bash scripts/run_benchmark.sh
# equivalent: bash scripts/benchmark/run_overnight_v2.sh
```

This runs annotation prep, union matrices, methods, limma DA, cross-domain metrics, permutations, ceilings, calibration, and assembly. **Steps involving permutations/ceilings can take hours.** Step-by-step reference: **[`scripts/benchmark/README.md`](../scripts/benchmark/README.md)**.

**Primary outputs:**

- `reports/benchmark_master/benchmark_results/comparison_summary.csv`
- `reports/benchmark_master/benchmark_results/disconnect_scores.csv`
- `reports/benchmark_master/logs/overnight_v2_*.log`

```bash
test -f reports/benchmark_master/benchmark_results/comparison_summary.csv && \
  head -5 reports/benchmark_master/benchmark_results/comparison_summary.csv
```

---

## 5) Optional — tables, LaTeX, figures, slides

| Goal | Command (repo root) |
|------|---------------------|
| Markdown bundle for papers / LLMs | `Rscript --vanilla scripts/benchmark/build_final_benchmark_tables.R` → `reports/benchmark_master/final_tables/BENCHMARK_NUMBERS_MASTER_FOR_LLM.md` |
| Report-style figures | `Rscript --vanilla scripts/presentation/generate_report_figures.R` |
| Slide pack | **`scripts/presentation/prepare_all.sh`** and **[`BENCHMARK_V2_AND_PRESENTATION.md`](BENCHMARK_V2_AND_PRESENTATION.md)** |

**Legacy Python chain** (debug only): `scripts/run_all.py` — prefer **`run_overnight_v2.sh`** for full runs.

---

## 6) Running on your own data (paths and limits)

### Same tasks, your file locations

1. Edit **`configs/preprocessing/default.yaml`**: `data_sources.cptac.*.gene_matrix`, `data_sources.ccle.gene_matrix`, and `sample_info` as needed.
2. For **bridge** harmonization, align **`configs/methods/bridge_aware.yaml`** `bridge_extraction` paths to your **`msstats_input.tsv`** files if they differ from defaults.
3. Re-run **`run_overnight_v2.sh`** (or **`python scripts/run_preprocessing.py --config … --output-dir …`** for matrix-only experiments).

### Custom studies or new biological contrasts

- **New CPTAC study IDs in YAML:** add blocks under `data_sources.cptac`; wire studies into **`configs/tasks/breast_subtype.yaml`** or **`breast_vs_lung.yaml`** as appropriate.
- **New benchmark task** (new contrast name / new metadata rules): YAML is not enough — extend **`scripts/run_preprocessing.py`**, **`src/harmonize/preprocessing/metadata.py`** (`build_sample_meta`), and **`scripts/benchmark/regenerate_methods_union.py`**. Follow **[HANDOFF_CHECKLIST.md](HANDOFF_CHECKLIST.md)**.

### Config reference

- **[config_system_overview.md](config_system_overview.md)** — YAML layout and which runners read which files.

---

## Documentation map

| Doc | Role |
|-----|------|
| **This file** | **Single run guide** (you are here) |
| [`README.md`](../README.md) | Repo home: question, findings, layout, full doc table |
| [`CLEAN_CLONE_REPRODUCIBILITY.md`](CLEAN_CLONE_REPRODUCIBILITY.md) | Verify script, commit/Zenodo policy |
| [`HANDOFF_CHECKLIST.md`](HANDOFF_CHECKLIST.md) | Lab checklist, extend methods/tasks |
| [`LAB_ONBOARDING.md`](LAB_ONBOARDING.md) | Mirror paths, env vars |
| [`data/PIPELINE_README.md`](../data/PIPELINE_README.md) | PDC → matrix pipeline (long form) |
| [`../pipeline/psm_to_gene_matrix/README.md`](../pipeline/psm_to_gene_matrix/README.md) | Same pipeline — **entry-point table** (shell + R) |
| [`INFERENCE_BASELINES.md`](INFERENCE_BASELINES.md) | MSstatsTMT vs limma — where each runs |
| [`NAMING_AND_PATHS.md`](NAMING_AND_PATHS.md) | Why paths are scattered; naming rules for new files |
| [`scripts/benchmark/README.md`](../scripts/benchmark/README.md) | Overnight steps 0–n |
| [`METHODS.md`](METHODS.md) | Raw, bridge, Celligner definitions |
| [`END_TO_END_TECHNICAL_REPORT.md`](END_TO_END_TECHNICAL_REPORT.md) | Paper-length narrative |

---

## See also

- **[`docs/README.md`](README.md)** — full index of topic notes and older reports.
