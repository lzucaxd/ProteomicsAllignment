# Cleanup log (handoff refactor)

Conservative changes only; anything uncertain stays in place or goes to **`archive/`** / **`notebooks/exploratory/`** with this log updated.

## 2026-04-21 — Manifest hygiene (earlier commit `ea69be5`)

| Action | Target | Reason |
|--------|--------|--------|
| **Removed from git** | `data/manifests/PDC_file_manifest_*.csv` (dated exports) | PDC download URLs **expire**; not reproducible in version control |
| **Added** | `data/manifests/README.md`, `example_pdc_file_manifest.csv` | Document how to export manifests; show CSV shape |
| **Gitignore** | `data/manifests/PDC_file_manifest_*.csv`, `PDC*_pdc_file_manifest.csv` | Keep local exports private to each clone |
| **Gitignore** | `data/processed_union/`, `data/MSI VS MSS/` | Optional local / ad-hoc paths |
| **Updated** | `run_pipeline_per_manifest.sh`, `run_batch_studies.sh`, `run_study_PDC000127.sh`, `run_study_PDC000153.sh`, check scripts | Conventional per-study manifest names + clearer errors |

## 2026-04-21 — Handoff documentation + runners (this pass)

| Action | Target | Reason |
|--------|--------|--------|
| **Added** | Root `REPO_AUDIT.md`, `REPO_LAYOUT_PLAN.md`, `HANDOFF.md`, `HANDOFF_SANITY_CHECK.md`, `CLEANUP_LOG.md` | Lab handoff deliverables requested |
| **Added** | `environment/README.md` | Single place for Python + R setup |
| **Added** | `data/manifests/EXPECTED_INPUTS.md` | Expected inputs without versioning large files |
| **Added** | `reports/final_report/README.md`, `FIGURE_MANIFEST.md` | Index to final tables/figures without moving binaries |
| **Added** | `scripts/run_benchmark.sh`, `scripts/run_diagnostics.sh`, `scripts/run_methods.sh` | Stable names wrapping existing entry points |
| **Added** | `scripts/README_RUNNERS.md` | Table of runner scripts |
| **Added** | `configs/benchmark/README.md` | Points to `default.yaml` as canonical benchmark config |
| **Added** | `notebooks/exploratory/README.md` | Policy for exploratory work |
| **Moved** | `data/gene_matrix_exploration.ipynb` → `notebooks/exploratory/` | Not referenced by pipeline; reduce `data/` clutter |
| **Updated** | Root `README.md`, `docs/README.md`, `docs/HANDOFF_CHECKLIST.md`, `archive/README.md`, `.gitignore` | Cross-links; ignore editor swap files |

## Not done (intentionally — would break paths or science)

- **No mass move** of `data/pdc_psm_to_msstatsTMT_protein_matrix.R` or other `data/` drivers into `src/` (downstream paths and muscle memory).
- **No rename** of `reports/benchmark_master/` to `results/` (widespread R `repo-root` assumptions).
- **No deletion** of untracked local `reports/` trees on disk (may be paper artifacts); document policy in `HANDOFF.md`.

## 2026-04-21 — Pipeline front door + inference map (doc-only)

| Action | Target | Reason |
|--------|--------|--------|
| **Added** | `pipeline/README.md`, `pipeline/psm_to_gene_matrix/README.md` | Single “manifest → matrix” map without moving `data/` drivers |
| **Added** | `data/scripts/README.md` | Clarify `data/scripts/` ≠ reproducible PDC→matrix path |
| **Added** | `scripts/exploratory/README.md` | Policy: new one-offs here, not under `data/scripts/` |
| **Added** | `docs/INFERENCE_BASELINES.md` | Prominent MSstatsTMT vs limma (matrix build vs benchmark) |
| **Updated** | Root `README.md`, `HANDOFF.md`, `docs/README.md`, `docs/HOW_TO_RUN_EVERYTHING.md`, `data/PIPELINE_README.md`, `data/README.md`, `scripts/README.md`, `scripts/preprocessing/README.md`, `REPO_AUDIT.md`, `reports/final_report/README.md`, `docs/DOC_LINK_CHECK.txt`, `scripts/run_all.py` | Cross-links |
| **Not done** | `git rm` of `data/scripts/` | Many `reports/*.md` and `PROJECT_REPORT.md` still reference paths; archive+grep would be a separate PR |

## 2026-04-21 — Naming / path map (doc-only)

| Action | Target | Reason |
|--------|--------|--------|
| **Added** | `docs/NAMING_AND_PATHS.md` | Single place: mental model, historical quirks, conventions for new files |
| **Updated** | `README.md` (layout tree + link), `CONTRIBUTING.md`, `HANDOFF.md`, `REPO_LAYOUT_PLAN.md`, `docs/README.md`, `docs/DOC_LINK_CHECK.txt` | Point contributors at the map |

## 2026-04-21 — Track annotation assets + CPTAC subtype mapping path

| Action | Target | Reason |
|--------|--------|--------|
| **Added** | `data/annotations/` + `data/annotations/cptac/PDC000120/gene_matrix_subtype_mapping.csv` | Version **subtype / column mapping** outside gitignored `data/results/` |
| **Updated** | `configs/preprocessing/*.yaml`, `configs/tasks/breast_subtype.yaml`, `src/harmonize/utils/paths.py`, R benchmark scripts, `build_PDC000120_subtype_mapping.py` | Point preprocessing at tracked mapping; script writes annotations + mirrors under `results/` when present |
| **Added / tracked** | `data/biospecimen/*`, `data/biospecimen/README.md`, `data/ccle/README.md`, `data/cptac_samples/` study README | Push biospecimen + CCLE label provenance the user asked to share |

## 2026-04-21 — START_HERE for git clones + README pipeline emphasis

| Action | Target | Reason |
|--------|--------|--------|
| **Added** | Root **`START_HERE.md`** | Obvious GitHub / `git clone` entry; leads with data pipeline + run commands |
| **Updated** | `README.md`, `CONTRIBUTING.md`, `HANDOFF.md`, `docs/README.md`, `docs/HOW_TO_RUN_EVERYTHING.md`, `pipeline/README.md`, `data/README.md`, `docs/DOC_LINK_CHECK.txt` | Point first-time readers at `START_HERE.md`; rename README section to “Data pipeline and benchmark” |

## 2026-04-21 — CPTAC/CCLE sample + annotation docs for data pipeline

| Action | Target | Reason |
|--------|--------|--------|
| **Updated** | `pipeline/psm_to_gene_matrix/README.md`, `START_HERE.md` §0, `data/manifests/EXPECTED_INPUTS.md`, `README.md` (clone + pipeline bullets) | Explicit **manifests, `sample_files_msstats_tmt.csv`, `.sample.txt`, mirror**, **CCLE peptide vs pre-built matrix**, vs **benchmark** `data/annotations/` tables |

## 2026-04-21 — Single doc: annotation from samples + simpler entry points

| Action | Target | Reason |
|--------|--------|--------|
| **Added** | **`docs/ANNOTATION_FROM_SAMPLES.md`** | One narrative for CPTAC **`sample.txt` → annotation** and **CCLE Sheet2 + converter → annotation**; same R driver with `--msstats_input_dir` |
| **Slimmed** | `START_HERE.md` §0, `pipeline/psm_to_gene_matrix/README.md` inputs section | Point to the new doc instead of duplicating tables |
| **Updated** | `README.md`, `HANDOFF.md`, `docs/README.md`, `docs/HOW_TO_RUN_EVERYTHING.md`, `docs/NAMING_AND_PATHS.md`, `pipeline/README.md`, `data/ccle_peptide/README_CCLE.md`, `docs/DOC_LINK_CHECK.txt` | Cross-links + “small doc spine” line on README |

## Future (optional)

- Refresh **`docs/REPO_AUDIT.md`** tree snapshot on a schedule or after refactors.
- If a paper freeze needs a Zenodo bundle, **copy** (not move) selected figures into `reports/final_report/exports/` and list them in `FIGURE_MANIFEST.md`.
