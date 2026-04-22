# Proteomics alignment benchmark

**A calibrated benchmark** for evaluating **cross-dataset harmonization** of clinical (**CPTAC**) and preclinical (**CCLE**) TMT proteomics on a **shared gene space**, with explicit **null** and **ceiling** calibration.

> **Folder name:** the repository directory is spelled `ProteomicsAllignment` on disk (historical); documentation uses “alignment” in prose.

> **Run the full pipeline (start here):** **[`docs/HOW_TO_RUN_EVERYTHING.md`](docs/HOW_TO_RUN_EVERYTHING.md)** — single guide for data, install, PSM → `gene_matrix.csv`, overnight benchmark, optional tables/slides, and **custom paths / your own matrices**.

## The question

When we harmonize tumor and cell-line proteomes, do we improve **per-gene fold-change agreement** across domains, or mainly improve **PCA geometry**?

## The finding (illustrative numbers)

On **breast vs lung**, **raw** domain separation on PC1 is extreme (**domain R² ≈ 0.98**), while **Celligner** drives that toward mixing (**domain R² ≈ 0.003**) yet can **invert** cross-domain FC correlation for the same genes (see `comparison_summary.csv` and `disconnect_scores.csv`). **Bridge shift** preserves CPTAC–CCLE FC correlation for subtype while reducing domain dominance on PCs; see **`docs/METHODS.md`** for definitions.

---

## Running everything (end-to-end)

**Full walkthrough (single page):** **[`docs/HOW_TO_RUN_EVERYTHING.md`](docs/HOW_TO_RUN_EVERYTHING.md)** — data acquisition, install, verify, PSM → `gene_matrix.csv`, **`run_overnight_v2.sh`**, optional tables/figures, and **pointing configs at your matrices**.

**Minimal path** (after manifests, `sample_files_msstats_tmt.csv`, and a venv with `pip install -r requirements.txt` + **`pip install -e .`**):

```bash
cd /your/actual/path/ProteomicsAllignment
source .venv/bin/activate
cd data && ./run_pipeline_per_manifest.sh    # or ./run_batch_studies.sh with CPTAC_LOCAL_MIRROR
cd .. && bash scripts/benchmark/run_overnight_v2.sh
```

Use your real clone path, not a placeholder like `/path/to/ProteomicsAllignment`. Repro checks and commit policy: **`docs/CLEAN_CLONE_REPRODUCIBILITY.md`**.

---

## Methods evaluated (benchmark)

| Repo ID | Plot label | Role |
|---------|--------------|------|
| `raw` | Raw | No correction |
| `bridge_shift` | Bridge shift | Per-gene offset from **bridge channel** summaries in `msstats_input.tsv` (Norm rows); **not** domain medians of sample abundances |
| `bridge_scale` | Bridge shift+scale | Bridge medians + **MAD**-based rescaling |
| `celligner` | Celligner | cPCA + MNN in PC space; **no per-gene preservation guarantee** |

---

## Repository layout (abbrev.)

```
configs/                 # YAML: preprocessing, tasks, methods
data/                    # Raw + MSstatsTMT outputs (mostly gitignored) + processed union
scripts/
  preprocessing/       # Docs index; executables still under data/ (see below)
  benchmark/             # run_overnight_v2.sh + metrics + calibration
  methods/               # Method drivers (Celligner, bridge-aware R, …)
  presentation/          # Figures / slide exports
src/harmonize/           # Python package (preprocessing + benchmark helpers)
reports/benchmark_master/# Results, diagnostics, meeting figures, final_tables/
docs/                    # Design + handoff + audit
```

**Preprocessing executables** (`run_pipeline_per_manifest.sh`, `pdc_manifest_downloader.py`, `pdc_psm_to_msstatsTMT_protein_matrix.R`, …) live under **`data/`** today; see **`scripts/preprocessing/README.md`** for the map.

---

## Documentation map

| Doc | Audience |
|-----|----------|
| **`docs/HOW_TO_RUN_EVERYTHING.md`** | **Central run guide:** clone → data → matrices → benchmark → custom paths |
| **`data/PIPELINE_README.md`** | PDC manifest → download → MSstatsTMT → `gene_matrix.csv` (authoritative) |
| **`data/manifests/README.md`** | Where to get PDC manifests (PSM / Text); what to save locally |
| **`docs/LAB_ONBOARDING.md`** | `CPTAC_LOCAL_MIRROR`, clone layout, env vars |
| **`scripts/preprocessing/README.md`** | Preprocessing narrative + path map |
| **`scripts/benchmark/README.md`** | Overnight benchmark steps 0–12 |
| **`data/README.md`** | What lives under `data/`, regeneration |
| **`docs/METHODS.md`** | Raw, bridge shift, bridge shift+scale, Celligner |
| **`docs/END_TO_END_TECHNICAL_REPORT.md`** | Paper-style full narrative |
| **`docs/BENCHMARK_V2_AND_PRESENTATION.md`** | Slides + paths checklist |
| **`docs/CLEAN_CLONE_REPRODUCIBILITY.md`** | Clone → install → verify → run (canonical) |
| **`docs/HANDOFF_CHECKLIST.md`** | Lab reproduction checklist |
| **`docs/REPO_AUDIT.md`** | Machine-generated tree snapshot |
| **`PROJECT_REPORT.md`** | Inventory-style overview |
| **`CONTRIBUTING.md`** | Conventions |

---

## What not to commit

Large downloads (`data/pdc_psm/`, `data/results/`, …), **dated PDC manifest CSVs** (expiring URLs), `.venv/`, and local mirrors are **gitignored** (see root `.gitignore`). The repo keeps **`data/manifests/README.md`** and **`example_pdc_file_manifest.csv`** only. Regenerated benchmark CSVs under `reports/benchmark_master/` are a **team policy** choice.

---

## Citation

Zamfira, L.-A. (2025–2026). *A calibrated benchmark for cross-dataset harmonization of clinical and preclinical cancer proteomics.* Vitek Lab, Northeastern University.

## License

See [`LICENSE`](LICENSE) (MIT).
