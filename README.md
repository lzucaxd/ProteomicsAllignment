# Proteomics alignment benchmark

**A calibrated benchmark** for evaluating **cross-dataset harmonization** of clinical (**CPTAC**) and preclinical (**CCLE**) TMT proteomics on a **shared gene space**, with explicit **null** and **ceiling** calibration.

> **Folder name:** the repository directory is spelled `ProteomicsAllignment` on disk (historical); documentation uses “alignment” in prose.

## The question

When we harmonize tumor and cell-line proteomes, do we improve **per-gene fold-change agreement** across domains, or mainly improve **PCA geometry**?

## The finding (illustrative numbers)

On **breast vs lung**, **raw** domain separation on PC1 is extreme (**domain R² ≈ 0.98**), while **Celligner** drives that toward mixing (**domain R² ≈ 0.003**) yet can **invert** cross-domain FC correlation for the same genes (see `comparison_summary.csv` and `disconnect_scores.csv`). **Bridge shift** preserves CPTAC–CCLE FC correlation for subtype while reducing domain dominance on PCs; see **`docs/METHODS.md`** for definitions.

## Quick start

```bash
# Python deps (repo root)
pip install -r requirements.txt
# Optional venv:
#   python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt

# R deps (limma, MSstatsTMT, data.table, ggplot2, …)
Rscript install_r_packages.R

# Harmonization benchmark (after gene_matrix.csv exist — see data/README.md)
export PYTHONPATH="${PWD}/src${PYTHONPATH:+:$PYTHONPATH}"
bash scripts/benchmark/run_overnight_v2.sh

# Canonical results
test -f reports/benchmark_master/benchmark_results/comparison_summary.csv && \
  head -5 reports/benchmark_master/benchmark_results/comparison_summary.csv
```

**Preprocessing (PDC PSM → gene matrices)** is documented in [`scripts/preprocessing/README.md`](scripts/preprocessing/README.md) and runs primarily from **`data/`** (historical layout).

## Documentation map

| Doc | Audience |
|-----|----------|
| [`scripts/preprocessing/README.md`](scripts/preprocessing/README.md) | PSM → `gene_matrix.csv` |
| [`scripts/benchmark/README.md`](scripts/benchmark/README.md) | Overnight benchmark steps |
| [`data/README.md`](data/README.md) | What lives under `data/`, provenance |
| [`docs/METHODS.md`](docs/METHODS.md) | Raw, bridge shift, bridge shift+scale, Celligner |
| [`docs/END_TO_END_TECHNICAL_REPORT.md`](docs/END_TO_END_TECHNICAL_REPORT.md) | Paper-style full narrative |
| [`docs/BENCHMARK_V2_AND_PRESENTATION.md`](docs/BENCHMARK_V2_AND_PRESENTATION.md) | Slides + paths checklist |
| [`docs/HANDOFF_CHECKLIST.md`](docs/HANDOFF_CHECKLIST.md) | Lab reproduction checklist |
| [`docs/REPO_AUDIT.md`](docs/REPO_AUDIT.md) | Machine-generated tree snapshot |
| [`PROJECT_REPORT.md`](PROJECT_REPORT.md) | Inventory-style overview |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | Conventions |

## Methods evaluated (benchmark)

| Repo ID | Plot label | Role |
|---------|--------------|------|
| `raw` | Raw | No correction |
| `bridge_shift` | Bridge shift | Per-gene offset from **bridge channel** summaries in `msstats_input.tsv` (Norm rows); **not** domain medians of sample abundances |
| `bridge_scale` | Bridge shift+scale | Bridge medians + **MAD**-based rescaling |
| `celligner` | Celligner | cPCA + MNN in PC space; **no per-gene preservation guarantee** |

## Repository layout (abbrev.)

```
configs/                 # YAML: preprocessing, tasks, methods
data/                    # Raw + MSstatsTMT outputs (mostly gitignored) + processed union
scripts/
  preprocessing/         # Documentation + future home for PSM scripts (see README)
  benchmark/             # run_overnight_v2.sh + metrics + calibration
  methods/               # Method drivers (Celligner, bridge-aware R, …)
  presentation/          # Figures / slide exports
src/harmonize/           # Python package (preprocessing + benchmark helpers)
reports/benchmark_master/# Results, diagnostics, meeting figures
docs/                    # Design + handoff + audit
archive/                 # Legacy / duplicate outputs (moved, not deleted)
```

## Optional: Celligner environment

See original setup in `pyproject.toml` / `install_r_packages.R` and:

```bash
pip install -e ".[celligner]"
# mnnpy build notes in CONTRIBUTING / PROJECT_REPORT
```

## Citation

Zamfira, L.-A. (2025–2026). *A calibrated benchmark for cross-dataset harmonization of clinical and preclinical cancer proteomics.* Vitek Lab, Northeastern University.

## License

See [`LICENSE`](LICENSE) (MIT).
