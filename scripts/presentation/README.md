# Presentation materials assembly

Builds **`presentation_materials/`** at the repo root: slide-ready numbers, copied figures, profile/QC PDFs, and **`backup_slides/backup_content.md`** for Q&A.

## Prerequisites

- Repo root, with **`reports/benchmark_master/benchmark_results/comparison_summary.csv`** (run `scripts/benchmark/run_overnight_v2.sh` first).
- **`data/processed/`** union + method matrices for profile and bridge scripts.
- R packages: `data.table`, `ggplot2`, `tidyr`, **`cowplot`** (see `install_r_packages.R`).

## Run everything

```bash
cd /path/to/ProteomicsAllignment
./scripts/presentation/prepare_all.sh
```

## Run steps individually

```bash
Rscript scripts/presentation/extract_slide_numbers.R    # → checks/slide_numbers.txt + tables/
Rscript scripts/presentation/extract_marker_panel.R
Rscript scripts/presentation/generate_profile_plots.R   # needs cowplot
Rscript scripts/presentation/bridge_analysis.R
Rscript scripts/presentation/assumption_checks.R
Rscript scripts/presentation/fc_se_summary.R
Rscript scripts/presentation/celligner_check.R
Rscript scripts/presentation/plot_marker_agreement_by_method.R   # → figures/marker_agreement/*.pdf (needs representation_level_da marker_summary)
Rscript scripts/presentation/check_all_assumptions.R             # → figures/assumptions/*.pdf + tables/assumption_summary.csv (needs cowplot + benchmark_results DA)
Rscript scripts/presentation/generate_backup_doc.R
Rscript scripts/presentation/final_checks.R
```

Paths resolve via **`scripts/benchmark/harmonize_paths.R`** (set **`PROTEOMICS_ALIGNMENT_ROOT`** if needed).

## Output layout

| Path | Contents |
|------|-----------|
| `presentation_materials/checks/slide_numbers.txt` | Printed metrics + logs from several steps |
| `presentation_materials/tables/` | `comparison_summary_full.csv`, marker panel, FC/SE, bridge tables |
| `presentation_materials/figures/` | Profiles, QQ, residuals, SA plots, `meeting/`, `structure/` |
| `presentation_materials/backup_slides/backup_content.md` | Consolidated backup tables |

Add `presentation_materials/` to `.gitignore` locally if you do not want the bundle in git.
