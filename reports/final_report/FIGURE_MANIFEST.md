# Figure manifest (handoff)

Mapping from **common report / slide needs** to **canonical generator + output location**. Paths are relative to the **repository root**. Many figure trees are **gitignored**; regenerate on a machine that has run the benchmark or presentation scripts.

## Benchmark / PCA-style panels

| Concept | Typical outputs | How to rebuild |
|---------|-----------------|----------------|
| PCA panels (breast subtype, BvL, …) | `presentation_materials/figures/report/pca_panel_*.png` | `Rscript scripts/presentation/generate_report_figures.R` |
| MSstats / profile plots | `presentation_materials/figures/msstats_profiles/*.png` | `scripts/presentation/generate_profile_plots.R` (see script header) |

## Benchmark diagnostics (under `reports/benchmark_master/`)

| Concept | Location | Notes |
|---------|----------|--------|
| Overnight logs | `reports/benchmark_master/logs/overnight_v2_*.log` | gitignored |
| Diagnostics folder | `reports/benchmark_master/diagnostics/` | QC from pipeline steps |
| Meeting / slide exports | `reports/benchmark_master/meeting/` | From `generate_meeting_figures.R` |

## Subtype / marker exploratory plots (legacy / lab)

Some older PNGs and gene lists live under **`reports/`** root (e.g. `venn_*.png`, `subtype_*.csv`). Treat as **exploratory** unless your paper explicitly cites them; prefer regenerating from `scripts/presentation/` and `scripts/benchmark/` for a clean provenance chain.

## LaTeX / paper

If the paper uses generated `.tex` from `final_tables/`, list the exact filenames in your supplement and keep the generating R script pinned to a **git tag**.

## Vendor / Celligner

Notebooks under `models/celligner-master/` are **upstream**; not part of this repo’s figure manifest unless you opt into tracking them.
