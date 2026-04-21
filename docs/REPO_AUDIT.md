# Repository audit (generated for lab handoff)

Generated: 2026-04-21T02:42:12Z  
Repo root: /Users/zamfiraluca/Desktop/ProteomicsAllignment

**Handoff note:** duplicate slide/report PNGs previously under `presentation_materials/figures/report_figures/` were **moved** to `archive/report_figures_legacy/` (see `archive/README.md`). Canonical regenerated figures: `presentation_materials/figures/report/` via `scripts/presentation/generate_report_figures.R`.

## 1. Directory structure (maxdepth 3)
.
./.venv
./.venv-diagnostics
./.venv-diagnostics/bin
./.venv-diagnostics/include
./.venv-diagnostics/include/python3.13
./.venv-diagnostics/lib
./.venv-diagnostics/lib/python3.13
./.venv-diagnostics/share
./.venv-diagnostics/share/man
./.venv-venn
./.venv-venn/bin
./.venv-venn/include
./.venv-venn/include/python3.13
./.venv-venn/lib
./.venv-venn/lib/python3.13
./.venv-venn/share
./.venv-venn/share/man
./.venv/.mpl
./.venv/bin
./.venv/include
./.venv/include/python3.13
./.venv/include/site
./.venv/lib
./.venv/lib/python3.13
./.venv/share
./.venv/share/man
./configs
./configs/benchmark
./configs/methods
./configs/preprocessing
./configs/tasks
./data
./data/.venv
./data/.venv/bin
./data/.venv/include
./data/.venv/lib
./data/MSI VS MSS
./data/biospecimen
./data/ccle
./data/ccle_peptide
./data/ccle_sum
./data/cptac_samples
./data/cptac_samples/PDC000127
./data/manifests
./data/pdc_psm
./data/pdc_psm/PDC000120
./data/pdc_psm/PDC000153
./data/pdc_psm/PDC000234
./data/pdc_psm/PDC000270
./data/pdc_psm/PDC000325
./data/pdc_psm/PDC000327
./data/pdc_psm/PDC000329
./data/pdc_psm/PDC000446
./data/pdc_psm/PDC000489
./data/processed
./data/processed/methods
./data/processed/union
./data/processed_union
./data/results
./data/results/CCLE
./data/results/CCLE_corrected
./data/results/CCLE_qc_test
./data/results/PDC000120
./data/results/PDC000127
./data/results/PDC000153
./data/results/PDC000204
./data/results/PDC000221
./data/results/PDC000234
./data/results/PDC000270
./data/scripts
./docs
./models
./models/celligner-master
./models/celligner-master/R
./models/celligner-master/celligner
./models/celligner-master/celligner.egg-info
./models/celligner-master/docs
./models/celligner-master/man
./models/celligner-master/mnnpy
./models/celligner-master/runs
./presentation_materials
./presentation_materials/backup_slides
./presentation_materials/checks
./presentation_materials/figures
./presentation_materials/figures/assumptions
./presentation_materials/figures/main
./presentation_materials/figures/marker_agreement
./presentation_materials/figures/marker_profiles
./presentation_materials/figures/meeting
./presentation_materials/figures/msstats_profiles
./presentation_materials/figures/report
./presentation_materials/figures/report_figures
./presentation_materials/figures/structure
./presentation_materials/main_slides
./presentation_materials/tables
./reports
./reports/benchmark_master
./reports/benchmark_master/benchmark_results
./reports/benchmark_master/benchmark_results_union
./reports/benchmark_master/calibration
./reports/benchmark_master/cancer_type_structure
./reports/benchmark_master/celligner_all
./reports/benchmark_master/celligner_full
./reports/benchmark_master/celligner_task
./reports/benchmark_master/celligner_trial
./reports/benchmark_master/diagnostics
./reports/benchmark_master/logs
./reports/benchmark_master/marker_profiles
./reports/benchmark_master/meeting
./reports/benchmark_master/methods
./reports/benchmark_master/native_domain_da
./reports/benchmark_master/representation_level_da
./reports/benchmark_master/tasks
./reports/benchmark_v1
./reports/benchmark_v1/diagnostics
./reports/benchmark_v1/diagnostics_feedback
./reports/benchmark_v1/report_update
./reports/bridge_qc
./reports/figures
./reports/presentation_subtype_benchmark
./reports/presentation_subtype_benchmark/00_overview
./reports/presentation_subtype_benchmark/01_methods_and_assumptions
./reports/presentation_subtype_benchmark/02_cptac_subtype
./reports/presentation_subtype_benchmark/03_ccle_subtype
./reports/presentation_subtype_benchmark/04_qc_support
./reports/presentation_subtype_benchmark/05_overlap_and_summary
./reports/presentation_subtype_benchmark/06_alignment_benchmark
./reports/presentation_subtype_benchmark/07_meeting_notes
./results
./results/PDC000120
./results/PDC000120/diagnostics
./scripts
./scripts/benchmark
./scripts/methods
./scripts/presentation
./src
./src/harmonize
./src/harmonize/benchmark
./src/harmonize/methods
./src/harmonize/preprocessing
./src/harmonize/reporting
./src/harmonize/utils

## 2. Scripts under scripts/
scripts/benchmark/assemble_comparison_table.py
scripts/benchmark/benchmark_comparison_summary.R
scripts/benchmark/benchmark_runner.R
scripts/benchmark/bridge_aware_correction.R
scripts/benchmark/build_shared_union_v2.sh
scripts/benchmark/calibration_figures.R
scripts/benchmark/ccle_plex_batch_correction.R
scripts/benchmark/celligner_marker_and_agreement_analysis.R
scripts/benchmark/celligner_union_task.py
scripts/benchmark/compute_cross_domain_metrics.R
scripts/benchmark/compute_disconnect_scores.R
scripts/benchmark/compute_intersection_masks.R
scripts/benchmark/compute_stratified_fc.R
scripts/benchmark/diagnose_ccle_tissue_structure.py
scripts/benchmark/diagnose_subtype_sign.R
scripts/benchmark/diagnostics.R
scripts/benchmark/evaluation_helpers.R
scripts/benchmark/extract_bridge_summaries.R
scripts/benchmark/generate_meeting_figures.R
scripts/benchmark/generate_volcanos.R
scripts/benchmark/harmonize_paths.R
scripts/benchmark/install_celligner_benchmark.sh
scripts/benchmark/native_domain_da.R
scripts/benchmark/polished_profile_plots.R
scripts/benchmark/preflight_diagnostics.R
scripts/benchmark/process_ccle_annotations_v2.R
scripts/benchmark/regenerate_methods_union.py
scripts/benchmark/run_all_benchmarks.R
scripts/benchmark/run_all_limma_da.R
scripts/benchmark/run_celligner_all_data.py
scripts/benchmark/run_celligner_full.py
scripts/benchmark/run_celligner_task.py
scripts/benchmark/run_concordance_ceilings.R
scripts/benchmark/run_fast_calibration.R
scripts/benchmark/run_full_benchmark.R
scripts/benchmark/run_overnight.sh
scripts/benchmark/run_overnight_v2.sh
scripts/benchmark/run_permutation_nulls.R
scripts/benchmark/run_polished_profile_plots.R
scripts/benchmark/run_sample_profile_plots.R
scripts/benchmark/run_structure_batch.py
scripts/benchmark/sample_profile_plots.R
scripts/benchmark/subset_strategies.R
scripts/benchmark/task_breast_subtype.R
scripts/benchmark/task_breast_vs_lung.R
scripts/benchmark/trial_celligner_subsample.py
scripts/methods/method_interface.R
scripts/methods/run_bridge_aware_representation.R
scripts/methods/run_celligner_representation.R
scripts/methods/run_celligner_representation.py
scripts/methods/run_raw_representation.R
scripts/presentation/assumption_checks.R
scripts/presentation/bridge_analysis.R
scripts/presentation/celligner_check.R
scripts/presentation/check_all_assumptions.R
scripts/presentation/extract_all_numbers.R
scripts/presentation/extract_everything.R
scripts/presentation/extract_marker_panel.R
scripts/presentation/extract_slide_numbers.R
scripts/presentation/fc_se_summary.R
scripts/presentation/final_checks.R
scripts/presentation/generate_backup_doc.R
scripts/presentation/generate_profile_plots.R
scripts/presentation/generate_report_figures.R
scripts/presentation/plot_marker_agreement_by_method.R
scripts/presentation/plot_msstats_yeast_profiles.R
scripts/presentation/prepare_all.sh
scripts/presentation/presentation_paths.R
scripts/run_all.py
scripts/run_benchmark.py
scripts/run_meeting_exports.py
scripts/run_methods.py
scripts/run_native_baselines.py
scripts/run_preprocessing.py

## 3. Configs
configs/benchmark/default.yaml
configs/methods/bridge_aware.yaml
configs/methods/celligner.yaml
configs/methods/raw.yaml
configs/preprocessing/default.yaml
configs/preprocessing/union.yaml
configs/tasks/breast_subtype.yaml
configs/tasks/breast_vs_lung.yaml

## 4. Sample CSV under data/processed/ (first 40)
data/processed/methods/bridge_scale/transformed_breast_subtype.csv
data/processed/methods/bridge_scale/transformed_breast_vs_lung.csv
data/processed/methods/bridge_shift/transformed_breast_subtype.csv
data/processed/methods/bridge_shift/transformed_breast_vs_lung.csv
data/processed/methods/celligner/transformed_breast_subtype.csv
data/processed/methods/celligner/transformed_breast_vs_lung.csv
data/processed/methods/raw/transformed_breast_subtype.csv
data/processed/methods/raw/transformed_breast_vs_lung.csv
data/processed/feature_meta_breast_vs_lung.csv
data/processed/shared_gene_matrix_breast_vs_lung.csv
data/processed/sample_meta_breast_subtype.csv
data/processed/ccle_breast_subtype_annotation_processed.csv
data/processed/union/feature_meta_breast_vs_lung.csv
data/processed/union/shared_gene_matrix_breast_vs_lung.csv
data/processed/union/sample_meta_breast_subtype.csv
data/processed/union/feature_meta_breast_subtype.csv
data/processed/union/shared_gene_matrix_breast_subtype.csv
data/processed/union/sample_meta_breast_vs_lung.csv
data/processed/feature_meta_breast_subtype.csv
data/processed/shared_gene_matrix_breast_subtype.csv
data/processed/sample_meta_breast_vs_lung.csv

## 5. Sample CSV under reports/ (first 40)
reports/subtype_marker_sanity.csv
reports/subtype_marker_sanity_extended.csv
reports/benchmark_v1/shared_feature_table.csv
reports/benchmark_v1/marker_panel_master.csv
reports/subtype_canonical_markers_in_DA.csv
reports/subtype_marker_panel.csv
reports/presentation_subtype_benchmark/06_alignment_benchmark/alignment_tasks_table.csv
reports/presentation_subtype_benchmark/02_cptac_subtype/DA_MSstatsTMT_Luminal_vs_Basal_marker_sanity.csv
reports/presentation_subtype_benchmark/03_ccle_subtype/canonical_markers_slide_table_10genes.csv
reports/presentation_subtype_benchmark/03_ccle_subtype/canonical_markers_check.csv
reports/presentation_subtype_benchmark/03_ccle_subtype/canonical_markers_benchmark_by_gene.csv
reports/presentation_subtype_benchmark/03_ccle_subtype/canonical_markers_cptac_11protein_rows.csv
reports/presentation_subtype_benchmark/05_overlap_and_summary/da_summary_three_way_reference.csv
reports/presentation_subtype_benchmark/05_overlap_and_summary/subtype_marker_panel.csv
reports/presentation_subtype_benchmark/05_overlap_and_summary/subtype_cptac_ccle_summary.csv
reports/presentation_subtype_benchmark/05_overlap_and_summary/da_summary_table_primary_CCLE_corrected.csv
reports/subtype_cptac_ccle_summary.csv
reports/subtype_summary_table.csv
reports/benchmark_master/benchmark_results/comparison_summary_tiered.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/residual_corr_matrix_ccle.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/destruction_summary_bridge_scale.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/biology_destruction_bridge_scale.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/null_distribution.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/concordance_ceiling_ccle.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/marker_sanity_summary_bridge_scale_cptac_breast_vs_lung.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/concordance_ceiling_cptac.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/ceiling_summary_cptac.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/observed_metrics.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/ceiling_summary.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/marker_sanity_bridge_scale_cptac_breast_vs_lung.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/residual_dependence_ccle.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/residual_corr_matrix_cptac.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/observed_vs_null_summary.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/marker_sanity_summary_bridge_scale_ccle_breast_vs_lung.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/marker_sanity_bridge_scale_ccle_breast_vs_lung.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/residual_dependence_cptac.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/ceiling_summary_ccle.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/calibration/destruction_grid_bridge_scale.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/representation_da/ccle/da_limma_result.csv
reports/benchmark_master/benchmark_results/bridge_scale/breast_vs_lung/representation_da/cross_domain_metrics.csv

## 6. Files > 50MB (excluding .git; first 50)
./.venv/lib/python3.13/site-packages/llvmlite/binding/libllvmlite.dylib
./data/processed_union/shared_gene_matrix_breast_vs_lung.csv
./data/ccle_peptide/ccle_protein_quant_with_peptides_14745.tsv
./data/ccle_peptide/Table_S2_Protein_Quant_Normalized.xlsx
./data/ccle_sum/Table_S2_Protein_Quant_Normalized (1).xlsx
./data/results/CCLE/protein_summary.tsv
./data/results/CCLE/msstats_input.tsv
./data/results/PDC000234/msstats_input.tsv
./data/results/PDC000234/parsed_psm_long.tsv
./data/results/PDC000120/protein_summary.tsv
./data/results/PDC000120/msstats_input.tsv
./data/results/PDC000120/parsed_psm_long.tsv
./data/results/CCLE_qc_test/msstats_input.tsv
./data/results/CCLE_corrected/protein_summary.tsv
./data/results/CCLE_corrected/protein_matrix_wide.csv
./data/results/CCLE_corrected/msstats_input.tsv
./data/results/CCLE_corrected/gene_matrix.csv
./data/results/PDC000153/protein_summary.tsv
./data/results/PDC000153/protein_matrix_wide.csv
./data/results/PDC000153/msstats_input.tsv
./data/results/PDC000153/gene_matrix.csv
./data/results/PDC000153/parsed_psm_long.tsv
./data/results/PDC000270/parsed_psm_long.tsv
./data/processed/methods/bridge_scale/transformed_breast_vs_lung.csv
./data/processed/methods/bridge_shift/transformed_breast_vs_lung.csv
./data/processed/methods/celligner/transformed_breast_vs_lung.csv
./data/processed/methods/raw/transformed_breast_vs_lung.csv
./data/processed/shared_gene_matrix_breast_vs_lung.csv
./data/processed/union/shared_gene_matrix_breast_vs_lung.csv
./reports/benchmark_master/methods/bridge_aware/bridge_aware_shift_scale_matrix.csv
./reports/benchmark_master/methods/bridge_aware/bridge_aware_shift_only_matrix.csv
./reports/benchmark_master/celligner_all/celligner_aligned_matrix.csv
./reports/benchmark_master/celligner_full/celligner_aligned_matrix.csv
./reports/benchmark_master/celligner_task/breast_vs_lung/celligner_aligned_matrix.csv

## 7. Root .gitignore
```
# OS
.DS_Store

# Python
.venv/
.venv-diagnostics/
.venv-venn/
__pycache__/
*.py[cod]
*$py.class
.Python
*.egg-info/
.eggs/

# Large/generated data (keep code only)
data/pdc_psm/
data/results/

# Presentation bundle (regenerate: ./scripts/presentation/prepare_all.sh)
presentation_materials/

# Logs
*.log
MSstats_groupComparison_log_*.log
MSstats_dataProcess_log_*.log

# R default plot output (often created at cwd)
Rplots.pdf

# Tooling / local env
.pytest_cache/
.mypy_cache/
.ruff_cache/
.ipynb_checkpoints/

# Vendored Celligner: local builds (regenerate with pip install ./models/celligner-master/mnnpy)
models/celligner-master/mnnpy/build/
models/celligner-master/**/__pycache__/

# Optional: uncomment if you never commit generated benchmark artifacts
# data/processed/
# reports/benchmark_master/benchmark_results/
# reports/benchmark_master/logs/```

## Notes
- Audit command pipeline may omit tail if `head` closes early; re-run `find . -size +50M` locally for full large-file list.
