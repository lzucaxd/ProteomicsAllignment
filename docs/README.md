# Documentation index

## Start here

| Doc | Use when you need… |
|-----|---------------------|
| **[HOW_TO_RUN_EVERYTHING.md](HOW_TO_RUN_EVERYTHING.md)** | **Run the full pipeline:** data → matrices → benchmark → custom paths. |
| **[CLEAN_CLONE_REPRODUCIBILITY.md](CLEAN_CLONE_REPRODUCIBILITY.md)** | Fresh clone: verify script, data policy, what to commit for papers. |
| **[HANDOFF_CHECKLIST.md](HANDOFF_CHECKLIST.md)** | Lab pickup: reproduce, extend methods/tasks, path checks. |
| **[REPO_AUDIT.md](REPO_AUDIT.md)** | Machine-generated tree snapshot (refresh after major moves). |
| **[METHODS.md](METHODS.md)** | Harmonization methods (raw, bridge, Celligner) for papers. |
| **[PREPROCESSING.md](PREPROCESSING.md)** | Index → PSM / MSstatsTMT preprocessing docs. |
| **[BENCHMARK.md](BENCHMARK.md)** | Index → overnight benchmark docs. |
| **[LAB_ONBOARDING.md](LAB_ONBOARDING.md)** | Clone on any machine: `CPTAC_LOCAL_MIRROR`, `PROTEOMICS_ALIGNMENT_ROOT`, Python/R. |
| **[END_TO_END_TECHNICAL_REPORT.md](END_TO_END_TECHNICAL_REPORT.md)** | Paper / supplement: full pipeline, methods, metrics, figures, reproducibility. |
| **[BENCHMARK_V2_AND_PRESENTATION.md](BENCHMARK_V2_AND_PRESENTATION.md)** | Slides + v2 benchmark: steps, paths, checklist. |
| [how_to_run_end_to_end.md](how_to_run_end_to_end.md) | Short pointer → **HOW_TO_RUN_EVERYTHING.md**. |

Repo root **[README.md](../README.md)** — PDC setup, Celligner extras, layout.  
**[data/PIPELINE_README.md](../data/PIPELINE_README.md)** — PDC → MSstatsTMT → matrix.

## Topic notes (design detail)

| Doc | Topic |
|-----|--------|
| [benchmark_overview.md](benchmark_overview.md) | Benchmark philosophy and evaluation levels. |
| [benchmark_task_definitions.md](benchmark_task_definitions.md) | Task definitions. |
| [benchmark_subsampling_strategy.md](benchmark_subsampling_strategy.md) | Subsampling / cohort rules. |
| [preprocessing_overview.md](preprocessing_overview.md) | Shared-space preprocessing. |
| [config_system_overview.md](config_system_overview.md) | YAML config layout. |
| [representation_level_inference_overview.md](representation_level_inference_overview.md) | Limma / representation-level DA. |
| [native_domain_inference_overview.md](native_domain_inference_overview.md) | Native (MSstatsTMT / design-aware) inference. |
| [structure_metrics_overview.md](structure_metrics_overview.md) | PCA / structure metrics. |
| [matching_metrics_overview.md](matching_metrics_overview.md) | Matching metrics (if used). |
| [method_interface_overview.md](method_interface_overview.md) | Method interfaces. |
| [diagnostics_master_summary.md](diagnostics_master_summary.md) | Diagnostics scope. |
| [meeting_outputs_overview.md](meeting_outputs_overview.md) | Meeting figures: v2 paths vs legacy export. |

## Repository reports (outside `docs/`)

| Path | Role |
|------|------|
| [PROJECT_REPORT.md](../PROJECT_REPORT.md) | Long-form project inventory (PDC, PDC000120, CCLE, `data/scripts/`). |
| [reports/benchmark_master/benchmark_pipeline_overview.md](../reports/benchmark_master/benchmark_pipeline_overview.md) | Older R `run_benchmark()` layout; v2 shell orchestration is primary. |
| [reports/BENCHMARKING_AND_SLIDES_REPORT.md](../reports/BENCHMARKING_AND_SLIDES_REPORT.md) | Earlier subtype / slides context. |
