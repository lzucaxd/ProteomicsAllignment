# Config System Overview

## Design

All parameters are in YAML files under `configs/`. Python loads them
via `harmonize.utils.config`. No hardcoded paths in pipeline code.

## Config Structure

```
configs/
  preprocessing/default.yaml    Data sources, filtering thresholds
  tasks/breast_subtype.yaml     Task A definition
  tasks/breast_vs_lung.yaml     Task B definition
  methods/raw.yaml              Raw method parameters
  methods/bridge_aware.yaml     Bridge-aware parameters
  methods/celligner.yaml        Celligner parameters
  benchmark/default.yaml        Which tasks, methods, levels to run
```

## Key Config Fields

### preprocessing/default.yaml
- `data_sources`: Paths to CPTAC studies and CCLE
- `filtering`: `min_prevalence`, `min_sd`, `min_obs_frac`
- `imputation`: `strategy` (within_domain_gene_median)
- `output_dir`: Where processed matrices go

### tasks/*.yaml
- `task_name`, `contrast`, `block_order`
- `markers`, `expected_marker_directions`
- `cptac`: Study IDs, subset strategy, subtype mapping
- `ccle`: Cell line lists or tissue filters
- `native_domain_inference`, `representation_level_inference`

### methods/*.yaml
- `method_name`, `display_name`
- `feature_level`, `transforms_values`, `value_scale`
- `parameters`: Method-specific (min_obs_frac, mode, estimators, etc.)
- `r_script` or `python_script`: Path to implementation
- `output_dir`

### benchmark/default.yaml
- `tasks`: List of task names to evaluate
- `methods`: List of method names to evaluate
- `levels`: Toggle each evaluation level on/off
- Output directory paths for each level

## Usage in Code

```python
from harmonize.utils.config import load_config, load_task_config, load_method_config

bench_cfg = load_config("configs/benchmark/default.yaml")
task_cfg = load_task_config("breast_subtype")
method_cfg = load_method_config("bridge_aware")
```

## Customization

To run a custom benchmark:
1. Edit or copy the relevant YAML files
2. Pass `--config <your_config.yaml>` to the runner scripts
