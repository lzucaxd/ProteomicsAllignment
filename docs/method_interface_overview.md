# Method Interface Overview

## Contract

Every harmonization method must produce a `MethodResult` containing:

| Field | Type | Description |
|-------|------|-------------|
| `matrix` | DataFrame (genes x samples) | The harmonized representation |
| `sample_meta` | DataFrame | sample_id, domain, condition, study_id |
| `feature_meta` | DataFrame | gene, included, exclusion_reason, bridge_corrected |
| `method_name` | str | Machine-readable identifier |
| `display_name` | str | Human-readable label |
| `notes` | str | Method-specific notes and caveats |
| `qc_paths` | list[str] | Paths to QC output files |
| `transforms_values` | bool | Whether the method changes the value scale |
| `value_scale` | str | "log2_abundance", "z_scored", "log2_abundance_calibrated" |

## Python Interface

```python
from harmonize.methods.base import MethodInterface, MethodResult

class MyMethod(MethodInterface):
    name = "my_method"
    display_name = "My Method"

    def run(self, cptac_matrix, ccle_matrix, sample_meta, config, outdir):
        # ... harmonization logic ...
        return MethodResult(matrix=..., sample_meta=..., ...)
```

## R Interface

Existing R methods implement the same contract via `scripts/methods/method_interface.R`:

```r
make_method_result(matrix, sample_meta, feature_meta, method_name, notes, qc_paths)
save_method_result(result, outdir)
```

Python wrappers call R scripts via subprocess and load the saved outputs.

## Registered Methods

| Name | Python class | R script | Transforms values? |
|------|-------------|----------|-------------------|
| `raw` | `RawMethod` | `run_raw_representation.R` | No |
| `bridge_aware` | `BridgeAwareMethod` | `bridge_aware_correction.R` | Yes (calibrated) |
| `celligner` | `CellignerMethod` | `run_celligner_all_data.py` | Yes (z-scored) |

## Adding a New Method

1. Create a Python class inheriting from `MethodInterface`
2. Register it in `src/harmonize/methods/registry.py`
3. Create a YAML config in `configs/methods/`
4. The benchmark will automatically discover and evaluate it
