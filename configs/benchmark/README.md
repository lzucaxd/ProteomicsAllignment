# Benchmark configuration

**Canonical file:** **`default.yaml`**

This directory holds YAML consumed by Python helpers (`scripts/run_benchmark.py`, `harmonize` benchmark code). The **overnight shell benchmark** (`scripts/benchmark/run_overnight_v2.sh`) orchestrates R + Python steps directly; keep task/method lists in sync with that script when you add new methods.

| File | Role |
|------|------|
| `default.yaml` | Tasks, methods, evaluation levels, structure/matching metric knobs |

To run with a custom copy:

```bash
python scripts/run_benchmark.py --config configs/benchmark/your_copy.yaml
```

For the **primary lab workflow**, prefer **`bash scripts/run_benchmark.sh`** (full overnight v2).
