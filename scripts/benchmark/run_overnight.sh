#!/usr/bin/env bash
# Overnight benchmark: Steps 0–11 (union matrices, per-method nulls, assembly, figures).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
export PYTHONPATH="${REPO}/src${PYTHONPATH:+:$PYTHONPATH}"
if [[ ! -f "${REPO}/renv/activate.R" ]]; then
  export R_PROFILE_USER="/dev/null"
fi

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${REPO}/reports/benchmark_master/logs"
DIAG_DIR="${REPO}/reports/benchmark_master/diagnostics"
mkdir -p "$LOG_DIR" "$DIAG_DIR"

LOG="${LOG_DIR}/overnight_${TIMESTAMP}.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== Benchmark overnight run: ${TIMESTAMP} ==="
echo "Repo: ${REPO}"
echo "Log: ${LOG}"
echo "Intersection coverage threshold: 0.30 (see preflight gene_coverage_audit_*.csv)"
echo "Random seed: 42 (permutation nulls, ceilings)"

die() { echo "FATAL: $*" >&2; exit 1; }

echo ""
echo "[Step 0] Pre-flight diagnostics"
Rscript "${REPO}/scripts/benchmark/preflight_diagnostics.R" --repo-root "$REPO" \
  --processed-dir data/processed_union || die "Step 0 failed"

echo ""
echo "[Step 1] Intersection gene masks"
Rscript "${REPO}/scripts/benchmark/compute_intersection_masks.R" --repo-root "$REPO" || die "Step 1 failed"

echo ""
echo "[Step 2] Regenerate method matrices from union"
# Set --skip-celligner if Celligner deps are unavailable (bridge + raw still run).
python3 "${REPO}/scripts/benchmark/regenerate_methods_union.py" --repo "$REPO" \
  --processed-dir "${REPO}/data/processed_union" \
  --out-methods-root "${REPO}/data/processed/methods" || die "Step 2 failed"

echo ""
echo "[Step 2b] Structure metrics (fixed-basis PCA from raw)"
python3 "${REPO}/scripts/benchmark/run_structure_batch.py" --repo "$REPO" || die "Step 2b failed"

echo ""
echo "[Step 3] limma DA (all methods × tasks)"
Rscript "${REPO}/scripts/benchmark/run_all_limma_da.R" --repo-root "$REPO" \
  --methods-root data/processed/methods \
  --meta-dir data/processed_union \
  --results-root reports/benchmark_master/benchmark_results || die "Step 3 failed"

echo ""
echo "[Step 4] Cross-domain metrics (union + intersection)"
Rscript "${REPO}/scripts/benchmark/compute_cross_domain_metrics.R" --repo-root "$REPO" \
  --intersection-dir data/processed_union || die "Step 4 failed"

echo ""
echo "[Step 5] Stratified FC"
Rscript "${REPO}/scripts/benchmark/compute_stratified_fc.R" --repo-root "$REPO" || die "Step 5 failed"

echo ""
echo "[Step 6] Volcano plots (raw)"
Rscript "${REPO}/scripts/benchmark/generate_volcanos.R" --repo-root "$REPO" || die "Step 6 failed"

echo ""
echo "[Step 7] Permutation nulls (~hours; raw null reused for bridge_shift)"
Rscript "${REPO}/scripts/benchmark/run_permutation_nulls.R" --repo-root "$REPO" \
  --n-perm 1000 --seed 42 || die "Step 7 failed"

echo ""
echo "[Step 8] Concordance ceilings (CPTAC)"
Rscript "${REPO}/scripts/benchmark/run_concordance_ceilings.R" --repo-root "$REPO" \
  --n-splits 200 --seed 42 || die "Step 8 failed"

echo ""
echo "[Step 9] Fast calibration"
Rscript "${REPO}/scripts/benchmark/run_fast_calibration.R" --repo-root "$REPO" || die "Step 9 failed"

echo ""
echo "[Step 10] Assemble comparison table"
python3 "${REPO}/scripts/benchmark/assemble_comparison_table.py" --repo "$REPO" || die "Step 10 failed"

echo ""
echo "[Step 11] Meeting figures"
Rscript "${REPO}/scripts/benchmark/generate_meeting_figures.R" --repo-root "$REPO" || die "Step 11 failed"

echo ""
echo "=== DONE $(date) ==="
echo "Diagnostics: ${DIAG_DIR}"
echo "Results: ${REPO}/reports/benchmark_master/benchmark_results/"
echo "Figures: ${REPO}/reports/benchmark_master/meeting/figures/"
echo "Log: ${LOG}"
