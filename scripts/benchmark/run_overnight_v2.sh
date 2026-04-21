#!/usr/bin/env bash
# Overnight benchmark v2: expanded CCLE subtype (25 lines, CAL120 merged) + BvL breast from v2 + split-half CCLE ceiling.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
if [[ ! -f "${REPO}/renv/activate.R" ]]; then
  export R_PROFILE_USER="/dev/null"
fi
export PYTHONPATH="${REPO}/src${PYTHONPATH:+:$PYTHONPATH}"
if [[ -x "${REPO}/.venv/bin/python3" ]]; then
  PYTHON="${REPO}/.venv/bin/python3"
else
  PYTHON="python3"
fi

UNION_DIR="data/processed/union"
META_DIR="${UNION_DIR}"
INTER_DIR="data/processed"
METHODS_ROOT="data/processed/methods"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${REPO}/reports/benchmark_master/logs"
DIAG_DIR="${REPO}/reports/benchmark_master/diagnostics"
mkdir -p "$LOG_DIR" "$DIAG_DIR" "${REPO}/data/processed"

LOG="${LOG_DIR}/overnight_v2_${TIMESTAMP}.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== Benchmark overnight v2: ${TIMESTAMP} ==="
echo "Union dir: ${UNION_DIR} | Intersection dir: ${INTER_DIR}"

die() { echo "FATAL: $*" >&2; exit 1; }

echo ""
echo "[Step 0] Annotation processing + validation"
Rscript "${REPO}/scripts/benchmark/process_ccle_annotations_v2.R" --repo-root "$REPO" || die "Step 0"

echo ""
echo "[Step 1] Build shared gene matrices -> ${UNION_DIR}"
"${PYTHON}" "${REPO}/scripts/run_preprocessing.py" --repo-root "$REPO" --output-dir "$UNION_DIR" || die "Step 1"

echo ""
echo "[Step 1b] Gene coverage audit + intersection lists"
Rscript "${REPO}/scripts/benchmark/preflight_diagnostics.R" --repo-root "$REPO" \
  --processed-dir "$UNION_DIR" || die "preflight"
Rscript "${REPO}/scripts/benchmark/compute_intersection_masks.R" --repo-root "$REPO" \
  --intersection-out-dir "$UNION_DIR" || die "intersection masks"

echo ""
echo "[Step 2] Regenerate method matrices"
"${PYTHON}" "${REPO}/scripts/benchmark/regenerate_methods_union.py" --repo "$REPO" \
  --processed-dir "${REPO}/${UNION_DIR}" \
  --out-methods-root "${REPO}/${METHODS_ROOT}" || die "Step 2"

echo ""
echo "[Step 2b] Structure metrics"
"${PYTHON}" "${REPO}/scripts/benchmark/run_structure_batch.py" --repo "$REPO" \
  --meta-dir "${REPO}/${META_DIR}" || die "Step 2b"

echo ""
echo "[Step 3] limma DA (16 runs)"
Rscript "${REPO}/scripts/benchmark/run_all_limma_da.R" --repo-root "$REPO" \
  --methods-root "$METHODS_ROOT" \
  --meta-dir "$META_DIR" \
  --results-root reports/benchmark_master/benchmark_results || die "Step 3"

echo ""
echo "[Step 4] Cross-domain metrics"
Rscript "${REPO}/scripts/benchmark/compute_cross_domain_metrics.R" --repo-root "$REPO" \
  --intersection-dir "$INTER_DIR" || die "Step 4"

echo ""
echo "[Step 5] Stratified FC"
Rscript "${REPO}/scripts/benchmark/compute_stratified_fc.R" --repo-root "$REPO" || die "Step 5"

echo ""
echo "[Step 6] Volcano plots"
Rscript "${REPO}/scripts/benchmark/generate_volcanos.R" --repo-root "$REPO" || die "Step 6"

echo ""
echo "[Step 7] Permutation nulls (slow)"
Rscript "${REPO}/scripts/benchmark/run_permutation_nulls.R" --repo-root "$REPO" \
  --meta-dir "$META_DIR" \
  --intersection-dir "$INTER_DIR" \
  --n-perm 1000 --seed 42 || die "Step 7"

echo ""
echo "[Step 8] Concordance ceilings (CPTAC + CCLE split-half)"
Rscript "${REPO}/scripts/benchmark/run_concordance_ceilings.R" --repo-root "$REPO" \
  --meta-dir "$META_DIR" \
  --intersection-dir "$INTER_DIR" \
  --ccle-split-half \
  --n-splits 200 --seed 42 || die "Step 8"

echo ""
echo "[Step 9] Fast calibration"
Rscript "${REPO}/scripts/benchmark/run_fast_calibration.R" --repo-root "$REPO" \
  --meta-dir "$META_DIR" || die "Step 9"

echo ""
echo "[Step 10] Disconnect scores"
Rscript "${REPO}/scripts/benchmark/compute_disconnect_scores.R" --repo-root "$REPO" || die "Step 10"

echo ""
echo "[Step 11] Assemble comparison table"
"${PYTHON}" "${REPO}/scripts/benchmark/assemble_comparison_table.py" --repo "$REPO" || die "Step 11"

echo ""
echo "[Step 12] Meeting figures"
Rscript "${REPO}/scripts/benchmark/generate_meeting_figures.R" --repo-root "$REPO" || die "Step 12"

echo ""
echo "=== DONE $(date) ==="
echo "Log: $LOG"
