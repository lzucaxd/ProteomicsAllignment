#!/usr/bin/env bash
# Assemble presentation_materials/ (numbers, figures, backup markdown).
# Run from repository root: ./scripts/presentation/prepare_all.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

echo "=== Preparing presentation materials in $REPO_ROOT/presentation_materials ==="
mkdir -p presentation_materials/{main_slides,backup_slides,tables,figures,checks}

RSCRIPT=(Rscript)
if [[ -x "$REPO_ROOT/.venv/bin/python" ]]; then :; fi

echo "[1/10] Extracting slide numbers..."
set +e
Rscript scripts/presentation/extract_slide_numbers.R \
  > presentation_materials/checks/slide_numbers.txt 2>&1
set -e
if ! grep -q "SLIDE 3" presentation_materials/checks/slide_numbers.txt 2>/dev/null; then
  echo "(extract_slide_numbers.R may have failed — check comparison_summary.csv)" >> presentation_materials/checks/slide_numbers.txt
fi
cat presentation_materials/checks/slide_numbers.txt

echo "[2/10] Extracting marker panel..."
Rscript scripts/presentation/extract_marker_panel.R \
  >> presentation_materials/checks/slide_numbers.txt 2>&1 || true

echo "[3/10] Generating profile plots (requires cowplot)..."
Rscript scripts/presentation/generate_profile_plots.R 2>&1 | tail -20

echo "[4/10] Bridge analysis..."
Rscript scripts/presentation/bridge_analysis.R \
  >> presentation_materials/checks/slide_numbers.txt 2>&1 || true

echo "[5/10] Assumption checks (QQ, residuals, SA)..."
Rscript scripts/presentation/assumption_checks.R 2>&1 | tail -15

echo "[6/10] FC and SE summary..."
Rscript scripts/presentation/fc_se_summary.R \
  >> presentation_materials/checks/slide_numbers.txt 2>&1 || true

echo "[7/10] Celligner status check..."
Rscript scripts/presentation/celligner_check.R \
  >> presentation_materials/checks/slide_numbers.txt 2>&1 || true

echo "[7.5/10] Marker agreement figures (per method, CPTAC vs CCLE)..."
mkdir -p presentation_materials/figures/marker_agreement
Rscript scripts/presentation/plot_marker_agreement_by_method.R --repo-root "$REPO_ROOT"

echo "[8/10] Copying meeting + marker + structure figures..."
mkdir -p presentation_materials/figures/meeting
mkdir -p presentation_materials/figures/marker_profiles
mkdir -p presentation_materials/figures/structure
mkdir -p presentation_materials/figures/marker_agreement
if [[ -d reports/benchmark_master/meeting/figures ]]; then
  cp -R reports/benchmark_master/meeting/figures/. presentation_materials/figures/meeting/ 2>/dev/null || true
fi
if [[ -d reports/benchmark_master/marker_profiles ]]; then
  cp -R reports/benchmark_master/marker_profiles/. presentation_materials/figures/marker_profiles/ 2>/dev/null || true
fi
if [[ -d presentation_materials/figures/marker_agreement/breast_subtype ]]; then :; fi
for method in raw bridge_shift bridge_scale celligner; do
  for task in breast_subtype breast_vs_lung; do
    src="reports/benchmark_master/benchmark_results/${method}/${task}/structure"
    if [[ -d "$src" ]]; then
      mkdir -p "presentation_materials/figures/structure/${method}/${task}"
      cp "$src"/*.pdf "$src"/*.png presentation_materials/figures/structure/"${method}/${task}"/ 2>/dev/null || true
    fi
  done
done

echo "[9/10] Generating backup_content.md..."
Rscript scripts/presentation/generate_backup_doc.R 2>&1

echo "[10/10] Final checklist..."
Rscript scripts/presentation/final_checks.R

echo ""
echo "=== Done. Output: $REPO_ROOT/presentation_materials/ ==="
find presentation_materials -type f 2>/dev/null | sort | head -60
echo "..."
echo "Open presentation_materials/backup_slides/backup_content.md for backup Q&A."
