#!/usr/bin/env bash
# Install Python + R pieces needed for models/celligner-master (Celligner) in the benchmark venv.
# Run from anywhere: ./scripts/benchmark/install_celligner_benchmark.sh
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"

if [[ -x "${REPO}/.venv/bin/pip" ]]; then
  PIP="${REPO}/.venv/bin/pip"
  PY="${REPO}/.venv/bin/python3"
else
  PIP="${PIP:-pip3}"
  PY="${PYTHON:-python3}"
fi

echo "== Celligner benchmark deps =="
echo "pip: $PIP"
echo "python: $PY"
"$PY" -V

if "$PY" -c 'import sys; raise SystemExit(0 if sys.version_info < (3, 13) else 1)'; then
  :
else
  echo "WARNING: Python 3.13+ often lacks wheels for scanpy/rpy2. If installs fail, recreate .venv with Python 3.11 or 3.12."
fi

"$PIP" install --upgrade pip setuptools wheel
"$PIP" install Cython "numpy>=1.23" "scipy>=1.9" "scikit-learn>=1.1" "pandas>=1.5" "matplotlib>=3.6"

echo "== mnnpy (from Celligner vendor tree) =="
"$PIP" install "${REPO}/models/celligner-master/mnnpy"

echo "== rpy2 + scanpy stack =="
"$PIP" install "rpy2>=3.5.4" "scanpy>=1.9" "anndata>=0.8" "umap-learn>=0.5" python-igraph || {
  echo "FAILED: try Python 3.11/3.12 venv or install binary deps for igraph (brew install igraph)."
  exit 1
}

echo "== Verify imports =="
PYTHONPATH="${REPO}/models/celligner-master${PYTHONPATH:+:$PYTHONPATH}" "$PY" -c "
import sys
sys.path.insert(0, '${REPO}/models/celligner-master/mnnpy')
sys.path.insert(0, '${REPO}/models/celligner-master')
import mnnpy
from celligner import Celligner
print('Celligner OK:', Celligner)
"

echo "== R limma (optional but required at runtime for DE steps) =="
if command -v Rscript &>/dev/null; then
  Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install("limma", ask=FALSE, update=FALSE)'
  echo "limma installed via BiocManager."
else
  echo "Rscript not in PATH. Install R, then in R: BiocManager::install(\"limma\")"
fi

echo "Done. Re-run method matrices: PYTHONPATH=src ${PY} scripts/benchmark/regenerate_methods_union.py (or full overnight v2)."
