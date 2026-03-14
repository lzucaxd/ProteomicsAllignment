#!/bin/bash
# Run PDC PSM -> protein matrix (avoids .Rprofile/renv "cannot open renv/activate.R" error).
# Usage: ./run_pdc_psm_matrix.sh
#        ./run_pdc_psm_matrix.sh --psm_dir pdc_psm --out_matrix my_matrix.csv

cd "$(dirname "$0")"
Rscript --no-init-file pdc_psm_to_protein_matrix.R "$@"
