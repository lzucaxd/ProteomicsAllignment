# Shared paths for scripts/presentation/*.R
# Run from repo root: Rscript scripts/presentation/<script>.R

.local_ca <- commandArgs(trailingOnly = FALSE)
.local_fl <- .local_ca[startsWith(.local_ca, "--file=")]
if (length(.local_fl)) {
  PRES_SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", .local_fl[1L])))
} else {
  PRES_SCRIPT_DIR <- normalizePath(file.path(getwd(), "scripts", "presentation"), mustWork = FALSE)
}
.hp <- file.path(dirname(PRES_SCRIPT_DIR), "benchmark", "harmonize_paths.R")
if (!file.exists(.hp)) {
  .hp <- normalizePath(file.path(getwd(), "scripts", "benchmark", "harmonize_paths.R"), mustWork = FALSE)
}
if (!file.exists(.hp)) {
  stop("Cannot find harmonize_paths.R; run from repository root or set PROTEOMICS_ALIGNMENT_ROOT")
}
source(.hp)
REPO <- harmonize_repo_root()
PRES_OUT <- file.path(REPO, "presentation_materials")

pres_ensure_dirs <- function() {
  subdirs <- c("main_slides", "backup_slides", "tables", "figures", "checks",
               "figures/meeting", "figures/marker_profiles", "figures/structure",
               "figures/msstats_profiles", "figures/report")
  for (d in subdirs) {
    dir.create(file.path(PRES_OUT, d), recursive = TRUE, showWarnings = FALSE)
  }
}
