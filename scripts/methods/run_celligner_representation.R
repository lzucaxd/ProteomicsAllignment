#!/usr/bin/env Rscript
# =============================================================================
# Method 2 — Celligner-Aligned Representation (R wrapper)
# =============================================================================
# Calls the Python Celligner wrapper via system() for integration with the
# R-based benchmark runner. Returns the same interface as other methods.
#
# Usage:
#   source("scripts/methods/method_interface.R")
#   source("scripts/methods/run_celligner_representation.R")
#   result <- run_celligner_representation(
#     cptac_mat_path, ccle_mat_path, cptac_meta, ccle_meta, outdir
#   )
# =============================================================================

run_celligner_representation <- function(cptac_mat_path, ccle_mat_path,
                                          cptac_meta = NULL, ccle_meta = NULL,
                                          outdir = "reports/benchmark_master/methods/celligner_style",
                                          python_cmd = NULL,
                                          min_obs_frac = 0.5,
                                          impute = "median",
                                          compute_umap = TRUE) {

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(python_cmd)) {
    python_cmd <- if (file.exists("data/.venv/bin/python")) "data/.venv/bin/python" else "python3"
  }

  script <- file.path("scripts", "methods", "run_celligner_representation.py")
  if (!file.exists(script)) stop("Celligner Python wrapper not found: ", script)

  cmd_parts <- c(
    shQuote(python_cmd), shQuote(script),
    "--cptac_matrix", shQuote(cptac_mat_path),
    "--ccle_matrix", shQuote(ccle_mat_path),
    "--outdir", shQuote(outdir),
    "--min_obs_frac", min_obs_frac,
    "--impute", impute
  )
  if (!compute_umap) cmd_parts <- c(cmd_parts, "--no_umap")

  # Save metadata for Python to read
  if (!is.null(cptac_meta) && is.data.frame(cptac_meta)) {
    meta_path <- file.path(outdir, ".cptac_meta_tmp.csv")
    fwrite(as.data.table(cptac_meta), meta_path)
    cmd_parts <- c(cmd_parts, "--cptac_meta", shQuote(meta_path))
  }
  if (!is.null(ccle_meta) && is.data.frame(ccle_meta)) {
    meta_path <- file.path(outdir, ".ccle_meta_tmp.csv")
    fwrite(as.data.table(ccle_meta), meta_path)
    cmd_parts <- c(cmd_parts, "--ccle_meta", shQuote(meta_path))
  }

  cmd <- paste(cmd_parts, collapse = " ")
  message("Running Celligner Python wrapper:\n  ", cmd)
  ret <- system(cmd, intern = FALSE)

  if (ret != 0) {
    warning("Celligner Python wrapper exited with code ", ret,
            ". Check ", file.path(outdir, "method_notes.txt"), " for details.")
  }

  # Load results back into R
  mat_path <- file.path(outdir, "transformed_matrix.csv")
  if (!file.exists(mat_path)) {
    stop("Celligner did not produce transformed_matrix.csv in ", outdir)
  }
  mat_dt <- fread(mat_path, header = TRUE)
  gene_col <- names(mat_dt)[1]
  sample_cols <- setdiff(names(mat_dt), gene_col)
  mat <- as.matrix(mat_dt[, ..sample_cols])
  rownames(mat) <- mat_dt[[gene_col]]

  sample_meta <- fread(file.path(outdir, "sample_metadata.csv"))
  feature_meta <- fread(file.path(outdir, "feature_metadata.csv"))
  notes <- readLines(file.path(outdir, "method_notes.txt"))

  # Clean up temp files
  unlink(file.path(outdir, ".cptac_meta_tmp.csv"))
  unlink(file.path(outdir, ".ccle_meta_tmp.csv"))

  make_method_result(
    matrix       = mat,
    sample_meta  = as.data.frame(sample_meta),
    feature_meta = as.data.frame(feature_meta),
    method_name  = "celligner",
    method_notes = notes,
    qc_paths     = c(
      transformed_matrix = mat_path,
      notes = file.path(outdir, "method_notes.txt"),
      umap = file.path(outdir, "celligner_umap.csv"),
      de_genes = file.path(outdir, "celligner_de_genes.txt")
    )
  )
}
