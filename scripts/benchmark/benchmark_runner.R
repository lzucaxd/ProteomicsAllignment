#!/usr/bin/env Rscript
# =============================================================================
# Unified Benchmark Runner
# =============================================================================
#
# run_benchmark(matrix, sample_meta, feature_meta, representation_name,
#               task_name, outdir, ...)
#
# Dispatches to task-specific functions and saves standardized outputs.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

source_dir <- dirname(sys.frame(1)$ofile %||% "scripts/benchmark/benchmark_runner.R")
for (f in c("evaluation_helpers.R", "subset_strategies.R",
            "task_breast_subtype.R", "task_breast_vs_lung.R",
            "diagnostics.R")) {
  fp <- file.path(source_dir, f)
  if (file.exists(fp)) source(fp)
}

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
run_benchmark <- function(matrix,
                          sample_meta,
                          feature_meta,
                          representation_name,
                          task_name = c("breast_subtype", "breast_vs_lung"),
                          outdir = NULL,
                          marker_genes = NULL,
                          ...) {
  task_name <- match.arg(task_name)

  if (is.null(outdir)) {
    outdir <- file.path("reports", "benchmark_master", "representation_level_da",
                        representation_name, task_name)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.character(matrix) && file.exists(matrix)) {
    mat_dt <- fread(matrix, header = TRUE)
    gene_col <- names(mat_dt)[1]
    id_col <- if ("UniProtID" %in% names(mat_dt)) "UniProtID" else NULL
    sample_cols <- setdiff(names(mat_dt), c(gene_col, id_col))
    mat <- as.matrix(mat_dt[, ..sample_cols])
    rownames(mat) <- mat_dt[[gene_col]]
    matrix <- mat
  }

  if (is.character(sample_meta) && file.exists(sample_meta))
    sample_meta <- fread(sample_meta)
  if (is.character(feature_meta) && file.exists(feature_meta))
    feature_meta <- fread(feature_meta)

  result <- switch(task_name,
    breast_subtype = run_task_breast_subtype(
      matrix = matrix, sample_meta = sample_meta, feature_meta = feature_meta,
      representation_name = representation_name, outdir = outdir,
      marker_genes = marker_genes, ...
    ),
    breast_vs_lung = run_task_breast_vs_lung(
      matrix = matrix, sample_meta = sample_meta, feature_meta = feature_meta,
      representation_name = representation_name, outdir = outdir,
      marker_genes = marker_genes, ...
    )
  )

  # Save run metadata
  meta_lines <- c(
    paste("Representation:", representation_name),
    paste("Task:", task_name),
    paste("Date:", Sys.time()),
    paste("Output directory:", outdir),
    paste("Matrix dimensions:", nrow(matrix), "genes x", ncol(matrix), "samples"),
    paste("Inference type: representation-level limma"),
    paste("NOT native-domain TMT-aware inference")
  )
  writeLines(meta_lines, file.path(outdir, "run_metadata.txt"))

  result
}
