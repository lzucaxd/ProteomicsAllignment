#!/usr/bin/env Rscript
# =============================================================================
# Common method interface for CPTAC–CCLE benchmark representations
# =============================================================================
#
# All methods return a named list with identical structure:
#
#   $matrix        — numeric matrix (genes × samples), the transformed representation
#   $sample_meta   — data.frame with at least: sample_id, domain (CPTAC|CCLE), condition
#   $feature_meta  — data.frame with at least: gene, included (logical), exclusion_reason
#   $method_name   — character, e.g. "raw", "bridge_shift", "bridge_shift_scale", "celligner"
#   $method_notes  — character vector of free-text notes
#   $qc_paths      — named character vector of paths to QC outputs
#
# Benchmark runner calls:
#   result <- run_<method>_representation(cptac_mat, ccle_mat, sample_meta, ...)
#
# and then passes result$matrix, result$sample_meta, result$feature_meta
# to the evaluation layer unchanged.
#
# All matrices are on a log2 abundance scale unless documented otherwise.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Helper: load a gene matrix CSV → numeric matrix (genes in rows, samples in cols)
# ---------------------------------------------------------------------------
load_gene_matrix <- function(path) {
  dt <- fread(path, header = TRUE)
  gene_col <- names(dt)[1]
  id_col <- if ("UniProtID" %in% names(dt)) "UniProtID" else NULL
  sample_cols <- setdiff(names(dt), c(gene_col, id_col))
  mat <- as.matrix(dt[, ..sample_cols])
  rownames(mat) <- dt[[gene_col]]
  mat
}

# ---------------------------------------------------------------------------
# Helper: intersect gene features across two matrices
# ---------------------------------------------------------------------------
intersect_features <- function(mat_a, mat_b, min_obs_frac = 0.1) {
  shared <- intersect(rownames(mat_a), rownames(mat_b))
  if (length(shared) == 0) stop("No shared gene features between the two matrices.")
  mat_a <- mat_a[shared, , drop = FALSE]
  mat_b <- mat_b[shared, , drop = FALSE]
  obs_a <- rowMeans(!is.na(mat_a))
  obs_b <- rowMeans(!is.na(mat_b))
  keep <- obs_a >= min_obs_frac & obs_b >= min_obs_frac
  list(
    mat_a = mat_a[keep, , drop = FALSE],
    mat_b = mat_b[keep, , drop = FALSE],
    genes_shared = rownames(mat_a)[keep],
    genes_dropped = shared[!keep],
    drop_reason = ifelse(obs_a[!keep] < min_obs_frac, "low_obs_cptac",
                         ifelse(obs_b[!keep] < min_obs_frac, "low_obs_ccle", "low_obs_both"))
  )
}

# ---------------------------------------------------------------------------
# Helper: build combined matrix + sample metadata from CPTAC + CCLE
# ---------------------------------------------------------------------------
combine_domains <- function(cptac_mat, ccle_mat, cptac_meta, ccle_meta) {
  genes <- intersect(rownames(cptac_mat), rownames(ccle_mat))
  combined <- cbind(cptac_mat[genes, , drop = FALSE], ccle_mat[genes, , drop = FALSE])
  meta <- rbind(cptac_meta, ccle_meta)
  list(matrix = combined, sample_meta = meta, genes = genes)
}

# ---------------------------------------------------------------------------
# Helper: build the standard return object
# ---------------------------------------------------------------------------
make_method_result <- function(matrix, sample_meta, feature_meta, method_name,
                               method_notes = character(), qc_paths = character()) {
  list(
    matrix       = matrix,
    sample_meta  = sample_meta,
    feature_meta = feature_meta,
    method_name  = method_name,
    method_notes = method_notes,
    qc_paths     = qc_paths
  )
}

# ---------------------------------------------------------------------------
# Helper: save method result to disk (benchmark contract)
# ---------------------------------------------------------------------------
save_method_result <- function(result, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  mat_dt <- data.table(Gene = rownames(result$matrix), result$matrix)
  fwrite(mat_dt, file.path(outdir, "transformed_matrix.csv"))
  fwrite(as.data.table(result$sample_meta), file.path(outdir, "sample_metadata.csv"))
  fwrite(as.data.table(result$feature_meta), file.path(outdir, "feature_metadata.csv"))
  writeLines(result$method_notes, file.path(outdir, "method_notes.txt"))
  message("Method '", result$method_name, "' outputs saved to ", outdir)
  invisible(outdir)
}
