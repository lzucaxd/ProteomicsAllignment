#!/usr/bin/env Rscript
# =============================================================================
# Marker Sanity Rate
# =============================================================================
# Generalized marker direction agreement check for any method x domain x task.
# Extends validate_contrast_direction() to produce a standard metric row
# suitable for inclusion in the benchmark comparison table.
#
# Usage (CLI):
#   Rscript marker_sanity.R \
#     --da-result da_limma_result.csv \
#     --markers "ESR1,PGR,GATA3,FOXA1,EGFR,KRT5,KRT17,FOXC1" \
#     --expected-signs "1,1,1,1,-1,-1,-1,-1" \
#     --method raw --domain CPTAC --task breast_subtype \
#     --outdir /path
#
# expected_sign: +1 means the gene should have positive logFC under the same
# contrast as limma_da_wrapper.R (contrast_b - contrast_a). Breast subtype:
# Luminal - Basal => +1 = higher in Luminal.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

compute_marker_sanity <- function(da_result, expected_directions,
                                   method = "unknown", domain = "unknown",
                                   task = "unknown") {
  # ---------------------------------------------------------------------------
  # Check fraction of marker genes with correct logFC sign.
  #
  # Args:
  #   da_result           : data.table with at least gene, logFC columns
  #   expected_directions : data.table with gene, expected_sign (+1/-1)
  #   method, domain, task: labels for the output row
  #

  # Returns list($marker_table, $summary)
  # ---------------------------------------------------------------------------
  da <- as.data.table(da_result)
  exp <- as.data.table(expected_directions)

  merged <- merge(da[, .(gene, logFC)], exp, by = "gene")

  if (nrow(merged) == 0) {
    summary_dt <- data.table(
      method = method, domain = domain, task = task,
      n_markers_tested = 0L, n_correct = 0L,
      marker_sanity_rate = NA_real_
    )
    return(list(marker_table = merged, summary = summary_dt))
  }

  merged[, observed_sign := sign(logFC)]
  merged[, correct := (observed_sign == expected_sign)]

  n_tested  <- nrow(merged)
  n_correct <- sum(merged$correct)
  sanity_rate <- n_correct / n_tested

  marker_table <- copy(merged)
  marker_table$method <- method
  marker_table$domain <- domain
  marker_table$task   <- task

  summary_dt <- data.table(
    method = method, domain = domain, task = task,
    n_markers_tested = n_tested,
    n_correct = n_correct,
    marker_sanity_rate = sanity_rate
  )

  list(marker_table = marker_table, summary = summary_dt)
}


# ── CLI entry point ──────────────────────────────────────────────────────────
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  parse_arg <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) == 0) return(NULL)
    args[idx + 1]
  }

  da_path   <- parse_arg("--da-result")
  mk_str    <- parse_arg("--markers")
  es_str    <- parse_arg("--expected-signs")
  method    <- parse_arg("--method") %||% "unknown"
  domain    <- parse_arg("--domain") %||% "unknown"
  task      <- parse_arg("--task") %||% "unknown"
  outdir    <- parse_arg("--outdir")

  if (is.null(da_path) || is.null(mk_str) || is.null(es_str) || is.null(outdir))
    stop("Required: --da-result, --markers, --expected-signs, --outdir")

  da <- fread(da_path)
  mk_genes  <- trimws(strsplit(mk_str, ",")[[1]])
  exp_signs <- as.integer(trimws(strsplit(es_str, ",")[[1]]))
  exp_dirs  <- data.table(gene = mk_genes, expected_sign = exp_signs)

  res <- compute_marker_sanity(da, exp_dirs, method, domain, task)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fname <- paste0("marker_sanity_", method, "_", tolower(domain), "_", task, ".csv")
  fwrite(res$marker_table, file.path(outdir, fname))
  fwrite(res$summary, file.path(outdir, paste0("marker_sanity_summary_", method, "_",
                                                 tolower(domain), "_", task, ".csv")))
  cat("Marker sanity:", res$summary$marker_sanity_rate, "\n")
  cat("Saved to:", outdir, "\n")
}
