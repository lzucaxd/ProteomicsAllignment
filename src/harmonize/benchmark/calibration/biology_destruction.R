#!/usr/bin/env Rscript
# =============================================================================
# Biology Destruction Check
# =============================================================================
# Compares native-domain (or raw) DA results against post-harmonization DA
# results to quantify gene loss and fold-change shrinkage. A method that
# destroys biology would lose significant genes or dramatically shrink their
# effect sizes.
#
# Usage (CLI):
#   Rscript biology_destruction.R \
#     --native-da native_da.csv --method-da method_da.csv \
#     --method bridge_shift --outdir /path
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

check_biology_destruction <- function(native_da, method_da, method_name = "method",
                                       p_thresholds = c(0.01, 0.05, 0.1),
                                       fc_thresholds = c(0.25, 0.5, 1.0)) {
  # ---------------------------------------------------------------------------
  # Compare native-domain DA (reference) against post-harmonization DA.
  #
  # Args:
  #   native_da  : data.table with columns gene, logFC, adj.P.Val
  #   method_da  : data.table with columns gene, logFC, adj.P.Val
  #   method_name: label for the harmonization method
  #
  # Returns list($per_gene, $summary, $sensitivity_grid)
  # ---------------------------------------------------------------------------
  native_da <- as.data.table(native_da)
  method_da <- as.data.table(method_da)

  # Merge on shared genes
  shared <- merge(
    native_da[, .(gene, native_logFC = logFC, native_adjP = adj.P.Val)],
    method_da[, .(gene, method_logFC = logFC, method_adjP = adj.P.Val)],
    by = "gene"
  )

  if (nrow(shared) == 0) {
    return(list(
      per_gene = data.table(),
      summary = data.table(method = method_name, error = "no shared genes"),
      sensitivity_grid = data.table()
    ))
  }

  # Per-gene metrics
  shared[, same_direction := sign(native_logFC) == sign(method_logFC)]
  shared[, fc_ratio := ifelse(native_logFC == 0, NA_real_,
                               abs(method_logFC) / abs(native_logFC))]
  shared[, fc_shrinkage := 1 - pmin(fc_ratio, 1, na.rm = TRUE)]

  per_gene <- copy(shared)
  per_gene$method <- method_name

  # Sensitivity grid across threshold combinations
  grid_rows <- list()
  for (p_thr in p_thresholds) {
    for (fc_thr in fc_thresholds) {
      # Significant in native
      sig_native <- shared[native_adjP < p_thr & abs(native_logFC) > fc_thr]
      n_sig_native <- nrow(sig_native)

      if (n_sig_native == 0) {
        grid_rows[[length(grid_rows) + 1]] <- data.table(
          method = method_name, p_threshold = p_thr, fc_threshold = fc_thr,
          n_sig_native = 0, n_retained = 0, retention_rate = NA_real_,
          mean_fc_shrinkage = NA_real_, median_fc_ratio = NA_real_,
          direction_agreement = NA_real_
        )
        next
      }

      # Still significant after harmonization
      sig_retained <- sig_native[method_adjP < p_thr & abs(method_logFC) > fc_thr]
      n_retained <- nrow(sig_retained)

      grid_rows[[length(grid_rows) + 1]] <- data.table(
        method = method_name,
        p_threshold = p_thr,
        fc_threshold = fc_thr,
        n_sig_native = n_sig_native,
        n_retained = n_retained,
        retention_rate = n_retained / n_sig_native,
        mean_fc_shrinkage = mean(sig_native$fc_shrinkage, na.rm = TRUE),
        median_fc_ratio = median(sig_native$fc_ratio, na.rm = TRUE),
        direction_agreement = mean(sig_native$same_direction, na.rm = TRUE)
      )
    }
  }

  sensitivity_grid <- rbindlist(grid_rows)

  # Default summary at p < 0.05, |FC| > 0.5
  default_row <- sensitivity_grid[p_threshold == 0.05 & fc_threshold == 0.5]
  if (nrow(default_row) == 0) default_row <- sensitivity_grid[1]

  summary_dt <- data.table(
    method = method_name,
    n_shared_genes = nrow(shared),
    overall_fc_correlation = cor(shared$native_logFC, shared$method_logFC,
                                 use = "pairwise.complete.obs"),
    overall_direction_agreement = mean(shared$same_direction, na.rm = TRUE),
    default_n_sig_native = default_row$n_sig_native,
    default_retention_rate = default_row$retention_rate,
    default_mean_fc_shrinkage = default_row$mean_fc_shrinkage,
    default_median_fc_ratio = default_row$median_fc_ratio
  )

  list(per_gene = per_gene, summary = summary_dt, sensitivity_grid = sensitivity_grid)
}


# ── CLI entry point ──────────────────────────────────────────────────────────
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  parse_arg <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) == 0) return(NULL)
    args[idx + 1]
  }

  native_path <- parse_arg("--native-da")
  method_path <- parse_arg("--method-da")
  method_name <- parse_arg("--method") %||% "method"
  outdir      <- parse_arg("--outdir")

  if (is.null(native_path) || is.null(method_path) || is.null(outdir))
    stop("Required: --native-da, --method-da, --outdir")

  native_da <- fread(native_path)
  method_da <- fread(method_path)

  cat("Checking biology destruction for:", method_name, "\n")
  res <- check_biology_destruction(native_da, method_da, method_name)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(res$per_gene, file.path(outdir, paste0("biology_destruction_", method_name, ".csv")))
  fwrite(res$summary, file.path(outdir, paste0("destruction_summary_", method_name, ".csv")))
  fwrite(res$sensitivity_grid, file.path(outdir, paste0("destruction_grid_", method_name, ".csv")))
  cat("Saved to:", outdir, "\n")
  print(res$summary)
}
