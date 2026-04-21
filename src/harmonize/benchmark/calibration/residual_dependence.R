#!/usr/bin/env Rscript
# =============================================================================
# Residual Dependence Diagnostic
# =============================================================================
# Quantifies residual correlation structure in limma model residuals. High
# residual correlation implies non-independence of samples (e.g. from shared
# TMT plexes), reducing effective sample size.
#
# Usage (CLI):
#   Rscript residual_dependence.R \
#     --matrix input.csv --meta meta.csv \
#     --contrast-a Basal --contrast-b Luminal \
#     --outdir /path
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

compute_residual_dependence <- function(matrix, sample_meta,
                                         contrast_a, contrast_b,
                                         domain = NULL) {
  # ---------------------------------------------------------------------------
  # Fit ~0 + condition via limma, extract residuals, compute pairwise sample
  # correlation, and report mean |correlation|, effective sample size.
  #
  # If domain is non-NULL, subset to that domain first.
  #
  # Returns list($residual_corr_matrix, $summary)
  # ---------------------------------------------------------------------------
  meta <- as.data.table(sample_meta)

  if (!is.null(domain)) {
    dom_u <- toupper(domain)
    if ("domain_col" %in% names(meta)) {
      meta <- meta[toupper(domain_col) == dom_u]
    } else if ("domain" %in% names(meta)) {
      meta <- meta[toupper(domain) == dom_u]
    }
  }

  meta <- meta[condition %in% c(contrast_a, contrast_b)]
  sids <- intersect(meta$sample_id, colnames(matrix))
  if (length(sids) < 6)
    stop("Too few samples (", length(sids), ") for residual dependence analysis")

  mat <- matrix[, sids, drop = FALSE]
  groups <- meta[match(sids, sample_id), condition]

  # Filter genes
  keep <- rowSums(!is.na(mat)) >= ncol(mat) * 0.5
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 50)
    stop("Too few genes (", nrow(mat), ") after filtering")

  # Fit model
  gf <- factor(groups, levels = c(contrast_a, contrast_b))
  design <- model.matrix(~ 0 + gf)
  colnames(design) <- levels(gf)

  fit <- lmFit(mat, design)
  residuals <- residuals(fit, mat)

  # Pairwise sample correlation on residuals
  resid_cor <- cor(residuals, use = "pairwise.complete.obs")

  # Extract off-diagonal elements
  n <- ncol(resid_cor)
  off_diag <- resid_cor[upper.tri(resid_cor)]

  mean_abs_corr <- mean(abs(off_diag), na.rm = TRUE)
  mean_corr     <- mean(off_diag, na.rm = TRUE)
  median_corr   <- median(off_diag, na.rm = TRUE)

  # Effective sample size (Kish-like approximation)
  # n_eff = n / (1 + (n-1) * mean_abs_rho)
  mean_abs_rho <- mean_abs_corr
  n_eff <- n / (1 + (n - 1) * mean_abs_rho)

  summary_dt <- data.table(
    domain = if (!is.null(domain)) domain else "all",
    n_samples = n,
    n_genes = nrow(mat),
    mean_abs_residual_corr = mean_abs_corr,
    mean_residual_corr = mean_corr,
    median_residual_corr = median_corr,
    max_abs_residual_corr = max(abs(off_diag), na.rm = TRUE),
    effective_n = round(n_eff, 1),
    effective_n_ratio = round(n_eff / n, 3)
  )

  list(
    residual_corr_matrix = resid_cor,
    summary = summary_dt
  )
}


# ── CLI entry point ──────────────────────────────────────────────────────────
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  parse_arg <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) == 0) return(NULL)
    args[idx + 1]
  }

  matrix_path <- parse_arg("--matrix")
  meta_path   <- parse_arg("--meta")
  contrast_a  <- parse_arg("--contrast-a")
  contrast_b  <- parse_arg("--contrast-b")
  domain_arg  <- parse_arg("--domain")
  outdir      <- parse_arg("--outdir")

  if (is.null(matrix_path) || is.null(meta_path) || is.null(outdir))
    stop("Required: --matrix, --meta, --outdir")
  if (is.null(contrast_a) || is.null(contrast_b))
    stop("Required: --contrast-a and --contrast-b")

  mat  <- as.matrix(fread(matrix_path), rownames = 1)
  meta <- fread(meta_path)

  dom_label <- if (!is.null(domain_arg)) domain_arg else "all"
  cat("Computing residual dependence for domain:", dom_label, "\n")

  res <- compute_residual_dependence(mat, meta, contrast_a, contrast_b,
                                      domain = domain_arg)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(res$residual_corr_matrix, keep.rownames = TRUE),
         file.path(outdir, paste0("residual_corr_matrix_", tolower(dom_label), ".csv")))
  fwrite(res$summary,
         file.path(outdir, paste0("residual_dependence_", tolower(dom_label), ".csv")))
  cat("Saved to:", outdir, "\n")
  print(res$summary)
}
