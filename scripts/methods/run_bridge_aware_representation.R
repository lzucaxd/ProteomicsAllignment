#!/usr/bin/env Rscript
# =============================================================================
# Method 1 — Bridge-Aware Per-Protein Representation
# =============================================================================
# Aligns CPTAC and CCLE using bridge/reference channel information per protein.
#
# Two modes:
#   "shift_only"      — subtract per-protein domain-specific bridge median,
#                        re-center both domains to the global bridge median
#   "shift_and_scale" — additionally rescale by per-protein bridge MAD ratio
#
# Usage:
#   source("scripts/methods/method_interface.R")
#   source("scripts/methods/run_bridge_aware_representation.R")
#   result <- run_bridge_aware_representation(
#     cptac_mat, ccle_mat, cptac_meta, ccle_meta,
#     bridge_cptac, bridge_ccle, outdir, mode = "shift_only"
#   )
#
# Bridge inputs:
#   bridge_cptac — numeric vector or matrix (genes × bridge_samples) of CPTAC
#                  bridge/pool channel abundances (log2 scale)
#   bridge_ccle  — same for CCLE
#
# If bridge inputs are data.frames with a Gene/gene column, the first column is
# used as gene names and the rest as bridge observations.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Extract per-protein bridge statistics
# ---------------------------------------------------------------------------
compute_bridge_stats <- function(bridge_input, label) {
  if (is.data.frame(bridge_input) || is.data.table(bridge_input)) {
    bridge_input <- as.data.table(bridge_input)
    gene_col <- names(bridge_input)[1]
    genes <- as.character(bridge_input[[gene_col]])
    mat <- as.matrix(bridge_input[, -1, with = FALSE])
    rownames(mat) <- genes
  } else if (is.matrix(bridge_input)) {
    mat <- bridge_input
  } else if (is.numeric(bridge_input) && !is.null(names(bridge_input))) {
    mat <- matrix(bridge_input, ncol = 1, dimnames = list(names(bridge_input), label))
  } else {
    stop("bridge_", label, " must be a named numeric vector, matrix, or data.frame")
  }
  data.frame(
    gene = rownames(mat),
    bridge_median = apply(mat, 1, median, na.rm = TRUE),
    bridge_mad = apply(mat, 1, mad, na.rm = TRUE),
    bridge_n_obs = rowSums(!is.na(mat)),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
run_bridge_aware_representation <- function(cptac_mat, ccle_mat,
                                             cptac_meta, ccle_meta,
                                             bridge_cptac, bridge_ccle,
                                             outdir = "reports/benchmark_master/methods/bridge_aware",
                                             mode = c("shift_only", "shift_and_scale"),
                                             min_obs_frac = 0.1,
                                             min_bridge_obs = 2,
                                             min_bridge_mad = 0.01) {
  mode <- match.arg(mode)

  # Intersect features
  shared <- intersect_features(cptac_mat, ccle_mat, min_obs_frac = min_obs_frac)
  cptac_sh <- shared$mat_a
  ccle_sh <- shared$mat_b
  genes <- rownames(cptac_sh)

  # Bridge statistics
  bs_cptac <- compute_bridge_stats(bridge_cptac, "cptac")
  bs_ccle <- compute_bridge_stats(bridge_ccle, "ccle")

  # Align bridge stats to shared genes
  rownames(bs_cptac) <- bs_cptac$gene
  rownames(bs_ccle) <- bs_ccle$gene
  bridge_genes <- intersect(genes, intersect(bs_cptac$gene, bs_ccle$gene))

  # Determine which genes have usable bridge information
  usable <- bridge_genes[
    bs_cptac[bridge_genes, "bridge_n_obs"] >= min_bridge_obs &
    bs_ccle[bridge_genes, "bridge_n_obs"] >= min_bridge_obs
  ]

  if (mode == "shift_and_scale") {
    usable <- usable[
      bs_cptac[usable, "bridge_mad"] >= min_bridge_mad &
      bs_ccle[usable, "bridge_mad"] >= min_bridge_mad
    ]
  }

  no_bridge <- setdiff(genes, usable)

  # Compute offsets and scaling factors
  offset_cptac <- bs_cptac[usable, "bridge_median"]
  offset_ccle <- bs_ccle[usable, "bridge_median"]
  global_bridge_mean <- (offset_cptac + offset_ccle) / 2

  # --- shift_only: X_aligned = X - domain_bridge_median + global_bridge_mean ---
  cptac_aligned <- cptac_sh
  ccle_aligned <- ccle_sh

  cptac_aligned[usable, ] <- sweep(cptac_sh[usable, , drop = FALSE], 1, offset_cptac) +
                              global_bridge_mean
  ccle_aligned[usable, ] <- sweep(ccle_sh[usable, , drop = FALSE], 1, offset_ccle) +
                             global_bridge_mean

  scale_factors <- rep(NA_real_, length(usable))
  names(scale_factors) <- usable

  if (mode == "shift_and_scale") {
    mad_cptac <- bs_cptac[usable, "bridge_mad"]
    mad_ccle <- bs_ccle[usable, "bridge_mad"]
    target_mad <- (mad_cptac + mad_ccle) / 2

    sf_cptac <- target_mad / mad_cptac
    sf_ccle <- target_mad / mad_ccle

    # Apply shift first (already done), then scale around bridge center
    cptac_aligned[usable, ] <- sweep(cptac_aligned[usable, , drop = FALSE],
                                      1, global_bridge_mean) * sf_cptac + global_bridge_mean
    ccle_aligned[usable, ] <- sweep(ccle_aligned[usable, , drop = FALSE],
                                     1, global_bridge_mean) * sf_ccle + global_bridge_mean
    scale_factors <- sf_cptac / sf_ccle
  }

  # Genes without bridge info: pass through unmodified
  # (documented in feature_meta as included but not bridge-corrected)

  combined <- combine_domains(cptac_aligned, ccle_aligned, cptac_meta, ccle_meta)

  # Feature metadata
  feature_meta <- data.frame(
    gene = c(genes, shared$genes_dropped),
    included = c(rep(TRUE, length(genes)), rep(FALSE, length(shared$genes_dropped))),
    bridge_corrected = c(genes %in% usable, rep(FALSE, length(shared$genes_dropped))),
    exclusion_reason = c(
      ifelse(genes %in% usable, NA_character_,
             ifelse(genes %in% bridge_genes, "bridge_obs_or_mad_too_low", "no_bridge_data")),
      shared$drop_reason
    ),
    stringsAsFactors = FALSE
  )

  # QC: offset and scale distributions
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  qc_dt <- data.table(
    gene = usable,
    offset_cptac = offset_cptac,
    offset_ccle = offset_ccle,
    shift_applied = offset_cptac - offset_ccle,
    scale_factor = scale_factors
  )
  qc_path <- file.path(outdir, "bridge_correction_details.csv")
  fwrite(qc_dt, qc_path)

  # QC summary plots
  qc_plot_path <- file.path(outdir, "bridge_offset_distribution.png")
  tryCatch({
    png(qc_plot_path, width = 1000, height = 500, res = 100)
    par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
    hist(qc_dt$shift_applied, breaks = 60, col = "steelblue",
         main = "Per-protein bridge offset\n(CPTAC median - CCLE median)",
         xlab = "Bridge offset (log2)")
    abline(v = median(qc_dt$shift_applied, na.rm = TRUE), col = "red", lwd = 2)
    if (mode == "shift_and_scale") {
      hist(log2(qc_dt$scale_factor), breaks = 60, col = "darkorange",
           main = "Per-protein scale factor ratio\n(log2, CPTAC/CCLE bridge MAD)",
           xlab = "log2(scale ratio)")
      abline(v = 0, col = "red", lwd = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "Scale correction\nnot applied\n(shift_only mode)", cex = 1.2)
    }
    dev.off()
  }, error = function(e) warning("QC plot failed: ", conditionMessage(e)))

  notes <- c(
    paste0("Method: bridge_aware (", mode, ")"),
    paste0("Date: ", Sys.time()),
    "",
    "Description:",
    if (mode == "shift_only") c(
      "  Per-protein shift correction using bridge/reference channel medians.",
      "  For each gene g with usable bridge data:",
      "    X_aligned[g, s] = X[g, s] - bridge_median_domain[g] + global_bridge_mean[g]",
      "  where global_bridge_mean = (bridge_median_cptac + bridge_median_ccle) / 2."
    ) else c(
      "  Per-protein shift + scale correction using bridge channel medians and MADs.",
      "  Shift: same as shift_only.",
      "  Scale: after shifting, rescale around global_bridge_mean using MAD ratio:",
      "    X_final[g, s] = (X_shifted[g, s] - center) * (target_MAD / domain_MAD) + center"
    ),
    "",
    paste0("Genes with bridge correction: ", length(usable)),
    paste0("Genes without bridge (passed through): ", length(no_bridge)),
    paste0("Genes excluded (low obs): ", length(shared$genes_dropped)),
    paste0("Median offset (CPTAC - CCLE): ",
           round(median(qc_dt$shift_applied, na.rm = TRUE), 4)),
    if (mode == "shift_and_scale")
      paste0("Median log2(scale ratio): ",
             round(median(log2(qc_dt$scale_factor), na.rm = TRUE), 4))
    else NULL,
    "",
    "Assumptions:",
    "  - Bridge/reference channels measure the same pooled standard in both domains",
    "  - Per-protein bridge medians are stable estimators of systematic domain offset",
    "  - The shift (and optionally scale) is approximately constant across samples within a domain",
    "",
    "Scale: log2 abundance (same as input; bridge-centered)",
    "",
    "Genes without bridge data:",
    "  Passed through WITHOUT correction. Feature metadata marks these genes",
    "  (bridge_corrected = FALSE). Downstream analyses may want to exclude them."
  )

  result <- make_method_result(
    matrix       = combined$matrix,
    sample_meta  = combined$sample_meta,
    feature_meta = feature_meta,
    method_name  = paste0("bridge_", gsub("_and_", "_", mode)),
    method_notes = notes,
    qc_paths     = c(correction_details = qc_path, offset_plot = qc_plot_path)
  )
  save_method_result(result, outdir)
  result
}
