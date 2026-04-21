#!/usr/bin/env Rscript
# =============================================================================
# Bridge-Aware Per-Protein Correction Engine
# =============================================================================
# Implements shift-only and shift+scale bridge-aware harmonization.
# Reads bridge summaries, applies corrections, produces QC outputs.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

.local_args <- commandArgs(trailingOnly = FALSE)
.local_file <- .local_args[startsWith(.local_args, "--file=")]
if (length(.local_file)) {
  .local_bench <- dirname(normalizePath(sub("^--file=", "", .local_file[1L])))
} else {
  .local_bench <- normalizePath(file.path(getwd(), "scripts", "benchmark"), mustWork = FALSE)
}
source(file.path(.local_bench, "harmonize_paths.R"))
REPO <- harmonize_repo_root()
BRIDGE_DIR <- file.path(REPO, "reports/benchmark_master/methods/bridge_aware")

# ─── Core correction function ────────────────────────────────────────────
run_bridge_aware_representation <- function(
  cptac_mat,            # genes × samples matrix
  ccle_mat,             # genes × samples matrix
  cptac_bridge,         # data.table with GeneSymbol, bridge_median, bridge_mad, bridge_n
  ccle_bridge,          # same
  config = list()
) {
  # Parse config with defaults
  mode             <- config$mode             %||% "shift_only"
  center_est       <- config$center_estimator %||% "median"
  scale_est        <- config$scale_estimator  %||% "mad"
  min_bridge_obs   <- config$min_bridge_obs   %||% 3L
  min_spread       <- config$min_spread       %||% 0.01
  max_offset       <- config$max_offset       %||% 10.0
  reference_domain <- config$reference_domain %||% "CPTAC"

  cat("  Bridge-aware correction config:\n")
  cat("    mode:", mode, "\n")
  cat("    center:", center_est, "\n")
  cat("    scale:", scale_est, "\n")
  cat("    min_bridge_obs:", min_bridge_obs, "\n")

  # Select center and spread columns
  center_col <- paste0("bridge_", center_est)
  spread_col <- paste0("bridge_", scale_est)

  # Shared genes across gene matrices and bridge summaries
  shared_genes <- Reduce(intersect, list(
    rownames(cptac_mat),
    rownames(ccle_mat),
    cptac_bridge[!is.na(GeneSymbol), GeneSymbol],
    ccle_bridge[!is.na(GeneSymbol), GeneSymbol]
  ))
  cat("    Shared genes (matrix ∩ bridge):", length(shared_genes), "\n")

  # Build per-gene correction table
  cb <- cptac_bridge[GeneSymbol %in% shared_genes]
  cb <- cb[!duplicated(GeneSymbol)]
  setkey(cb, GeneSymbol)

  eb <- ccle_bridge[GeneSymbol %in% shared_genes]
  eb <- eb[!duplicated(GeneSymbol)]
  setkey(eb, GeneSymbol)

  offsets <- data.table(gene = shared_genes)
  offsets[, cptac_center := cb[gene, get(center_col)]]
  offsets[, ccle_center := eb[gene, get(center_col)]]
  offsets[, cptac_n := cb[gene, bridge_n]]
  offsets[, ccle_n := eb[gene, bridge_n]]

  # Determine which proteins can be corrected
  offsets[, sufficient := (cptac_n >= min_bridge_obs) & (ccle_n >= min_bridge_obs)]
  offsets[, offset := cptac_center - ccle_center]
  offsets[, abs_offset := abs(offset)]
  # Extreme = offset deviates more than max_offset from the median offset
  # (not absolute value — the bulk shift is expected)
  med_off <- median(offsets[sufficient == TRUE, offset], na.rm = TRUE)
  offsets[, deviation_from_median := abs(offset - med_off)]
  offsets[, extreme := deviation_from_median > max_offset]

  # Flags
  offsets[, status := "corrected"]
  offsets[sufficient == FALSE, status := "insufficient_bridge"]
  offsets[is.na(offset), status := "missing_bridge"]

  # Genes in matrices but not in bridge
  all_shared_mat <- intersect(rownames(cptac_mat), rownames(ccle_mat))
  no_bridge <- setdiff(all_shared_mat, shared_genes)

  extra_rows <- data.table(
    gene = no_bridge,
    cptac_center = NA_real_, ccle_center = NA_real_,
    cptac_n = 0L, ccle_n = 0L,
    sufficient = FALSE, offset = NA_real_, abs_offset = NA_real_,
    extreme = FALSE, status = "no_bridge_data"
  )
  offsets <- rbind(offsets, extra_rows, fill = TRUE)

  # Scale factors (for shift+scale mode)
  scale_dt <- NULL
  if (mode == "shift_and_scale") {
    offsets[, cptac_spread := cb[gene, get(spread_col)]]
    offsets[, ccle_spread := eb[gene, get(spread_col)]]
    offsets[, scale_factor := ifelse(ccle_spread > min_spread & cptac_spread > min_spread,
                                     cptac_spread / ccle_spread, NA_real_)]
    offsets[is.na(scale_factor) & sufficient == TRUE, status := "unstable_spread_fallback_shift"]
    scale_dt <- offsets[, .(gene, cptac_spread, ccle_spread, scale_factor,
                            scale_applied = !is.na(scale_factor) & sufficient == TRUE)]
  }

  # ── Apply correction ──────────────────────────────────────────────────
  correctable <- offsets[status == "corrected", gene]

  # Start with CCLE matrix, apply corrections
  ccle_adj <- ccle_mat[all_shared_mat, , drop = FALSE]
  cptac_sub <- cptac_mat[all_shared_mat, , drop = FALSE]

  for (g in correctable) {
    delta <- offsets[gene == g, offset]
    if (is.na(delta)) next

    if (mode == "shift_and_scale") {
      sf <- offsets[gene == g, scale_factor]
      if (!is.na(sf)) {
        ccle_center_g <- offsets[gene == g, ccle_center]
        cptac_center_g <- offsets[gene == g, cptac_center]
        ccle_adj[g, ] <- (ccle_adj[g, ] - ccle_center_g) / offsets[gene == g, ccle_spread] *
                          offsets[gene == g, cptac_spread] + cptac_center_g
      } else {
        ccle_adj[g, ] <- ccle_adj[g, ] + delta
      }
    } else {
      ccle_adj[g, ] <- ccle_adj[g, ] + delta
    }
  }

  # Combine into one matrix
  combined <- cbind(cptac_sub, ccle_adj)

  # QC
  n_corrected <- sum(offsets$status == "corrected")
  n_insuff <- sum(offsets$status == "insufficient_bridge")
  n_no_bridge <- sum(offsets$status == "no_bridge_data")
  n_missing <- sum(offsets$status == "missing_bridge")
  n_extreme <- sum(offsets[status == "corrected", extreme], na.rm = TRUE)
  n_unstable <- sum(offsets$status == "unstable_spread_fallback_shift")

  qc <- list(
    n_shared_genes = length(all_shared_mat),
    n_with_bridge = length(shared_genes),
    n_corrected = n_corrected,
    n_insufficient_bridge = n_insuff,
    n_no_bridge_data = n_no_bridge,
    n_missing_bridge = n_missing,
    n_extreme_offset = n_extreme,
    n_unstable_spread = n_unstable,
    median_offset = median(offsets[status == "corrected", offset], na.rm = TRUE),
    mad_offset = mad(offsets[status == "corrected", offset], na.rm = TRUE),
    q05_offset = quantile(offsets[status == "corrected", offset], 0.05, na.rm = TRUE),
    q95_offset = quantile(offsets[status == "corrected", offset], 0.95, na.rm = TRUE)
  )

  notes <- sprintf(
    "Bridge-aware %s | %d corrected | %d skipped (insuff: %d, no_bridge: %d, missing: %d) | %d extreme | median offset: %.3f",
    mode, n_corrected, n_insuff + n_no_bridge + n_missing, n_insuff, n_no_bridge, n_missing,
    n_extreme, qc$median_offset
  )

  list(
    matrix = combined,
    offsets = offsets,
    scale_factors = scale_dt,
    qc = qc,
    notes = notes
  )
}

# Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# MAIN: Run bridge-aware correction
# =============================================================================
main <- function() {
  cat("\n", strrep("=", 70), "\n")
  cat("  BRIDGE-AWARE PER-PROTEIN CORRECTION\n")
  cat(strrep("=", 70), "\n\n")

  # ── Load bridge summaries ──────────────────────────────────────────────
  cat("Loading bridge summaries...\n")
  cptac_bridge <- fread(file.path(BRIDGE_DIR, "bridge_summary_cptac.tsv"))
  ccle_bridge <- fread(file.path(BRIDGE_DIR, "bridge_summary_ccle.tsv"))
  cat("  CPTAC bridge:", nrow(cptac_bridge), "proteins\n")
  cat("  CCLE bridge:", nrow(ccle_bridge), "proteins\n")
  cat("  CPTAC with GeneSymbol:", sum(!is.na(cptac_bridge$GeneSymbol)), "\n")
  cat("  CCLE with GeneSymbol:", sum(!is.na(ccle_bridge$GeneSymbol)), "\n")

  # ── Load gene matrices ────────────────────────────────────────────────
  cat("\nLoading gene matrices...\n")
  load_gm <- function(path) {
    dt <- fread(path, header = TRUE)
    id_cols <- intersect(c("GeneSymbol", "UniProtID", "Gene"), names(dt))
    scols <- setdiff(names(dt), id_cols)
    mat <- as.matrix(dt[, ..scols])
    rownames(mat) <- dt[[names(dt)[1]]]
    mat
  }

  # Load all available CPTAC studies and combine
  cptac_paths <- c(
    file.path(REPO, "data/results/PDC000120/gene_matrix.csv"),
    file.path(REPO, "data/results/PDC000153/gene_matrix.csv")
  )
  cptac_paths <- cptac_paths[file.exists(cptac_paths)]

  cptac_mats <- lapply(cptac_paths, load_gm)
  cat("  Loaded", length(cptac_mats), "CPTAC studies\n")
  for (i in seq_along(cptac_paths)) {
    cat("    ", basename(dirname(cptac_paths[i])), ":",
        nrow(cptac_mats[[i]]), "×", ncol(cptac_mats[[i]]), "\n")
  }

  # Combine CPTAC matrices:
  # - samples: concatenate columns across studies (union of sample IDs)
  # - genes: INTERSECTION across studies (no NA-filled union rows)
  #
  # Rationale: an outer gene-union introduces many genes that exist in only one
  # study (all-NaN in the other), which is usually *not* what we want for a
  # harmonized CPTAC reference matrix feeding bridge-aware correction.
  if (length(cptac_mats) == 1) {
    cptac_mat <- cptac_mats[[1]]
  } else {
    shared_genes_cptac <- Reduce(intersect, lapply(cptac_mats, rownames))
    shared_genes_cptac <- sort(shared_genes_cptac)
    all_samples <- unlist(lapply(cptac_mats, colnames))
    cptac_mat <- matrix(
      NA_real_,
      nrow = length(shared_genes_cptac),
      ncol = length(all_samples),
      dimnames = list(shared_genes_cptac, all_samples)
    )
    for (m in cptac_mats) {
      g <- intersect(rownames(m), shared_genes_cptac)
      cptac_mat[g, colnames(m)] <- m[g, , drop = FALSE]
    }
  }

  ccle_mat <- load_gm(file.path(REPO, "data/results/CCLE_corrected/gene_matrix.csv"))
  cat("  Combined CPTAC:", nrow(cptac_mat), "×", ncol(cptac_mat), "\n")
  cat("  CCLE:", nrow(ccle_mat), "×", ncol(ccle_mat), "\n")

  # ── Run shift-only ────────────────────────────────────────────────────
  cat("\n", strrep("-", 50), "\n  MODE: shift_only\n", strrep("-", 50), "\n")
  res_shift <- run_bridge_aware_representation(
    cptac_mat, ccle_mat, cptac_bridge, ccle_bridge,
    config = list(mode = "shift_only")
  )
  cat("\n  ", res_shift$notes, "\n")

  # Save outputs
  cat("\nSaving shift-only outputs...\n")
  out_mat <- data.table(GeneSymbol = rownames(res_shift$matrix), as.data.table(res_shift$matrix))
  fwrite(out_mat, file.path(BRIDGE_DIR, "bridge_aware_shift_only_matrix.csv"))
  fwrite(res_shift$offsets, file.path(BRIDGE_DIR, "bridge_offsets.tsv"), sep = "\t")

  # QC report
  qc <- res_shift$qc
  qc_lines <- c(
    "# Bridge-Aware Shift-Only QC Summary",
    "",
    sprintf("Total shared genes (matrix intersection): %d", qc$n_shared_genes),
    sprintf("Genes with bridge data in both domains: %d", qc$n_with_bridge),
    "",
    "## Correction Status",
    sprintf("- Corrected (shift applied): %d", qc$n_corrected),
    sprintf("- Insufficient bridge (< min_bridge_obs): %d", qc$n_insufficient_bridge),
    sprintf("- No bridge data: %d", qc$n_no_bridge_data),
    sprintf("- Missing bridge values: %d", qc$n_missing_bridge),
    "",
    "## Offset Distribution (corrected proteins only)",
    sprintf("- Median offset: %.4f", qc$median_offset),
    sprintf("- MAD of offsets: %.4f", qc$mad_offset),
    sprintf("- 5th percentile: %.4f", qc$q05_offset),
    sprintf("- 95th percentile: %.4f", qc$q95_offset),
    sprintf("- Extreme offsets (|Δ - median(Δ)| > 10): %d", qc$n_extreme_offset),
    "",
    "## Top 10 Largest Offsets",
    ""
  )
  top10 <- res_shift$offsets[status == "corrected"][order(-abs_offset)][1:min(10, .N)]
  for (i in seq_len(nrow(top10))) {
    qc_lines <- c(qc_lines, sprintf("- %s: Δ = %.4f (CPTAC=%.2f, CCLE=%.2f, n_cptac=%d, n_ccle=%d)",
                                      top10$gene[i], top10$offset[i],
                                      top10$cptac_center[i], top10$ccle_center[i],
                                      top10$cptac_n[i], top10$ccle_n[i]))
  }
  qc_lines <- c(qc_lines, "", "## Warnings", "")
  if (qc$n_extreme_offset > 0) {
    qc_lines <- c(qc_lines, sprintf("- %d proteins have extreme offsets (|Δ| > 10). Review these for biological plausibility.", qc$n_extreme_offset))
  }
  if (qc$n_no_bridge_data > 0) {
    qc_lines <- c(qc_lines, sprintf("- %d proteins in the gene matrix have no bridge data and were left uncorrected.", qc$n_no_bridge_data))
  }
  writeLines(qc_lines, file.path(BRIDGE_DIR, "bridge_aware_shift_only_qc.md"))

  # ── Run shift + scale ─────────────────────────────────────────────────
  cat("\n", strrep("-", 50), "\n  MODE: shift_and_scale\n", strrep("-", 50), "\n")
  res_scale <- run_bridge_aware_representation(
    cptac_mat, ccle_mat, cptac_bridge, ccle_bridge,
    config = list(mode = "shift_and_scale")
  )
  cat("\n  ", res_scale$notes, "\n")

  cat("\nSaving shift+scale outputs...\n")
  out_mat2 <- data.table(GeneSymbol = rownames(res_scale$matrix), as.data.table(res_scale$matrix))
  fwrite(out_mat2, file.path(BRIDGE_DIR, "bridge_aware_shift_scale_matrix.csv"))
  fwrite(res_scale$scale_factors, file.path(BRIDGE_DIR, "bridge_scale_factors.tsv"), sep = "\t")

  # Scale QC
  qc2 <- res_scale$qc
  sf_applied <- res_scale$scale_factors[scale_applied == TRUE]
  qc2_lines <- c(
    "# Bridge-Aware Shift+Scale QC Summary",
    "",
    sprintf("Total shared genes: %d", qc2$n_shared_genes),
    sprintf("Genes with bridge data: %d", qc2$n_with_bridge),
    "",
    "## Correction Status",
    sprintf("- Fully scaled (shift + scale): %d", nrow(sf_applied)),
    sprintf("- Shift-only fallback (unstable spread): %d", qc2$n_unstable_spread),
    sprintf("- Insufficient bridge: %d", qc2$n_insufficient_bridge),
    sprintf("- No bridge data: %d", qc2$n_no_bridge_data),
    "",
    "## Scale Factor Distribution (scaled proteins only)",
    if (nrow(sf_applied) > 0) c(
      sprintf("- Median scale factor: %.4f", median(sf_applied$scale_factor, na.rm = TRUE)),
      sprintf("- MAD of scale factors: %.4f", mad(sf_applied$scale_factor, na.rm = TRUE)),
      sprintf("- 5th percentile: %.4f", quantile(sf_applied$scale_factor, 0.05, na.rm = TRUE)),
      sprintf("- 95th percentile: %.4f", quantile(sf_applied$scale_factor, 0.95, na.rm = TRUE))
    ) else "- No proteins were scaled.",
    "",
    "## Top 10 Most Extreme Scale Factors",
    ""
  )
  if (nrow(sf_applied) > 0) {
    sf_applied[, deviation := abs(scale_factor - 1)]
    top10sf <- sf_applied[order(-deviation)][1:min(10, .N)]
    for (i in seq_len(nrow(top10sf))) {
      qc2_lines <- c(qc2_lines, sprintf("- %s: SF = %.4f (CPTAC spread=%.4f, CCLE spread=%.4f)",
                                          top10sf$gene[i], top10sf$scale_factor[i],
                                          top10sf$cptac_spread[i], top10sf$ccle_spread[i]))
    }
  }
  qc2_lines <- c(qc2_lines, "",
    "## Warnings",
    "",
    "Shift+scale correction is more aggressive than shift-only. Proteins with",
    "very different bridge spreads may be over-scaled. Review extreme scale factors.",
    if (qc2$n_unstable_spread > 0)
      sprintf("- %d proteins fell back to shift-only due to unstable spread.", qc2$n_unstable_spread)
  )
  writeLines(qc2_lines, file.path(BRIDGE_DIR, "bridge_aware_shift_scale_qc.md"))

  # ── Summary ───────────────────────────────────────────────────────────
  cat("\n", strrep("=", 70), "\n")
  cat("  BRIDGE-AWARE CORRECTION COMPLETE\n")
  cat(strrep("=", 70), "\n")
  cat("  Shift-only matrix:", nrow(res_shift$matrix), "×", ncol(res_shift$matrix), "\n")
  cat("  Shift+scale matrix:", nrow(res_scale$matrix), "×", ncol(res_scale$matrix), "\n")
  cat("  Outputs:", BRIDGE_DIR, "\n\n")
}

# =============================================================================
# CLI: apply bridge correction to a combined union task matrix (genes × samples)
# Usage:
#   Rscript bridge_aware_correction.R --union-task-matrix \
#     --repo ROOT --matrix path/to/matrix.csv --meta path/to/meta.csv \
#     --mode shift_only|shift_and_scale --out path/to/out.csv
# =============================================================================
run_union_task_matrix_cli <- function() {
  ca <- commandArgs(trailingOnly = TRUE)
  if (!("--union-task-matrix" %in% ca)) return(FALSE)

  parse <- function(flag) {
    i <- which(ca == flag)
    if (length(i) == 0) return(NULL)
    ca[i + 1L]
  }
  repo <- parse("--repo") %||% dirname(dirname(dirname(getwd())))
  matrix_path <- parse("--matrix")
  meta_path <- parse("--meta")
  mode <- parse("--mode") %||% "shift_only"
  out_path <- parse("--out")
  if (is.null(matrix_path) || is.null(meta_path) || is.null(out_path))
    stop("union-task-matrix requires --matrix, --meta, --out")

  BRIDGE_DIR <<- file.path(repo, "reports/benchmark_master/methods/bridge_aware")
  cat("Union-task bridge:", mode, "\n  matrix:", matrix_path, "\n  out:", out_path, "\n")

  cptac_bridge <- fread(file.path(BRIDGE_DIR, "bridge_summary_cptac.tsv"))
  ccle_bridge <- fread(file.path(BRIDGE_DIR, "bridge_summary_ccle.tsv"))

  dt <- fread(matrix_path, header = TRUE)
  id_cols <- intersect(c("GeneSymbol", "UniProtID", "Gene"), names(dt))
  scols <- setdiff(names(dt), id_cols)
  mat <- as.matrix(dt[, ..scols])
  rownames(mat) <- dt[[names(dt)[1]]]

  ann <- fread(meta_path)
  stopifnot(all(c("sample_id", "domain") %in% names(ann)))
  c_sid <- ann[toupper(domain) == "CPTAC", sample_id]
  e_sid <- ann[toupper(domain) == "CCLE", sample_id]
  c_cols <- intersect(c_sid, colnames(mat))
  e_cols <- intersect(e_sid, colnames(mat))
  cptac_mat <- mat[, c_cols, drop = FALSE]
  ccle_mat <- mat[, e_cols, drop = FALSE]

  res <- run_bridge_aware_representation(
    cptac_mat, ccle_mat, cptac_bridge, ccle_bridge,
    config = list(mode = mode)
  )
  out_dt <- data.table(GeneSymbol = rownames(res$matrix), as.data.table(res$matrix))
  fwrite(out_dt, out_path)
  cat("Wrote:", out_path, " (", nrow(res$matrix), "x", ncol(res$matrix), ")\n")
  TRUE
}

if ("--union-task-matrix" %in% commandArgs(trailingOnly = TRUE)) {
  ok <- run_union_task_matrix_cli()
  if (!ok) stop("union-task-matrix CLI failed")
} else {
  main()
}
