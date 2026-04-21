#!/usr/bin/env Rscript
# =============================================================================
# Diagnose Subtype Sign — standalone contrast-direction check on raw data
# =============================================================================
# Usage:
#   Rscript scripts/benchmark/diagnose_subtype_sign.R [--repo-root .]
#   Rscript ... [--matrix path/to/shared_gene_matrix_breast_subtype.csv] [--meta path/to/sample_meta_breast_subtype.csv]
#
# Loads the RAW shared gene matrix and sample metadata for breast_subtype,
# runs validate_contrast_direction() for CPTAC and CCLE, and prints a
# clear diagnostic before any downstream calibration can proceed.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ── Parse arguments ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
matrix_path_cli <- NULL
meta_path_cli <- NULL
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--matrix" && i < length(args)) matrix_path_cli <- args[i + 1]
  if (args[i] == "--meta" && i < length(args)) meta_path_cli <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

# ── Source helpers ───────────────────────────────────────────────────────────
source(file.path(repo_root, "src/harmonize/benchmark/calibration/contrast_validation.R"))
source(file.path(repo_root, "scripts/benchmark/evaluation_helpers.R"))

# ── Load raw shared gene matrix and metadata ─────────────────────────────────
matrix_path <- if (!is.null(matrix_path_cli)) {
  normalizePath(matrix_path_cli, mustWork = TRUE)
} else {
  file.path(repo_root, "data/processed/shared_gene_matrix_breast_subtype.csv")
}
meta_path <- if (!is.null(meta_path_cli)) {
  normalizePath(meta_path_cli, mustWork = TRUE)
} else {
  file.path(repo_root, "data/processed/sample_meta_breast_subtype.csv")
}

if (!file.exists(matrix_path)) stop("Raw matrix not found: ", matrix_path)
if (!file.exists(meta_path))   stop("Sample metadata not found: ", meta_path)

mat_raw <- as.matrix(fread(matrix_path), rownames = 1)
meta    <- fread(meta_path)

cat("\n", strrep("=", 60), "\n")
cat("  CONTRAST DIRECTION DIAGNOSTIC — breast_subtype\n")
cat(strrep("=", 60), "\n\n")

cat("Matrix:", nrow(mat_raw), "genes x", ncol(mat_raw), "samples\n")
cat("Metadata:", nrow(meta), "rows\n\n")

# ── Expected marker directions for Luminal_vs_Basal ──────────────────────────
# Convention: contrast = Luminal - Basal
#   positive logFC = higher in Luminal
expected_dirs <- data.table(
  gene = c("ESR1", "PGR", "GATA3", "FOXA1", "EGFR", "KRT5", "KRT17", "FOXC1"),
  expected_sign = c(1, 1, 1, 1, -1, -1, -1, -1)
)

# ── Run per-domain validation ────────────────────────────────────────────────
diag_results <- list()

for (dom in c("CPTAC", "CCLE")) {
  cat(strrep("-", 50), "\n")
  cat("Domain:", dom, "\n")
  cat(strrep("-", 50), "\n")

  dom_sids <- meta[toupper(domain) == dom, sample_id]
  dom_cols <- intersect(dom_sids, colnames(mat_raw))

  if (length(dom_cols) < 4) {
    cat("  SKIP: only", length(dom_cols), "samples\n\n")
    next
  }

  dom_mat    <- mat_raw[, dom_cols, drop = FALSE]
  dom_groups <- meta[match(dom_cols, sample_id), condition]

  if (length(unique(dom_groups)) < 2) {
    cat("  SKIP: only 1 condition level\n\n")
    next
  }

  # Filter genes with >50% NA in either group
  for (g in unique(dom_groups)) {
    g_cols <- which(dom_groups == g)
    na_frac <- rowMeans(is.na(dom_mat[, g_cols, drop = FALSE]))
    dom_mat[na_frac > 0.5, g_cols] <- NA
  }
  keep <- rowSums(!is.na(dom_mat)) >= ncol(dom_mat) * 0.5
  dom_mat <- dom_mat[keep, , drop = FALSE]

  cat("  Samples:", ncol(dom_mat), "  Genes:", nrow(dom_mat), "\n")
  cat("  Groups:", paste(table(dom_groups), collapse = " / "),
      "(", paste(names(table(dom_groups)), collapse = " / "), ")\n")

  val <- validate_contrast_direction(
    dom_mat, dom_groups, expected_dirs,
    contrast_label = paste0("Luminal_vs_Basal_", dom)
  )

  # Print marker table
  cat("\n  Marker Direction Table:\n")
  mt <- val$marker_table
  for (r in seq_len(nrow(mt))) {
    flag <- if (mt$correct[r]) "  OK" else "  XX"
    cat(sprintf("    %-8s  expected=%+d  observed_logFC=%+7.3f  sign=%+d  %s\n",
                mt$gene[r], mt$expected_sign[r], mt$logFC[r],
                mt$observed_sign[r], flag))
  }

  cat(sprintf("\n  Summary: %d / %d correct (%.1f%%)\n",
              val$summary$n_correct, val$summary$n_markers_tested,
              val$summary$fraction_correct * 100))
  cat("  Contrast direction:", val$summary$contrast_direction, "\n")
  cat("  Likely flipped:", val$likely_flipped, "\n\n")

  diag_results[[dom]] <- val
}

# ── Overall verdict ──────────────────────────────────────────────────────────
cat(strrep("=", 60), "\n")
cat("  OVERALL DIAGNOSTIC\n")
cat(strrep("=", 60), "\n\n")

any_flipped <- any(sapply(diag_results, function(x) isTRUE(x$likely_flipped)))

if (any_flipped) {
  cat("  *** LIKELY FLIPPED — STOP ***\n")
  cat("  At least one domain shows majority-incorrect marker directions.\n")
  cat("  Check: (1) contrast convention, (2) metadata condition labels,\n")
  cat("         (3) potential confounders, (4) subtype annotation source.\n")
  cat("  Do NOT proceed with downstream calibration until resolved.\n\n")
} else if (length(diag_results) == 0) {
  cat("  *** NO DOMAINS TESTED — check data paths ***\n\n")
} else {
  cat("  Contrast OK — marker directions are consistent in all domains.\n")
  cat("  Safe to proceed with downstream calibration.\n\n")
}

# ── Save results ─────────────────────────────────────────────────────────────
outdir <- file.path(repo_root, "reports/benchmark_master/diagnostics")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

all_markers <- rbindlist(lapply(names(diag_results), function(dom) {
  dt <- diag_results[[dom]]$marker_table
  dt$domain <- dom
  dt
}), fill = TRUE)

if (nrow(all_markers) > 0) {
  fwrite(all_markers, file.path(outdir, "subtype_sign_diagnostic.csv"))
  cat("Saved:", file.path(outdir, "subtype_sign_diagnostic.csv"), "\n")
}

summary_dt <- rbindlist(lapply(names(diag_results), function(dom) {
  s <- diag_results[[dom]]$summary
  data.table(
    domain = dom,
    contrast_label = s$contrast_label,
    contrast_direction = s$contrast_direction,
    n_markers_tested = s$n_markers_tested,
    n_correct = s$n_correct,
    fraction_correct = s$fraction_correct,
    likely_flipped = diag_results[[dom]]$likely_flipped
  )
}))

if (nrow(summary_dt) > 0) {
  fwrite(summary_dt, file.path(outdir, "subtype_sign_diagnostic_summary.csv"))
  cat("Saved:", file.path(outdir, "subtype_sign_diagnostic_summary.csv"), "\n")
}

cat("\nDone.\n")
