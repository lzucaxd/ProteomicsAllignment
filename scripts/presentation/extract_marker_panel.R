#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1]))) else normalizePath(file.path(getwd(), "scripts", "presentation"))
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

extract_markers <- function(task, markers_expected, repo) {
  da_cptac <- fread(file.path(repo, sprintf(
    "reports/benchmark_master/benchmark_results/raw/%s/representation_da/cptac/da_limma_result.csv", task
  )))
  da_ccle <- fread(file.path(repo, sprintf(
    "reports/benchmark_master/benchmark_results/raw/%s/representation_da/ccle/da_limma_result.csv", task
  )))

  results <- data.frame(
    gene = markers_expected$gene,
    expected_direction = markers_expected$direction,
    stringsAsFactors = FALSE
  )
  results$cptac_logFC <- NA_real_
  results$cptac_pval <- NA_real_
  results$cptac_sig <- NA
  results$ccle_logFC <- NA_real_
  results$ccle_pval <- NA_real_
  results$ccle_sig <- NA
  results$cptac_correct <- NA
  results$ccle_correct <- NA

  for (i in seq_len(nrow(results))) {
    g <- results$gene[i]
    row_c <- da_cptac[gene == g]
    if (nrow(row_c) > 0L) {
      results$cptac_logFC[i] <- round(row_c$logFC[1L], 3L)
      results$cptac_pval[i] <- signif(row_c$adj.P.Val[1L], 3L)
      results$cptac_sig[i] <- row_c$adj.P.Val[1L] < 0.05
    }
    row_e <- da_ccle[gene == g]
    if (nrow(row_e) > 0L) {
      results$ccle_logFC[i] <- round(row_e$logFC[1L], 3L)
      results$ccle_pval[i] <- signif(row_e$adj.P.Val[1L], 3L)
      results$ccle_sig[i] <- row_e$adj.P.Val[1L] < 0.05
    }
    exp_dir <- markers_expected$direction[i]
    if (!is.na(results$cptac_logFC[i])) {
      results$cptac_correct[i] <- if (identical(exp_dir, "up_luminal")) {
        results$cptac_logFC[i] > 0
      } else {
        results$cptac_logFC[i] < 0
      }
    }
    if (!is.na(results$ccle_logFC[i])) {
      results$ccle_correct[i] <- if (identical(exp_dir, "up_luminal")) {
        results$ccle_logFC[i] > 0
      } else {
        results$ccle_logFC[i] < 0
      }
    }
  }
  results
}

subtype_markers <- data.frame(
  gene = c("ESR1", "PGR", "GATA3", "FOXA1", "EGFR", "KRT5", "KRT17", "FOXC1"),
  direction = c(rep("up_luminal", 4L), rep("up_basal", 4L)),
  stringsAsFactors = FALSE
)

panel <- extract_markers("breast_subtype", subtype_markers, REPO)
cat("\n=== MARKER PANEL: BREAST SUBTYPE (raw representation DA) ===\n")
print(panel)

fwrite(panel, file.path(PRES_OUT, "tables/marker_panel_subtype.csv"))

cat(sprintf(
  "\nCPTAC: %d/%d correct direction\n",
  sum(panel$cptac_correct, na.rm = TRUE),
  sum(!is.na(panel$cptac_correct))
))
cat(sprintf(
  "CCLE: %d/%d correct direction\n",
  sum(panel$ccle_correct, na.rm = TRUE),
  sum(!is.na(panel$ccle_correct))
))
cat(sprintf(
  "CPTAC: %d/%d significant (FDR<0.05)\n",
  sum(panel$cptac_sig, na.rm = TRUE),
  sum(!is.na(panel$cptac_sig))
))
cat(sprintf(
  "CCLE: %d/%d significant (FDR<0.05)\n",
  sum(panel$ccle_sig, na.rm = TRUE),
  sum(!is.na(panel$ccle_sig))
))
