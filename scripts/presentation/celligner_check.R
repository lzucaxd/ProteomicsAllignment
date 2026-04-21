#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1]))) else normalizePath(file.path(getwd(), "scripts", "presentation"))
source(file.path(sd, "presentation_paths.R"))

cat("\n=== CELLIGNER STATUS CHECK ===\n")

log_dir <- file.path(REPO, "reports/benchmark_master/logs")
if (dir.exists(log_dir)) {
  log_files <- sort(Sys.glob(file.path(log_dir, "overnight_v2_*.log")))
  if (length(log_files) > 0L) {
    latest_log <- log_files[length(log_files)]
    log_content <- readLines(latest_log, warn = FALSE)
    celligner_lines <- grep("[Cc]elligner", log_content, value = TRUE)
    implemented <- grep("FULLY_IMPLEMENTED", log_content, value = TRUE)
    scaffold <- grep("[Ss]caffold|SCAFFOLDED", log_content, value = TRUE)

    cat("Latest log:", latest_log, "\n")
    if (length(implemented) > 0L) {
      cat("STATUS: FULLY_IMPLEMENTED (lines found)\n")
      cat("Safe to present Celligner results.\n")
    } else if (length(scaffold) > 0L) {
      cat("STATUS: SCAFFOLD MODE\n")
      cat("WARNING: Celligner may be placeholder — verify matrices.\n")
    } else {
      cat("STATUS: UNCLEAR — inspect log\n")
      if (length(celligner_lines) > 0L) cat(paste(celligner_lines, collapse = "\n"), "\n")
    }
  } else {
    cat("No overnight_v2_*.log in ", log_dir, "\n", sep = "")
  }
} else {
  cat("Log directory missing: ", log_dir, "\n", sep = "")
}

cell_da <- file.path(REPO, "reports/benchmark_master/benchmark_results/celligner/breast_subtype/representation_da/cptac/da_limma_result.csv")
raw_da <- file.path(REPO, "reports/benchmark_master/benchmark_results/raw/breast_subtype/representation_da/cptac/da_limma_result.csv")

if (file.exists(cell_da) && file.exists(raw_da)) {
  c_da <- fread(cell_da)
  r_da <- fread(raw_da)
  merged <- merge(c_da, r_da, by = "gene", suffixes = c("_cell", "_raw"))
  cor_val <- cor(merged$logFC_cell, merged$logFC_raw, use = "complete.obs")
  cat(sprintf("\nCelligner vs Raw logFC correlation (CPTAC, subtype): %.4f\n", cor_val))
  if (!is.na(cor_val) && cor_val > 0.999) {
    cat("WARNING: logFCs nearly identical → possible scaffold / no real transform.\n")
  } else {
    cat("logFCs differ → transformation likely changed values.\n")
  }
} else {
  cat("\nMissing celligner or raw CPTAC da_limma_result for subtype — skip correlation check.\n")
}
