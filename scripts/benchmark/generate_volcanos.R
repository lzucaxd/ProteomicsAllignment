#!/usr/bin/env Rscript
# Step 6: Volcano plots (ggplot2) — raw method, both tasks, CPTAC + CCLE.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
results_root <- "reports/benchmark_master/benchmark_results"
diag_dir <- "reports/benchmark_master/diagnostics"
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--results-root" && i < length(args)) results_root <- args[i + 1]
  if (args[i] == "--diag-dir" && i < length(args)) diag_dir <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)
diag_dir <- normalizePath(file.path(repo_root, diag_dir), mustWork = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

generate_volcano <- function(da_path, title, output_path) {
  da <- fread(da_path)
  n_sig <- sum(da$adj.P.Val < 0.05, na.rm = TRUE)

  p <- ggplot(da, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(alpha = 0.15, size = 0.5, color = "grey50") +
    geom_point(data = da[da$adj.P.Val < 0.05], color = "steelblue", alpha = 0.4, size = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(
      title = title,
      subtitle = paste0(n_sig, " / ", nrow(da), " significant (adj.p < 0.05)")
    ) +
    theme_minimal(base_size = 14)

  ggsave(output_path, p, width = 7, height = 5)
  cat("Saved:", output_path, "\n")
}

tasks <- c("breast_subtype", "breast_vs_lung")
domains <- c("cptac", "ccle")

for (task in tasks) {
  for (dom in domains) {
    da_path <- file.path(repo_root, results_root, "raw", task, "representation_da", dom, "da_limma_result.csv")
    if (!file.exists(da_path)) {
      warning("Missing ", da_path)
      next
    }
    out_pdf <- file.path(diag_dir, paste0("volcano_", dom, "_", task, ".pdf"))
    generate_volcano(da_path, paste0("Raw — ", task, " — ", toupper(dom)), out_pdf)
  }
}
cat("Step 6 complete.\n")
