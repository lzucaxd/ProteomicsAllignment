#!/usr/bin/env Rscript
# Step 1: Write intersection gene lists from Step 0a audit CSVs.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
union_copy_dir <- NA_character_
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--intersection-out-dir" && i < length(args)) union_copy_dir <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

diag_dir <- file.path(repo_root, "reports/benchmark_master/diagnostics")
out_dir <- file.path(repo_root, "data/processed_union")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
proc_dir <- file.path(repo_root, "data", "processed")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

for (task in c("breast_subtype", "breast_vs_lung")) {
  audit_path <- file.path(diag_dir, paste0("gene_coverage_audit_", task, ".csv"))
  if (!file.exists(audit_path)) {
    warning("Missing audit (run preflight first): ", audit_path)
    next
  }
  audit <- fread(audit_path)
  genes <- audit[category == "both_domains", gene]
  out_file <- file.path(out_dir, paste0("intersection_genes_", task, ".txt"))
  writeLines(genes, out_file)
  cat(task, ": intersection genes =", length(genes), "->", out_file, "\n")
  proc_file <- file.path(proc_dir, paste0("intersection_genes_", task, ".txt"))
  writeLines(genes, proc_file)
  cat("  also:", proc_file, "\n")
  if (!is.na(union_copy_dir) && nzchar(union_copy_dir)) {
    alt <- file.path(repo_root, union_copy_dir, paste0("intersection_genes_", task, ".txt"))
    dir.create(dirname(alt), recursive = TRUE, showWarnings = FALSE)
    writeLines(genes, alt)
    cat("  also:", alt, "\n")
  }
}
