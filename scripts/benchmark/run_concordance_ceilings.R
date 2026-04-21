#!/usr/bin/env Rscript
# Step 8: Concordance ceiling per method × task × domain (CPTAC + optional CCLE split-half).

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
methods_root <- "data/processed/methods"
meta_dir <- "data/processed/union"
results_root <- "reports/benchmark_master/benchmark_results"
inter_dir <- "data/processed"
n_splits <- 200L
seed <- 42L
ccle_split_half <- FALSE

for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--methods-root" && i < length(args)) methods_root <- args[i + 1]
  if (args[i] == "--meta-dir" && i < length(args)) meta_dir <- args[i + 1]
  if (args[i] == "--results-root" && i < length(args)) results_root <- args[i + 1]
  if (args[i] == "--intersection-dir" && i < length(args)) inter_dir <- args[i + 1]
  if (args[i] == "--n-splits" && i < length(args)) n_splits <- as.integer(args[i + 1])
  if (args[i] == "--seed" && i < length(args)) seed <- as.integer(args[i + 1])
  if (args[i] == "--ccle-split-half") ccle_split_half <- TRUE
}

repo_root <- normalizePath(repo_root, mustWork = TRUE)
ceil_r <- file.path(repo_root, "src/harmonize/benchmark/calibration/concordance_ceiling.R")

tasks <- list(
  list(name = "breast_subtype", meta = "sample_meta_breast_subtype.csv",
       ca = "Basal", cb = "Luminal"),
  list(name = "breast_vs_lung", meta = "sample_meta_breast_vs_lung.csv",
       ca = "Breast", cb = "Lung")
)
methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")

run_ceiling <- function(matf, metf, inter, outd, domain, ca, cb, force_split) {
  cmd <- c(
    ceil_r,
    "--matrix", matf,
    "--meta", metf,
    "--domain", domain,
    "--contrast-a", ca,
    "--contrast-b", cb,
    "--outdir", outd,
    "--n-splits", as.character(n_splits),
    "--seed", as.character(seed)
  )
  if (file.exists(inter)) cmd <- c(cmd, "--intersection-genes-file", inter)
  if (isTRUE(force_split)) cmd <- c(cmd, "--force-split-half")
  status <- system2("Rscript", cmd, stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L))
    stop("concordance_ceiling failed: ", domain)
}

for (tk in tasks) {
  for (m in methods) {
    matf <- file.path(repo_root, methods_root, m, paste0("transformed_", tk$name, ".csv"))
    metf <- file.path(repo_root, meta_dir, tk$meta)
    if (!file.exists(metf))
      metf <- file.path(repo_root, "data/processed_union", tk$meta)
    inter <- file.path(repo_root, inter_dir, paste0("intersection_genes_", tk$name, ".txt"))
    if (!file.exists(inter))
      inter <- file.path(repo_root, "data/processed_union", paste0("intersection_genes_", tk$name, ".txt"))
    outd <- file.path(repo_root, results_root, m, tk$name, "calibration")
    dir.create(outd, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(matf)) {
      warning("skip missing ", matf)
      next
    }

    cat("\n=== Concordance ceiling CPTAC:", m, tk$name, "===\n")
    run_ceiling(matf, metf, inter, outd, "CPTAC", tk$ca, tk$cb, FALSE)

    if (isTRUE(ccle_split_half)) {
      cat("=== Concordance ceiling CCLE (split-half):", m, tk$name, "===\n")
      run_ceiling(matf, metf, inter, outd, "CCLE", tk$ca, tk$cb, TRUE)
    }
  }
}

cat("\nStep 8 complete.\n")
