#!/usr/bin/env Rscript
# Step 3: limma DA for each method × task (CPTAC + CCLE per wrapper call).

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
methods_root <- "data/processed/methods"
meta_dir <- "data/processed/union"
results_root <- "reports/benchmark_master/benchmark_results"
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--methods-root" && i < length(args)) methods_root <- args[i + 1]
  if (args[i] == "--meta-dir" && i < length(args)) meta_dir <- args[i + 1]
  if (args[i] == "--results-root" && i < length(args)) results_root <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)
wrapper <- file.path(repo_root, "src/harmonize/benchmark/calibration/limma_da_wrapper.R")

tasks <- list(
  list(name = "breast_subtype", meta = "sample_meta_breast_subtype.csv",
       ca = "Basal", cb = "Luminal", cname = "Luminal_vs_Basal"),
  list(name = "breast_vs_lung", meta = "sample_meta_breast_vs_lung.csv",
       ca = "Breast", cb = "Lung", cname = "Lung_vs_Breast")
)
methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")

for (m in methods) {
  for (tk in tasks) {
    matf <- file.path(repo_root, methods_root, m, paste0("transformed_", tk$name, ".csv"))
    metf <- file.path(repo_root, meta_dir, tk$meta)
    if (!file.exists(metf))
      metf <- file.path(repo_root, "data/processed/union", tk$meta)
    if (!file.exists(metf))
      metf <- file.path(repo_root, "data/processed_union", tk$meta)
    if (!file.exists(matf)) {
      warning("Missing matrix, skip: ", matf)
      next
    }
    outd <- file.path(repo_root, results_root, m, tk$name, "representation_da")
    dir.create(outd, recursive = TRUE, showWarnings = FALSE)
    cat("\n=== limma:", m, tk$name, "===\n")
    status <- system2(
      "Rscript",
      c(wrapper,
        "--matrix", matf,
        "--meta", metf,
        "--contrast-a", tk$ca,
        "--contrast-b", tk$cb,
        "--contrast-name", tk$cname,
        "--outdir", outd),
      stdout = "", stderr = ""
    )
    if (!identical(as.integer(status), 0L)) stop("limma failed: ", m, " ", tk$name)
  }
}
cat("\nStep 3 complete.\n")
