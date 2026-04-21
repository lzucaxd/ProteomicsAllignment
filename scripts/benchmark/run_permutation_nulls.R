#!/usr/bin/env Rscript
# Step 7: Per-method permutation nulls (1000 perm); intersection genes for FC stats.
# Reuses raw null for bridge_shift (identical within-domain logFC).

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
methods_root <- "data/processed/methods"
meta_dir <- "data/processed/union"
results_root <- "reports/benchmark_master/benchmark_results"
inter_dir <- "data/processed"
n_perm <- 1000L
seed <- 42L

parse_flag <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1L]
}
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--methods-root" && i < length(args)) methods_root <- args[i + 1]
  if (args[i] == "--meta-dir" && i < length(args)) meta_dir <- args[i + 1]
  if (args[i] == "--results-root" && i < length(args)) results_root <- args[i + 1]
  if (args[i] == "--intersection-dir" && i < length(args)) inter_dir <- args[i + 1]
}
if (!is.null(parse_flag("--n-perm"))) n_perm <- as.integer(parse_flag("--n-perm"))
if (!is.null(parse_flag("--seed"))) seed <- as.integer(parse_flag("--seed"))

repo_root <- normalizePath(repo_root, mustWork = TRUE)
perm_r <- file.path(repo_root, "src/harmonize/benchmark/calibration/permutation_null.R")

tasks <- list(
  list(name = "breast_subtype", meta = "sample_meta_breast_subtype.csv",
       ca = "Basal", cb = "Luminal",
       mk = "ESR1,PGR,GATA3,FOXA1,EGFR,KRT5,KRT17,FOXC1",
       es = "1,1,1,1,-1,-1,-1,-1"),
  list(name = "breast_vs_lung", meta = "sample_meta_breast_vs_lung.csv",
       ca = "Breast", cb = "Lung",
       mk = "GATA3,FOXA1,ESR1,NKX2-1,SFTPB,NAPSA",
       es = "-1,-1,-1,1,1,1")
)

methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")

run_one <- function(method, tk) {
  matf <- file.path(repo_root, methods_root, method, paste0("transformed_", tk$name, ".csv"))
  metf <- file.path(repo_root, meta_dir, tk$meta)
  if (!file.exists(metf))
    metf <- file.path(repo_root, "data/processed_union", tk$meta)
  inter <- file.path(repo_root, inter_dir, paste0("intersection_genes_", tk$name, ".txt"))
  if (!file.exists(inter))
    inter <- file.path(repo_root, "data/processed_union", paste0("intersection_genes_", tk$name, ".txt"))
  outd <- file.path(repo_root, results_root, method, tk$name, "calibration")
  dir.create(outd, recursive = TRUE, showWarnings = FALSE)

  if (method == "bridge_shift") {
    src <- file.path(repo_root, results_root, "raw", tk$name, "calibration")
    for (fn in c("null_distribution.csv", "observed_metrics.csv", "observed_vs_null_summary.csv")) {
      sf <- file.path(src, fn)
      if (!file.exists(sf)) stop("bridge_shift reuse: missing ", sf)
      file.copy(sf, file.path(outd, fn), overwrite = TRUE)
    }
    cat("Reused raw null -> ", outd, "\n")
    return(invisible())
  }

  if (!file.exists(matf)) stop("Missing matrix ", matf)

  cmd <- c(
    perm_r,
    "--matrix", matf,
    "--meta", metf,
    "--contrast-a", tk$ca,
    "--contrast-b", tk$cb,
    "--outdir", outd,
    "--n-perm", as.character(n_perm),
    "--seed", as.character(seed),
    "--markers", tk$mk,
    "--expected-signs", tk$es
  )
  if (file.exists(inter)) cmd <- c(cmd, "--intersection-genes-file", inter)

  cat("\n=== Permutation null:", method, tk$name, "===\n")
  status <- system2("Rscript", cmd, stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L)) stop("permutation_null failed: ", method, " ", tk$name)
}

for (tk in tasks) {
  for (m in methods) run_one(m, tk)
}

cat("\nStep 7 complete (n_perm=", n_perm, ", seed=", seed, ").\n", sep = "")
