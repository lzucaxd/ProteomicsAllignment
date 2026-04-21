#!/usr/bin/env Rscript
# Step 9: Biology destruction (vs raw CPTAC), marker sanity, residual dependence.

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
results_root <- "reports/benchmark_master/benchmark_results"
methods_root <- "data/processed/methods"
meta_dir <- "data/processed/union"

for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--results-root" && i < length(args)) results_root <- args[i + 1]
  if (args[i] == "--methods-root" && i < length(args)) methods_root <- args[i + 1]
  if (args[i] == "--meta-dir" && i < length(args)) meta_dir <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

dest_r <- file.path(repo_root, "src/harmonize/benchmark/calibration/biology_destruction.R")
mark_r <- file.path(repo_root, "src/harmonize/benchmark/calibration/marker_sanity.R")
resid_r <- file.path(repo_root, "src/harmonize/benchmark/calibration/residual_dependence.R")

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

for (tk in tasks) {
  native_da <- file.path(repo_root, results_root, "raw", tk$name,
                         "representation_da", "cptac", "da_limma_result.csv")
  if (!file.exists(native_da)) stop("Missing native CPTAC DA: ", native_da)

  for (m in methods) {
    outcal <- file.path(repo_root, results_root, m, tk$name, "calibration")
    dir.create(outcal, recursive = TRUE, showWarnings = FALSE)
    meth_da <- file.path(repo_root, results_root, m, tk$name,
                         "representation_da", "cptac", "da_limma_result.csv")
    matf <- file.path(repo_root, methods_root, m, paste0("transformed_", tk$name, ".csv"))
    metf <- file.path(repo_root, meta_dir, tk$meta)
    if (!file.exists(metf))
      metf <- file.path(repo_root, "data/processed/union", tk$meta)
    if (!file.exists(metf))
      metf <- file.path(repo_root, "data/processed_union", tk$meta)

    if (m != "raw" && file.exists(meth_da)) {
      cat("\n=== Biology destruction CPTAC:", m, tk$name, "===\n")
      st <- system2("Rscript", c(dest_r,
        "--native-da", native_da,
        "--method-da", meth_da,
        "--method", m,
        "--outdir", outcal), stdout = "", stderr = "")
      if (!identical(as.integer(st), 0L)) stop("biology_destruction failed")
    }

    for (dom in c("CPTAC", "CCLE")) {
      da_dom <- file.path(repo_root, results_root, m, tk$name,
                          "representation_da", tolower(dom), "da_limma_result.csv")
      if (!file.exists(da_dom)) next
      cat("\n=== Marker sanity:", m, tk$name, dom, "===\n")
      st <- system2("Rscript", c(mark_r,
        "--da-result", da_dom,
        "--markers", tk$mk,
        "--expected-signs", tk$es,
        "--method", m,
        "--domain", dom,
        "--task", tk$name,
        "--outdir", outcal), stdout = "", stderr = "")
      if (!identical(as.integer(st), 0L)) stop("marker_sanity failed")
    }

    if (file.exists(matf) && file.exists(metf)) {
      for (dom in c("CPTAC", "CCLE")) {
        cat("\n=== Residual dependence:", m, tk$name, dom, "===\n")
        st <- system2("Rscript", c(resid_r,
          "--matrix", matf,
          "--meta", metf,
          "--contrast-a", tk$ca,
          "--contrast-b", tk$cb,
          "--domain", dom,
          "--outdir", outcal), stdout = "", stderr = "")
        if (!identical(as.integer(st), 0L))
          warning("residual_dependence failed (small n?) ", m, " ", tk$name, " ", dom)
      }
    }
  }
}

cat("\nStep 9 complete.\n")
