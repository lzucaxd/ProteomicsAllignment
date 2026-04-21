#!/usr/bin/env Rscript
# Step 5: Stratified FC by significance stratum (all methods × tasks).

suppressPackageStartupMessages(library(data.table))

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

methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
tasks <- c("breast_subtype", "breast_vs_lung")

stratify_one <- function(merged) {
  merged[, stratum := fcase(
    adj.P.Val_cptac < 0.05 & adj.P.Val_ccle < 0.05, "sig_both",
    adj.P.Val_cptac < 0.05, "sig_cptac_only",
    adj.P.Val_ccle < 0.05, "sig_ccle_only",
    default = "sig_neither"
  )]
  out <- merged[, {
    list(
      n_genes = .N,
      fc_correlation = if (.N > 2) cor(logFC_cptac, logFC_ccle, use = "complete.obs") else NA_real_,
      same_dir_fraction = mean(sign(logFC_cptac) == sign(logFC_ccle), na.rm = TRUE)
    )
  }, by = stratum]
  out
}

all_out <- list()

for (task in tasks) {
  task_rows <- list()
  for (m in methods) {
    fc_path <- file.path(repo_root, results_root, m, task, "representation_da", "fc_agreement.csv")
    if (!file.exists(fc_path)) next
    ag <- fread(fc_path)
    if (!all(c("logFC_cptac", "logFC_ccle", "adj.P.Val_cptac", "adj.P.Val_ccle") %in% names(ag))) {
      warning("fc_agreement missing adj.P columns for ", m, " ", task, " — re-run Step 4")
      next
    }
    st <- stratify_one(ag)
    st[, `:=`(method = m, task = task)]
    task_rows[[length(task_rows) + 1]] <- st
  }
  if (length(task_rows)) {
    dt <- rbindlist(task_rows)
    setcolorder(dt, c("method", "task", "stratum", setdiff(names(dt), c("method", "task", "stratum"))))
    all_out[[length(all_out) + 1]] <- dt
    outf <- file.path(diag_dir, paste0("fc_stratified_", task, ".csv"))
    fwrite(dt, outf)
    cat("\n=== Stratified FC:", task, "===\n")
    print(dt)
    cat("Saved:", outf, "\n")
  }
}

if (length(all_out)) {
  fwrite(rbindlist(all_out), file.path(diag_dir, "fc_stratified_all_tasks.csv"))
}
cat("\nStep 5 complete.\n")
