#!/usr/bin/env Rscript
# Step 10 (v2): Disconnect scores from raw vs method structure + intersection FC + CPTAC ceiling.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
results_root <- "reports/benchmark_master/benchmark_results"
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--results-root" && i < length(args)) results_root <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

tasks <- c("breast_subtype", "breast_vs_lung")
methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
out_csv <- file.path(repo_root, results_root, "disconnect_scores.csv")

rows <- list()

for (task in tasks) {
  sraw <- file.path(repo_root, results_root, "raw", task, "structure", "structure_summary.csv")
  if (!file.exists(sraw)) next
  raw_dt <- fread(sraw)
  raw_dr2 <- raw_dt$domain_r2_pc1[1]

  cm_raw_path <- file.path(repo_root, results_root, "raw", task, "representation_da", "cross_domain_metrics.csv")
  if (!file.exists(cm_raw_path)) next
  cm_raw <- fread(cm_raw_path)
  raw_fc <- cm_raw[gene_set == "intersection"]$fc_correlation[1]

  ceil_p <- file.path(repo_root, results_root, "raw", task, "calibration", "ceiling_summary_cptac.csv")
  if (!file.exists(ceil_p)) {
    alt <- file.path(repo_root, results_root, "raw", task, "calibration", "ceiling_summary.csv")
    ceil_p <- if (file.exists(alt)) alt else NA_character_
  }
  if (is.na(ceil_p) || !file.exists(ceil_p)) next
  ceilv <- fread(ceil_p)$ceiling_fc_correlation[1]
  denom <- ceilv - raw_fc
  if (!is.finite(denom) || abs(denom) < 1e-8) denom <- NA_real_

  for (m in methods) {
    if (m == "raw") next
    sm <- file.path(repo_root, results_root, m, task, "structure", "structure_summary.csv")
    cm_m <- file.path(repo_root, results_root, m, task, "representation_da", "cross_domain_metrics.csv")
    dest <- file.path(repo_root, results_root, m, task, "calibration",
                      paste0("destruction_summary_", m, ".csv"))
    if (!file.exists(sm) || !file.exists(cm_m)) next
    m_dr2 <- fread(sm)$domain_r2_pc1[1]
    m_fc <- fread(cm_m)[gene_set == "intersection"]$fc_correlation[1]
    ret <- NA_real_
    if (file.exists(dest)) {
      d <- fread(dest)
      if ("default_retention_rate" %in% names(d)) ret <- d$default_retention_rate[1]
    }
    geom_i <- if (is.finite(raw_dr2) && raw_dr2 != 0) (raw_dr2 - m_dr2) / raw_dr2 else NA_real_
    da_i <- if (is.finite(denom) && is.finite(raw_fc) && is.finite(m_fc)) (m_fc - raw_fc) / denom else NA_real_
    bio_cost <- if (is.finite(ret)) 1 - ret else NA_real_
    rows[[length(rows) + 1]] <- data.table(
      method = m, task = task,
      geom_improvement = geom_i,
      da_improvement = da_i,
      disconnect = geom_i - da_i,
      disconnect_score = geom_i - da_i,
      biology_cost = bio_cost,
      raw_domain_r2_pc1 = raw_dr2,
      method_domain_r2_pc1 = m_dr2,
      raw_fc_intersection = raw_fc,
      method_fc_intersection = m_fc,
      ceiling_cptac_fc = ceilv
    )
  }
}

if (length(rows)) {
  dt <- rbindlist(rows)
  fwrite(dt, out_csv)
  cat("Wrote:", out_csv, "\n")
  print(dt)
} else {
  cat("No disconnect rows computed (missing structure or cross_domain inputs).\n")
}
