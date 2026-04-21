#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1]))) else normalizePath(file.path(getwd(), "scripts", "presentation"))
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

compute_fc_se_summary <- function(da_path, method, task, domain) {
  da <- fread(da_path)
  da[, SE := ifelse(!is.na(t) & t != 0, abs(logFC / t), NA_real_)]
  sig <- da[!is.na(adj.P.Val) & adj.P.Val < 0.05]

  data.table(
    method = method,
    task = task,
    domain = domain,
    n_tested = nrow(da),
    n_significant = nrow(sig),
    mean_abs_logFC_all = round(mean(abs(da$logFC), na.rm = TRUE), 4L),
    median_abs_logFC_all = round(median(abs(da$logFC), na.rm = TRUE), 4L),
    median_SE_all = round(median(da$SE, na.rm = TRUE), 4L),
    mean_abs_logFC_sig = round(mean(abs(sig$logFC), na.rm = TRUE), 4L),
    median_abs_logFC_sig = round(median(abs(sig$logFC), na.rm = TRUE), 4L),
    median_SE_sig = round(median(sig$SE, na.rm = TRUE), 4L),
    iqr_logFC_sig = round(stats::IQR(sig$logFC, na.rm = TRUE), 4L),
    sd_logFC_sig = round(stats::sd(sig$logFC, na.rm = TRUE), 4L)
  )
}

methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
tasks <- c("breast_subtype", "breast_vs_lung")
domains <- c("cptac", "ccle")

all_summaries <- data.table()
for (method in methods) {
  for (task in tasks) {
    for (domain in domains) {
      da_path <- file.path(REPO, sprintf(
        "reports/benchmark_master/benchmark_results/%s/%s/representation_da/%s/da_limma_result.csv",
        method, task, domain
      ))
      if (file.exists(da_path)) {
        all_summaries <- rbind(all_summaries, compute_fc_se_summary(da_path, method, task, domain))
      }
    }
  }
}

cat("\n=== FOLD CHANGE AND SE SUMMARY ===\n")
print(all_summaries)

fwrite(all_summaries, file.path(PRES_OUT, "tables/fc_se_summary.csv"))

cat("\n=== FC/SE COMPARISON ACROSS METHODS (subtype, CPTAC) ===\n")
sub_cptac <- all_summaries[task == "breast_subtype" & domain == "cptac"]
if (nrow(sub_cptac) > 0L) {
  print(sub_cptac[, .(method, n_significant, median_abs_logFC_sig, median_SE_sig, iqr_logFC_sig)])
}
