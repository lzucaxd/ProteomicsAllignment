#!/usr/bin/env Rscript
# Bridge channel run medians: log2 scale from qc_bridge TSV (CPTAC).
# Raw linear intensities are not in this summary file; see diagnostics_day1_summary.md.
#
# Usage: Rscript --vanilla data/scripts/diagnostics_bridge_plots_benchmark_v1.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

repo <- getwd()
if (!file.exists(file.path(repo, "data", "results", "PDC000120", "qc_bridge", "bridge_norm_run_medians_by_channel.tsv")))
  repo <- normalizePath(file.path(repo, ".."))

path <- file.path(repo, "data", "results", "PDC000120", "qc_bridge", "bridge_norm_run_medians_by_channel.tsv")
if (!file.exists(path)) stop("Missing ", path)

out_dir <- file.path(repo, "reports", "benchmark_v1", "diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bn <- fread(path)
bn[, Channel := as.character(Channel)]
# Bridge / POOL reference
bridge <- bn[Channel == "131" | BioReplicate == "POOL" | grepl("POOL", BioReplicate)]

p1 <- ggplot(bridge, aes(median_log2_intensity)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  theme_bw() +
  labs(
    title = "CPTAC PDC000120 — bridge (Norm/POOL) run medians",
    subtitle = "log2 intensity + 1 scale from qc_bridge pipeline",
    x = "Median log2 intensity (per run)",
    y = "Count"
  )
ggsave(file.path(out_dir, "cptac_bridge_histogram.pdf"), p1, width = 7, height = 5, dpi = 200)
ggsave(file.path(out_dir, "cptac_bridge_histogram.png"), p1, width = 7, height = 5, dpi = 200)

p2 <- ggplot(bn, aes(x = factor(Channel), y = median_log2_intensity)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "CPTAC — run median log2 intensity by TMT channel",
    subtitle = "All channels including bridge (typically 131)",
    x = "Channel",
    y = "Median log2 intensity"
  )
ggsave(file.path(out_dir, "bridge_boxplots_log.pdf"), p2, width = 10, height = 5, dpi = 200)
ggsave(file.path(out_dir, "bridge_boxplots_log.png"), p2, width = 10, height = 5, dpi = 200)

message("Wrote bridge plots to ", out_dir)
