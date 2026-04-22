#!/usr/bin/env Rscript
# QC plots for MSstatsTMT bridge / Norm channels.
#
# IMPORTANT: protein_summary.tsv typically *drops* the reference TMT channel(s) after normalization,
# so there is no channel 131 / Norm in that file. Bridge abundances are read from msstats_input.tsv
# (PSM-level intensities, Condition == Norm) — the same input used for summarization.
#
# CPTAC (PDC000120): Norm = pooled reference on channel 131 (BioReplicate POOL).
# CCLE: mixture 0 uses multiple reference cell lines (all Condition == Norm in that plex).
#
# Usage (from repo root):
#   Rscript data/scripts/qc_bridge_norm_protein.R PDC000120
#   Rscript data/scripts/qc_bridge_norm_protein.R CCLE
#
# Outputs: data/results/{study}/qc_bridge/

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript qc_bridge_norm_protein.R <PDC000120|CCLE>")
}

study <- args[[1]]
ca <- commandArgs(trailingOnly = FALSE)
fa <- ca[startsWith(ca, "--file=")]
if (length(fa)) {
  script_path <- sub("^--file=", "", fa[[1]])
  root <- dirname(dirname(dirname(normalizePath(script_path))))
} else {
  root <- getwd()
}

res_dir <- file.path(root, "data", "results", study)
msstats_path <- file.path(res_dir, "msstats_input.tsv")
ann_path <- file.path(res_dir, "annotation_filled.csv")
out_dir <- file.path(res_dir, "qc_bridge")

if (!file.exists(msstats_path)) stop("Missing: ", msstats_path)
if (!file.exists(ann_path)) stop("Missing: ", ann_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Stream-filter Norm rows (column 8 = Condition in MSstatsTMT input)
awk_cmd <- sprintf(
  "awk -F'\\t' 'NR==1 || $8==\"Norm\"' %s",
  shQuote(msstats_path, type = "sh")
)

message("Loading Norm rows from msstats_input (via awk) ...")
bn <- fread(
  cmd = awk_cmd,
  sep = "\t",
  quote = "",
  showProgress = TRUE,
  select = c("Mixture", "TechRepMixture", "Run", "Channel", "BioReplicate", "Condition", "Intensity")
)
if (nrow(bn) == 0) stop("No Norm rows found in msstats_input.")

bn[, Channel := as.character(Channel)]
bn[, log2_intensity := log2(pmax(as.numeric(Intensity), 0, na.rm = TRUE) + 1)]

# Per-run medians (all Norm PSMs pooled per run)
run_one <- bn[, .(
  median_log2_intensity = as.numeric(median(log2_intensity, na.rm = TRUE)),
  n_psm = .N
), by = .(Run, Mixture)]

# Per-run, per-channel medians
run_med <- bn[, .(
  median_log2_intensity = as.numeric(median(log2_intensity, na.rm = TRUE)),
  n_psm = .N
), by = .(Run, Mixture, Channel, BioReplicate)]

fwrite(run_one, file.path(out_dir, "bridge_norm_run_median_pooled.tsv"), sep = "\t")
fwrite(run_med, file.path(out_dir, "bridge_norm_run_medians_by_channel.tsv"), sep = "\t")

# Annotation summary
ann <- fread(ann_path)
if (!"Condition" %in% names(ann)) stop("annotation_filled.csv missing Condition column")
ann_n <- ann[Condition == "Norm"]
fwrite(
  ann_n[, .N, by = .(Mixture, Channel, BioReplicate)][order(Mixture, Channel)],
  file.path(out_dir, "bridge_annotation_norm_counts.tsv"),
  sep = "\t"
)

# --- Plots (subsample for histogram if huge) ---
theme_set(theme_bw(base_size = 11))
set.seed(1)
bn_plot <- if (nrow(bn) > 500000L) bn[sample.int(.N, 500000L)] else copy(bn)

p1 <- ggplot(bn_plot, aes(log2_intensity)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "grey35", colour = NA, alpha = 0.85) +
  geom_density(colour = "steelblue", linewidth = 0.9) +
  labs(
    title = sprintf("%s — Norm/bridge PSM intensities (subsampled for histogram if n > 5e5)", study),
    subtitle = "log2(intensity + 1) from msstats_input.tsv; reference channel excluded from protein_summary by design",
    x = "log2(intensity + 1)",
    y = "Density"
  )
ggsave(file.path(out_dir, "qc_bridge_overall_histogram_density.png"), p1, width = 8, height = 4.5, dpi = 150)

p2 <- ggplot(run_one, aes(x = reorder(Run, median_log2_intensity), y = median_log2_intensity)) +
  geom_point(size = 0.65, alpha = 0.65, colour = "grey25") +
  geom_hline(
    yintercept = median(run_one$median_log2_intensity, na.rm = TRUE),
    linetype = 2,
    colour = "red3",
    alpha = 0.75
  ) +
  labs(
    title = sprintf("%s — median Norm abundance per run (pooled channels)", study),
    subtitle = "Each point is one run; y = median log2 PSM intensity across all Norm rows in that run",
    x = "Run (ordered by median)",
    y = "Median log2(intensity + 1)"
  ) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(out_dir, "qc_bridge_run_median_scatter.png"), p2, width = 9, height = 4.5, dpi = 150)

p2b <- ggplot(run_one, aes(x = study, y = median_log2_intensity)) +
  geom_boxplot(width = 0.35, fill = "lightsteelblue", alpha = 0.85, outlier.size = 0.8) +
  labs(
    title = sprintf("%s — spread of run-level median Norm intensities", study),
    x = NULL,
    y = "Median log2(intensity + 1) per run"
  )
ggsave(file.path(out_dir, "qc_bridge_run_median_boxplot.png"), p2b, width = 4, height = 4.5, dpi = 150)

p3 <- ggplot(run_one, aes(x = factor(Mixture), y = median_log2_intensity)) +
  geom_boxplot(outlier.size = 0.4, fill = "wheat", alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.35) +
  labs(
    title = sprintf("%s — run medians by TMT mixture", study),
    x = "Mixture",
    y = "Median log2(intensity + 1) (Norm pooled)"
  )
ggsave(file.path(out_dir, "qc_bridge_by_mixture_boxplot.png"), p3, width = 8, height = 4.5, dpi = 150)

p4 <- ggplot(run_med, aes(x = reorder(Channel, median_log2_intensity, median), y = median_log2_intensity)) +
  geom_boxplot(outlier.size = 0.4, fill = "mistyrose2", alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.35) +
  labs(
    title = sprintf("%s — run medians by TMT channel (Norm only)", study),
    x = "Channel",
    y = "Median log2(intensity + 1)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "qc_bridge_by_channel_boxplot.png"), p4, width = 8, height = 4.8, dpi = 150)

n_runs <- uniqueN(run_one$Run)
n_mix <- uniqueN(run_one$Mixture)
n_ch <- uniqueN(run_med$Channel)
n_psm <- nrow(bn)
med_run <- median(run_one$median_log2_intensity, na.rm = TRUE)
iqr_run <- IQR(run_one$median_log2_intensity, na.rm = TRUE)
mad_run <- mad(run_one$median_log2_intensity, na.rm = TRUE)

lines <- c(
  sprintf("Bridge / Norm QC — %s", study),
  "Source: msstats_input.tsv (PSM-level reporter intensities; Condition == Norm).",
  "Note: protein_summary.tsv omits the reference MS3 channel(s) used for normalization, so bridge QC uses the same input file as MSstatsTMT summarization.",
  "",
  "Design notes:",
  if (study == "PDC000120") {
    c(
      "- CPTAC prospective breast: Norm rows are pooled reference (BioReplicate POOL) on channel 131.",
      "- MSstatsTMT reference normalization uses this channel; summarized protein abundances are relative and no longer list 131 as a separate row."
    )
  } else {
    c(
      "- CCLE: mixture 0 is the multi-line reference plex; every channel is labeled Norm (bridge lines per README).",
      "- This is not a single pooled channel like CPTAC; expect broader cross-channel spread even when stable."
    )
  },
  "",
  sprintf("Counts:"),
  sprintf("- Norm PSM rows: %s", format(n_psm, big.mark = ",")),
  sprintf("- Distinct runs (Norm present): %d", n_runs),
  sprintf("- Distinct mixtures (Norm present): %d", n_mix),
  sprintf("- Distinct channels (Norm present): %d", n_ch),
  "",
  "Stability (run-level medians of log2 PSM intensities):",
  sprintf("- Global median of run medians: %.4f", med_run),
  sprintf("- IQR of run medians: %.4f", iqr_run),
  sprintf("- MAD of run medians: %.4f", mad_run),
  "",
  "Interpretation:",
  "- Tight run-level medians suggest a stable reference anchor; heavy tails or outliers in the run plot point to specific LC-MS runs or fractions.",
  "- Compare CPTAC vs CCLE qualitatively only: designs differ (single POOL vs multi-line reference plex).",
  ""
)

writeLines(lines, file.path(out_dir, "qc_bridge_summary.txt"))

message("Done. Outputs in ", out_dir)
