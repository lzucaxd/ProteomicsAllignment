#!/usr/bin/env Rscript
# =============================================================================
# Regenerate MSstatsTMT QC plots for CCLE with readable size and fonts.
# Loads peptide input, re-runs proteinSummarization to get the list structure
# required by dataProcessPlotsTMT(data = ...), then generates QCPlot.
#
# Usage (from data/):
#   Rscript --no-init-file scripts/ccle_regenerate_qc_plots.R --outdir results/CCLE
#
# Note: Summarization may take a long time on large peptide data (~30+ min).
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("MSstatsTMT", quietly = TRUE))
    BiocManager::install("MSstatsTMT", update = FALSE, ask = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
})
library(MSstatsTMT)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
outdir <- "results/CCLE"
i <- 1
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

peptide_path <- file.path(outdir, "msstats_input.tsv")
plots_dir <- file.path(outdir, "plots")
if (!file.exists(peptide_path)) stop("Not found: ", peptide_path)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading peptide-level data: ", peptide_path)
input_dt <- fread(peptide_path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
setDT(input_dt)
input_dt <- input_dt[!is.na(Intensity) & Intensity > 0]
has_norm <- any(tolower(input_dt$Condition) == "norm")

message("Running proteinSummarization (this may take a while)...")
quant_result <- proteinSummarization(
  input_dt,
  method = "msstats",
  global_norm = TRUE,
  reference_norm = has_norm,
  remove_norm_channel = TRUE,
  remove_empty_channel = TRUE,
  MBimpute = TRUE
)

if (!is.list(quant_result) || !"ProteinLevelData" %in% names(quant_result)) {
  stop("proteinSummarization did not return expected list structure (FeatureLevelData, ProteinLevelData).")
}

plot_prefix <- paste0(plots_dir, .Platform$file.sep)
message("Generating readable MSstatsTMT QC plot in ", plots_dir, " ...")
tryCatch({
  dataProcessPlotsTMT(
    data = quant_result,
    type = "QCPlot",
    which.Protein = "allonly",
    width = 20,
    height = 12,
    x.axis.size = 12,
    y.axis.size = 12,
    text.size = 4,
    legend.size = 9,
    address = plot_prefix
  )
  message("QCPlot saved.")
}, error = function(e) warning("MSstatsTMT QC plot failed: ", conditionMessage(e)))
# Fallback: simple QC PNGs
tryCatch({
  q <- as.data.table(quant_result$ProteinLevelData)
  if (nrow(q) > 0L && "Abundance" %in% names(q)) {
    png(file.path(plots_dir, "qc_abundance_hist.png"), width = 600, height = 500, res = 100)
    hist(q$Abundance, breaks = 80, main = "Protein abundance distribution", xlab = "Abundance", col = "steelblue")
    dev.off()
    message("Fallback QC PNG saved: qc_abundance_hist.png")
  }
}, error = function(e) warning("Fallback PNG failed: ", conditionMessage(e)))
message("Done. Plots in: ", plots_dir)
