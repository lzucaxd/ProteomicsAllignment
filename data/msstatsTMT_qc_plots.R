#!/usr/bin/env Rscript
# =============================================================================
# MSstatsTMT built-in QC plots only (for lab debugging)
# =============================================================================
#
# Calls MSstatsTMT::dataProcessPlotsTMT to generate only the package’s native
# QCPlot (box plots of log intensities across channels and MS runs) and
# optional ProfilePlots. No custom plots.
#
# Usage:
#   Rscript --no-init-file msstatsTMT_qc_plots.R --outdir results/PDC000120
#   Rscript --no-init-file msstatsTMT_qc_plots.R --outdir results/PDC000120 --n_profile_proteins 5
#
# Arguments:
#   --outdir             Directory containing msstats_input.tsv and protein_summary.tsv.
#   --n_profile_proteins Number of proteins for ProfilePlot (0 = QCPlot only). Default: 0.
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

# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
outdir <- "."
n_profile_proteins <- 0L
i <- 1
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--n_profile_proteins" && i < length(args)) {
    n_profile_proteins <- as.integer(args[i + 1])
    i <- i + 2
  } else {
    i <- i + 1
  }
}

peptide_path <- file.path(outdir, "msstats_input.tsv")
protein_path <- file.path(outdir, "protein_summary.tsv")
plots_dir <- file.path(outdir, "plots")

if (!file.exists(peptide_path)) stop("Not found: ", peptide_path)
if (!file.exists(protein_path)) stop("Not found: ", protein_path)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading peptide-level data: ", peptide_path)
data_peptide <- fread(peptide_path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
setDT(data_peptide)

message("Loading protein-level data: ", protein_path)
data_summarization <- fread(protein_path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
setDT(data_summarization)

# MSstatsTMT may expect ProteinName in summarization
if (!"Protein" %in% names(data_summarization) && "ProteinName" %in% names(data_summarization))
  data_summarization[, Protein := ProteinName]
if (!"ProteinName" %in% names(data_summarization) && "Protein" %in% names(data_summarization))
  data_summarization[, ProteinName := Protein]

# MSstatsTMT 2.x expects data = output of proteinSummarization (list with FeatureLevelData, ProteinLevelData).
# If you have only saved peptide + protein tables, use scripts/ccle_regenerate_qc_plots.R to re-run summarization and plot.
message("Generating MSstatsTMT QC plot (readable size)...")
tryCatch({
  dataProcessPlotsTMT(
    data = list(FeatureLevelData = data_peptide, ProteinLevelData = data_summarization),
    type = "QCPlot",
    which.Protein = "allonly",
    width = 20,
    height = 12,
    x.axis.size = 12,
    y.axis.size = 12,
    text.size = 4,
    legend.size = 9,
    address = plots_dir
  )
  message("  QCPlot saved to ", plots_dir)
}, error = function(e1) {
  tryCatch({
    dataProcessPlotsTMT(
      data.peptide = data_peptide,
      data.summarization = data_summarization,
      type = "QCPlot",
      which.Protein = "allonly",
      width = 20,
      height = 12,
      address = plots_dir
    )
    message("  QCPlot saved to ", plots_dir)
  }, error = function(e2) stop("QCPlot failed: ", conditionMessage(e1)))
})

if (n_profile_proteins > 0L) {
  prot_col <- if ("Protein" %in% names(data_summarization)) "Protein" else "ProteinName"
  prots <- unique(data_summarization[[prot_col]])
  prots <- head(prots[!is.na(prots) & nzchar(prots)], n_profile_proteins)
  if (length(prots) > 0L) {
    message("Generating MSstatsTMT ProfilePlots for ", length(prots), " proteins...")
    tryCatch({
      dataProcessPlotsTMT(
        data.peptide = data_peptide,
        data.summarization = data_summarization,
        type = "ProfilePlot",
        which.Protein = prots,
        width = 14,
        height = 7,
        address = plots_dir
      )
      message("  ProfilePlots saved to ", plots_dir)
    }, error = function(e) warning("ProfilePlot failed: ", conditionMessage(e)))
  }
}

message("Done. MSstatsTMT plots in: ", plots_dir)
