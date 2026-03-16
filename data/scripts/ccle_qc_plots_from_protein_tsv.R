#!/usr/bin/env Rscript
# Generate simple QC PNGs from existing protein_summary.tsv (no re-summarization).
# Use when MSstatsTMT PDF is corrupted or missing.
#
# Usage (from data/):
#   Rscript --no-init-file scripts/ccle_qc_plots_from_protein_tsv.R --outdir results/CCLE

args <- commandArgs(trailingOnly = TRUE)
outdir <- "results/CCLE"
i <- 1
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]
    i <- i + 2
  } else i <- i + 1
}

prot_path <- file.path(outdir, "protein_summary.tsv")
plots_dir <- file.path(outdir, "plots")
if (!file.exists(prot_path)) stop("Not found: ", prot_path)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(data.table))
q <- fread(prot_path, sep = "\t", header = TRUE)
if (nrow(q) == 0L || !"Abundance" %in% names(q)) stop("No Abundance column or empty file.")

run_col <- if ("Run" %in% names(q)) "Run" else NULL
png(file.path(plots_dir, "qc_abundance_boxplot.png"), width = 1200, height = 600, res = 100)
par(mar = c(8, 4, 2, 1))
if (length(run_col) && uniqueN(q[[run_col]]) <= 100L) {
  boxplot(Abundance ~ get(run_col), data = q, las = 2, main = "Protein abundance by run", ylab = "Abundance")
} else {
  hist(q$Abundance, breaks = 80, main = "Protein abundance distribution", xlab = "Abundance", col = "steelblue")
}
dev.off()

png(file.path(plots_dir, "qc_abundance_hist.png"), width = 600, height = 500, res = 100)
hist(q$Abundance, breaks = 80, main = "Protein abundance distribution", xlab = "Abundance", col = "steelblue")
dev.off()

message("QC PNGs saved in ", plots_dir, ": qc_abundance_boxplot.png, qc_abundance_hist.png")
