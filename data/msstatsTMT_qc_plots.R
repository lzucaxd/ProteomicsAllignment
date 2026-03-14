#!/usr/bin/env Rscript
# =============================================================================
# MSstatsTMT QC – publication-quality plots for large TMT datasets
# =============================================================================
#
# Reads msstats_input.tsv and protein_summary.tsv from --outdir and generates:
#   - QC run boxplot (log2 intensity per run, RunIndex on x-axis)
#   - Density plot of log2 intensities per run
#   - PCA of protein abundances
#   - Missing value heatmap (protein × sample)
#   - Profile plots for selected proteins (top 15 peptides/protein, no legend)
#
# Designed for 150+ runs and thousands of peptides: plots split by Mixture,
# at most 15 runs per plot, numeric RunIndex, high-resolution PDFs (40×10 in).
#
# Usage:
#   Rscript --no-init-file msstatsTMT_qc_plots.R --outdir results/PDC000120
#   Rscript --no-init-file msstatsTMT_qc_plots.R --outdir results/PDC000120 --n_proteins 10
#
# Arguments:
#   --outdir    Directory containing msstats_input.tsv and protein_summary.tsv.
#   --n_proteins  Number of proteins for profile plots. Default: 5.
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("tidyr", quietly = TRUE))
    install.packages("tidyr", repos = "https://cloud.r-project.org")
})
library(data.table)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
outdir <- "results/PDC000120"
n_proteins <- 5L
i <- 1
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--n_proteins" && i < length(args)) {
    n_proteins <- as.integer(args[i + 1])
    i <- i + 2
  } else {
    i <- i + 1
  }
}

peptide_path <- file.path(outdir, "msstats_input.tsv")
protein_path <- file.path(outdir, "protein_summary.tsv")
address <- file.path(outdir, "plots")
dir.create(address, recursive = TRUE, showWarnings = FALSE)

WIDTH <- 40
HEIGHT <- 10
DPI <- 300
RUNS_PER_PLOT <- 15L
TOP_PEPTIDES_PER_PROTEIN <- 15L

if (!file.exists(peptide_path)) stop("Not found: ", peptide_path)
if (!file.exists(protein_path)) stop("Not found: ", protein_path)

message("Loading peptide-level data: ", peptide_path)
pep <- fread(peptide_path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
setDT(pep)

message("Loading protein-level data: ", protein_path)
prot <- fread(protein_path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
setDT(prot)

# Ensure Protein column
if (!"Protein" %in% names(prot) || all(is.na(prot$Protein)) || all(prot$Protein == ""))
  prot[, Protein := ProteinName]

# ------------------------------------------------------------------------------
# Add log2 intensity and RunIndex (numeric) by Mixture
# ------------------------------------------------------------------------------
pep[, log2Intensity := log2(pmax(Intensity, 0.5))]
# Run index within each Mixture (1, 2, 3, ...)
setorder(pep, Mixture, Run)
pep[, RunIndex := match(Run, unique(Run)), by = Mixture]
# Chunk index for splitting (each chunk = up to 15 runs)
pep[, RunChunk := ceiling(RunIndex / RUNS_PER_PLOT), by = Mixture]

# Same RunIndex/RunChunk for protein-level
setorder(prot, Mixture, Run)
prot[, RunIndex := match(Run, unique(Run)), by = Mixture]
prot[, RunChunk := ceiling(RunIndex / RUNS_PER_PLOT), by = Mixture]

# Unique runs per mixture for mapping
run_order <- unique(pep[, .(Mixture, Run)])
run_order[, RunIndex := match(Run, unique(Run)), by = Mixture]
run_order[, RunChunk := ceiling(RunIndex / RUNS_PER_PLOT), by = Mixture]

message("Runs: ", uniqueN(pep$Run), " | Mixtures: ", uniqueN(pep$Mixture), " | Proteins: ", uniqueN(prot$Protein))

# ------------------------------------------------------------------------------
# 1. QC run boxplot (log2 intensity distribution per run)
# Split by Mixture and RunChunk; x = RunIndex (1..15 per plot)
# ------------------------------------------------------------------------------
message("Generating QC run boxplots...")
pep_plot <- pep[!is.na(log2Intensity) & is.finite(log2Intensity)]
# Subsample if huge (for speed)
if (nrow(pep_plot) > 5e5) pep_plot <- pep_plot[sample(.N, 5e5)]

chunks <- unique(pep_plot[, .(Mixture, RunChunk)])
for (j in seq_len(nrow(chunks))) {
  m <- chunks$Mixture[j]
  c <- chunks$RunChunk[j]
  sub <- pep_plot[Mixture == m & RunChunk == c]
  sub[, RunIndexLocal := RunIndex - (c - 1L) * RUNS_PER_PLOT]
  p <- ggplot(sub, aes(x = factor(RunIndexLocal), y = log2Intensity)) +
    geom_boxplot(outlier.alpha = 0.15, outlier.size = 0.5, fill = "steelblue", alpha = 0.8) +
    labs(
      title = sprintf("QC run boxplot: %s (runs %d–%d)", m, (c - 1L) * RUNS_PER_PLOT + 1L, min((c) * RUNS_PER_PLOT, max(sub$RunIndex))),
      x = "Run index", y = expression(log[2] ~ intensity)
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), panel.grid.minor = element_blank())
  f <- file.path(address, sprintf("QC_run_boxplot_%s_chunk%d.pdf", gsub("[^A-Za-z0-9._-]", "_", m), c))
  ggsave(f, p, width = WIDTH, height = HEIGHT, dpi = DPI)
  message("  ", f)
}

# ------------------------------------------------------------------------------
# 2. Density plot of log2 intensities per run (one curve per run, by chunk)
# ------------------------------------------------------------------------------
message("Generating density plots...")
for (j in seq_len(nrow(chunks))) {
  m <- chunks$Mixture[j]
  c <- chunks$RunChunk[j]
  sub <- pep_plot[Mixture == m & RunChunk == c]
  sub[, RunIndexLocal := RunIndex - (c - 1L) * RUNS_PER_PLOT]
  p <- ggplot(sub, aes(x = log2Intensity, color = factor(RunIndexLocal), group = RunIndexLocal)) +
    geom_density(linewidth = 0.6) +
    labs(
      title = sprintf("Log2 intensity density by run: %s (chunk %d)", m, c),
      x = expression(log[2] ~ intensity), y = "Density", color = "Run index"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
  f <- file.path(address, sprintf("QC_density_%s_chunk%d.pdf", gsub("[^A-Za-z0-9._-]", "_", m), c))
  ggsave(f, p, width = WIDTH, height = HEIGHT, dpi = DPI)
  message("  ", f)
}

# ------------------------------------------------------------------------------
# 3. PCA of protein abundances (sample = Run × Channel or BioReplicate)
# ------------------------------------------------------------------------------
message("Generating PCA plot...")
prot_wide <- dcast(prot, Protein ~ Run + Channel, value.var = "Abundance", fun.aggregate = median, na.rm = TRUE)
mat <- as.matrix(prot_wide[, -1, with = FALSE])
rownames(mat) <- prot_wide$Protein
mat[!is.finite(mat)] <- NA
# Impute NA/Inf with column median for PCA
mat_imp <- mat
for (k in seq_len(ncol(mat))) {
  v <- mat[, k]
  med <- median(v[is.finite(v)], na.rm = TRUE)
  if (!is.finite(med)) med <- 0
  mat_imp[!is.finite(v), k] <- med
}
mat_imp[!is.finite(mat_imp)] <- 0
# Drop constant columns (samples)
var_col <- apply(mat_imp, 2, sd)
mat_imp <- mat_imp[, var_col > 0, drop = FALSE]
if (ncol(mat_imp) < 2 || nrow(mat_imp) < 2) stop("Not enough variation for PCA.")
pca <- prcomp(t(mat_imp), scale. = TRUE, center = TRUE)
samp <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
# Label by Run (first part of "Run.Channel" column names) for coloring
samp$Run <- sub("\\.[^.]+$", "", samp$Sample)
run2mix <- run_order[, .(Mixture = Mixture[1]), by = Run]
samp <- as.data.frame(merge(samp, run2mix, by = "Run", all.x = TRUE))
samp$Mixture[is.na(samp$Mixture)] <- "Other"
p <- ggplot(samp, aes(PC1, PC2, color = Mixture)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "PCA of protein abundances (samples = Run × Channel)",
    x = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
    y = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2])
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", panel.grid.minor = element_blank())
ggsave(file.path(address, "PCA_protein_abundances.pdf"), p, width = WIDTH, height = HEIGHT, dpi = DPI)
message("  ", file.path(address, "PCA_protein_abundances.pdf"))

# ------------------------------------------------------------------------------
# 4. Missing value heatmap (protein × sample)
# Subsample proteins/columns if very large for readability
# ------------------------------------------------------------------------------
message("Generating missing value heatmap...")
miss <- is.na(mat) | !is.finite(mat)
# Limit to 500 proteins and 200 samples if huge for readable heatmap
np <- nrow(miss)
ns <- ncol(miss)
if (np > 500) {
  set.seed(1)
  pr <- sample(np, 500)
  miss <- miss[pr, , drop = FALSE]
}
if (ns > 200) {
  set.seed(1)
  sc <- sample(ns, 200)
  miss <- miss[, sc, drop = FALSE]
}
miss_df <- expand.grid(Protein = seq_len(nrow(miss)), Sample = seq_len(ncol(miss)))
miss_df$Missing <- as.vector(miss)
miss_df$Missing <- factor(ifelse(miss_df$Missing, "Missing", "Present"))
p <- ggplot(miss_df, aes(x = Sample, y = Protein, fill = Missing)) +
  geom_raster() +
  scale_fill_manual(values = c("Present" = "grey90", "Missing" = "darkred")) +
  labs(title = "Missing value pattern (protein × sample)", x = "Sample index", y = "Protein index") +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "top")
ggsave(file.path(address, "Missing_value_heatmap.pdf"), p, width = WIDTH, height = HEIGHT, dpi = DPI)
message("  ", file.path(address, "Missing_value_heatmap.pdf"))

# ------------------------------------------------------------------------------
# 5. Profile plots: selected proteins, top 15 peptides by median intensity,
#    no peptide legend, RunIndex on x, split by Mixture (max 15 runs per plot)
# ------------------------------------------------------------------------------
prots_in_pep <- unique(pep$ProteinName)
prots_sel <- head(unique(prot$Protein[!is.na(prot$Protein) & nzchar(prot$Protein) & prot$Protein %in% prots_in_pep]), n_proteins)
# Top 15 peptides per protein by median intensity
med_pep <- pep[, .(med = median(Intensity, na.rm = TRUE)), by = .(ProteinName, PSM)]
med_pep[, rank := frank(-med, ties.method = "first"), by = ProteinName]
top_psm <- med_pep[rank <= TOP_PEPTIDES_PER_PROTEIN, .(ProteinName, PSM)]
pep_top <- merge(pep, top_psm, by = c("ProteinName", "PSM"))

for (pr in prots_sel) {
  sub <- pep_top[ProteinName == pr]
  if (nrow(sub) == 0) next
  sub[, RunIndexLocal := RunIndex - (RunChunk - 1L) * RUNS_PER_PLOT]
  chunks_pr <- unique(sub[, .(Mixture, RunChunk)])
  for (j in seq_len(nrow(chunks_pr))) {
    m <- chunks_pr$Mixture[j]
    c <- chunks_pr$RunChunk[j]
    subc <- sub[Mixture == m & RunChunk == c]
    p <- ggplot(subc, aes(x = RunIndexLocal, y = log2Intensity, group = PSM)) +
      geom_line(linewidth = 0.4, alpha = 0.7) +
      labs(
        title = sprintf("Profile: %s | %s (runs %d–%d)", pr, m, (c - 1L) * RUNS_PER_PLOT + 1L, min(c * RUNS_PER_PLOT, max(subc$RunIndex))),
        x = "Run index", y = expression(log[2] ~ intensity)
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", panel.grid.minor = element_blank())
    safe <- gsub("[^A-Za-z0-9._-]", "_", pr)
    f <- file.path(address, sprintf("Profile_%s_%s_chunk%d.pdf", safe, gsub("[^A-Za-z0-9._-]", "_", m), c))
    ggsave(f, p, width = WIDTH, height = HEIGHT, dpi = DPI)
    message("  ", f)
  }
}

message("Done. All plots saved under: ", address)
