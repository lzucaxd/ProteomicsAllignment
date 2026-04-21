#!/usr/bin/env Rscript
# =============================================================================
# Tumor vs NAT differential abundance — CPTAC breast proteomics (PDC000120)
# MSstats + limma, then compare and volcano plots.
# =============================================================================
# Usage: from data/ directory:
#   Rscript scripts/DA_tumor_vs_NAT_CPTAC_breast.R
# Outputs written to results/PDC000120/ (DA_*.csv, volcano_*.pdf).
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("MSstats", quietly = TRUE))
    BiocManager::install("MSstats", update = FALSE, ask = FALSE)
  if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", update = FALSE, ask = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
})
library(MSstats)
library(limma)
library(data.table)

# Run from data/ directory (e.g. Rscript scripts/DA_tumor_vs_NAT_CPTAC_breast.R)
if (length(sys.frames()) > 0 && exists("ofile", sys.frame(1))) {
  DATA_DIR <- dirname(dirname(sys.frame(1)$ofile))
} else {
  DATA_DIR <- getwd()
  if (!file.exists(file.path(DATA_DIR, "results", "PDC000120", "gene_matrix.csv")))
    DATA_DIR <- file.path(getwd(), "data")
}
setwd(DATA_DIR)
RESULTS_DIR <- file.path(DATA_DIR, "results", "PDC000120")
BIOSPEC_PATH <- file.path(DATA_DIR, "biospecimen", "PDC_study_biospecimen_03162026_190026.csv")
GENE_MATRIX_PATH <- file.path(RESULTS_DIR, "gene_matrix.csv")
ANNOT_PATH <- file.path(RESULTS_DIR, "annotation_filled_corrected.csv")
OUT_DIR <- RESULTS_DIR
dir.create(OUT_DIR, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Step 1 — Data preparation
# -----------------------------------------------------------------------------
message("Step 1: Data preparation")
gm <- fread(GENE_MATRIX_PATH)
rownames(gm) <- gm$GeneSymbol
gm$UniProtID <- NULL
aliquot_ids <- setdiff(names(gm), "GeneSymbol")
# Aliquot IDs are the 144 sample columns
mat <- as.matrix(gm[, ..aliquot_ids])
rownames(mat) <- gm$GeneSymbol

bio <- fread(BIOSPEC_PATH)
# Trim BOM/whitespace from header
setnames(bio, trimws(gsub("^\uFEFF", "", names(bio))))
aliquot_col <- names(bio)[grepl("Aliquot.*Submitter|Submitter.*ID", names(bio), ignore.case = TRUE)][1]
type_col <- names(bio)[grepl("Sample.Type", names(bio), ignore.case = TRUE)][1]
if (is.na(aliquot_col)) aliquot_col <- "Aliquot Submitter ID"
if (is.na(type_col)) type_col <- "Sample Type"
aliquot_to_type <- setNames(trimws(bio[[type_col]]), trimws(bio[[aliquot_col]]))

# Remove RetroIR column (Sample Type = Not Reported or column name is "RetroIR")
drop <- "RetroIR"
if (drop %in% colnames(mat)) {
  mat <- mat[, setdiff(colnames(mat), drop), drop = FALSE]
  aliquot_ids <- colnames(mat)
}
message("  Working matrix: ", nrow(mat), " x ", ncol(mat))

# Assign group: Tumor (Primary Tumor) vs NAT (Solid Tissue Normal)
aliquot_to_group <- setNames(
  ifelse(aliquot_to_type[aliquot_ids] == "Primary Tumor", "Tumor",
         ifelse(aliquot_to_type[aliquot_ids] == "Solid Tissue Normal", "NAT", NA_character_)),
  aliquot_ids
)
# Drop any aliquot not in biospecimen or not Tumor/NAT
keep <- !is.na(aliquot_to_group)
mat <- mat[, keep, drop = FALSE]
aliquot_ids <- colnames(mat)
aliquot_to_group <- aliquot_to_group[aliquot_ids]
message("  Tumor: ", sum(aliquot_to_group == "Tumor"), ", NAT: ", sum(aliquot_to_group == "NAT"))
group <- factor(aliquot_to_group, levels = c("NAT", "Tumor"))

# -----------------------------------------------------------------------------
# Step 2 — Build MSstats long format
# -----------------------------------------------------------------------------
message("Step 2: Build MSstats input (long format)")
ann <- fread(ANNOT_PATH)
# Run for MSstats long format: use aliquot as Run to avoid internal join explosion.
aliquot_to_run <- setNames(aliquot_ids, aliquot_ids)
# Mixture (TMT plex) per aliquot for limma batch covariate
ann_sample <- ann[tolower(Condition) != "norm", .(BioReplicate, Mixture)]
ann_sample <- unique(ann_sample[, .(Mixture = Mixture[1L]), by = BioReplicate])
aliquot_to_mixture <- setNames(ann_sample$Mixture, ann_sample$BioReplicate)

# One row per (Protein, aliquot). Matrix is log2-intensity; MSstats will log2 again, so pass 2^intensity
proteins <- rownames(mat)
long_list <- lapply(seq_along(proteins), function(i) {
  pr <- proteins[i]
  intens <- mat[i, ]
  log2val <- as.numeric(intens)
  data.table(
    ProteinName = pr,
    PeptideSequence = "PROTEIN_SUM",
    PrecursorCharge = 0L,
    FragmentIon = "NA",
    ProductCharge = 0L,
    IsotopeLabelType = "L",
    Condition = aliquot_to_group[names(intens)],
    BioReplicate = names(intens),
    Run = aliquot_to_run[names(intens)],
    Intensity = 2^log2val
  )
})
msstats_long <- rbindlist(long_list)
# Remove NA Run (aliquot not in annotation)
msstats_long <- msstats_long[!is.na(Run)]
# Pass as data.frame to avoid data.table merge explosion inside MSstats
msstats_long <- as.data.frame(msstats_long)
message("  MSstats long format: ", nrow(msstats_long), " rows")

# -----------------------------------------------------------------------------
# Step 3 — MSstats groupComparison (fallback: limma ~ 0 + group if MSstats fails)
# -----------------------------------------------------------------------------
message("Step 3: MSstats dataProcess + groupComparison")
msstats_ok <- FALSE
tryCatch({
  processed <- dataProcess(
    msstats_long,
    normalization = FALSE,
    summaryMethod = "TMP",
    MBimpute = FALSE
  )
  contrast_mat <- matrix(0, nrow = 1, ncol = 2)
  rownames(contrast_mat) <- "Tumor-NAT"
  colnames(contrast_mat) <- c("Tumor", "NAT")
  contrast_mat[1, "Tumor"] <- 1
  contrast_mat[1, "NAT"] <- -1
  comparison <- groupComparison(contrast.matrix = contrast_mat, data = processed)
  msstats_res <- as.data.frame(comparison$ComparisonResult)
  msstats_res$ProteinName <- as.character(msstats_res$ProteinName)
  if (!"log2FC" %in% names(msstats_res) && "logFC" %in% names(msstats_res)) msstats_res$log2FC <- msstats_res$logFC
  if (!"adj.pvalue" %in% names(msstats_res) && "Adjusted.Pvalue" %in% names(msstats_res)) msstats_res$adj.pvalue <- msstats_res$Adjusted.Pvalue
  msstats_ok <- TRUE
  message("  MSstats comparison: ", nrow(msstats_res), " proteins")
}, error = function(e) {
  message("  MSstats failed (", conditionMessage(e), "). Using limma ~ 0 + group as stand-in.")
})

if (!msstats_ok) {
  design_simple <- model.matrix(~ 0 + group)
  colnames(design_simple) <- gsub("group", "", colnames(design_simple))
  fit_simple <- lmFit(mat, design_simple)
  contrast_simple <- makeContrasts(Tumor - NAT, levels = design_simple)
  fit_simple2 <- contrasts.fit(fit_simple, contrast_simple)
  fit_simple2 <- eBayes(fit_simple2)
  msstats_res <- topTable(fit_simple2, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
  msstats_res$ProteinName <- rownames(msstats_res)
  msstats_res$log2FC <- msstats_res$logFC
  msstats_res$adj.pvalue <- msstats_res$adj.P.Val
  msstats_res$pvalue <- msstats_res$P.Value
  message("  MSstats stand-in (limma no batch): ", nrow(msstats_res), " proteins")
}

# -----------------------------------------------------------------------------
# Step 4 — limma (same 143-column matrix, with/without Mixture)
# -----------------------------------------------------------------------------
message("Step 4: limma")
group <- factor(aliquot_to_group, levels = c("NAT", "Tumor"))
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
fit <- lmFit(mat, design)
contrast <- makeContrasts(Tumor - NAT, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
limma_res <- topTable(fit2, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
limma_res$GeneSymbol <- rownames(limma_res)
limma_res <- limma_res[, c("GeneSymbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# With Mixture (batch) covariate
mixture_per_aliquot <- aliquot_to_mixture[colnames(mat)]
mixture_per_aliquot[is.na(mixture_per_aliquot)] <- "Unknown"
mixture_f <- factor(mixture_per_aliquot)
design_mix <- model.matrix(~ 0 + group + mixture_f)
fit_mix <- lmFit(mat, design_mix)
# Contrast Tumor - NAT via makeContrasts (handles coefficient names)
cm <- makeContrasts("groupTumor - groupNAT", levels = design_mix)
fit_mix2 <- tryCatch(contrasts.fit(fit_mix, cm), error = function(e) NULL)
if (!is.null(fit_mix2)) {
  fit_mix2 <- eBayes(fit_mix2)
  limma_res <- topTable(fit_mix2, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
}
limma_res$GeneSymbol <- rownames(limma_res)
limma_res <- limma_res[, c("GeneSymbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# -----------------------------------------------------------------------------
# Step 5 — Compare MSstats vs limma
# -----------------------------------------------------------------------------
message("Step 5: Merge and compare MSstats vs limma")
msstats_res$GeneSymbol <- msstats_res$ProteinName
merged <- merge(
  msstats_res[, c("GeneSymbol", "log2FC", "adj.pvalue", "pvalue")],
  limma_res[, c("GeneSymbol", "logFC", "adj.P.Val", "P.Value")],
  by = "GeneSymbol",
  all = TRUE,
  suffixes = c("_MSstats", "_limma")
)
names(merged)[names(merged) == "log2FC"] <- "log2FC_MSstats"
names(merged)[names(merged) == "logFC"] <- "log2FC_limma"

# Correlation of log2FC
ok <- complete.cases(merged$log2FC_MSstats, merged$log2FC_limma)
if (sum(ok) > 2) {
  cor_fc <- cor(merged$log2FC_MSstats[ok], merged$log2FC_limma[ok], use = "pairwise.complete.obs")
  message("  log2FC correlation (MSstats vs limma): ", round(cor_fc, 4))
}

# Overlap of significant hits (FDR < 0.05, |log2FC| > 1)
sig_msstats <- merged$GeneSymbol[which(merged$adj.pvalue < 0.05 & abs(merged$log2FC_MSstats) > 1)]
sig_limma   <- merged$GeneSymbol[which(merged$adj.P.Val < 0.05 & abs(merged$log2FC_limma) > 1)]
overlap <- intersect(sig_msstats, sig_limma)
message("  Significant (FDR<0.05, |log2FC|>1): MSstats ", length(sig_msstats),
        ", limma ", length(sig_limma), ", overlap ", length(overlap))

# -----------------------------------------------------------------------------
# Step 6 — Volcano plots (MSstats and limma)
# -----------------------------------------------------------------------------
message("Step 6: Volcano plots")
MARKERS <- c("ESR1", "GATA3", "KRT18", "PTPRC", "COL1A1", "VIM")

volcano_plot <- function(res, log2FC_col, pval_col, title, out_path, gene_col = "GeneSymbol") {
  res <- as.data.frame(res)
  x <- res[[log2FC_col]]
  y <- -log10(res[[pval_col]])
  y[is.infinite(y)] <- max(y[is.finite(y)], na.rm = TRUE) + 0.5
  sig <- !is.na(res[[pval_col]]) & res[[pval_col]] < 0.05 & abs(x) > 1
  genes <- res[[gene_col]]
  top20 <- order(res[[pval_col]])[seq_len(min(20, sum(sig, na.rm = TRUE)))]
  label_genes <- genes[top20]
  label_genes <- label_genes[!is.na(label_genes) & nzchar(label_genes)]
  pdf(out_path, width = 8, height = 7)
  plot(x, y, pch = 20, col = ifelse(sig, "red", "gray50"),
       xlab = "log2FC", ylab = "-log10(adj.pvalue)", main = title, cex = 0.6)
  abline(h = -log10(0.05), lty = 2, col = "gray40")
  abline(v = c(-1, 1), lty = 2, col = "gray40")
  if (length(label_genes) > 0) {
    idx <- match(label_genes, genes)
    idx <- idx[!is.na(idx)]
    text(x[idx], y[idx], labels = genes[idx], pos = 4, cex = 0.5)
  }
  # Mark known markers
  for (m in MARKERS) {
    j <- match(m, genes)
    if (!is.na(j)) points(x[j], y[j], pch = 8, col = "blue", cex = 1.2)
  }
  legend("topright", legend = c("FDR<0.05, |log2FC|>1", "Known markers"), col = c("red", "blue"), pch = c(20, 8), bty = "n")
  dev.off()
}

volcano_plot(msstats_res, "log2FC", "adj.pvalue", "Tumor vs NAT (MSstats)", file.path(OUT_DIR, "volcano_MSstats.pdf"))
volcano_plot(limma_res, "logFC", "adj.P.Val", "Tumor vs NAT (limma + Mixture)", file.path(OUT_DIR, "volcano_limma.pdf"), "GeneSymbol")

# -----------------------------------------------------------------------------
# Step 7 — Save outputs
# -----------------------------------------------------------------------------
message("Step 7: Save outputs")
fwrite(msstats_res, file.path(OUT_DIR, "DA_MSstats_tumor_vs_NAT.csv"))
fwrite(limma_res, file.path(OUT_DIR, "DA_limma_tumor_vs_NAT.csv"))
fwrite(merged, file.path(OUT_DIR, "DA_merged_comparison.csv"))

# -----------------------------------------------------------------------------
# Final summary
# -----------------------------------------------------------------------------
message("\n========== Final summary ==========")
message("Total significant (FDR < 0.05, |log2FC| > 1):")
message("  MSstats: ", length(sig_msstats))
message("  limma:   ", length(sig_limma))
message("  Overlap: ", length(overlap))

msstats_res_dt <- as.data.table(msstats_res)
msstats_res_dt <- msstats_res_dt[order(-log2FC)]
message("\nTop 10 upregulated in Tumor (MSstats):")
print(msstats_res_dt[log2FC > 0, .(GeneSymbol, log2FC, adj.pvalue)][1:10])
message("\nTop 10 upregulated in NAT (MSstats):")
print(msstats_res_dt[log2FC < 0][order(log2FC)][1:10, .(GeneSymbol, log2FC, adj.pvalue)])

message("\nSanity check markers (expected: Tumor up = ESR1, GATA3, KRT18; NAT up = PTPRC, COL1A1, VIM):")
for (m in MARKERS) {
  r <- msstats_res[msstats_res$ProteinName == m, ]
  if (nrow(r) > 0) {
    dir <- if (r$log2FC > 0) "Tumor" else "NAT"
    sig <- if (r$adj.pvalue < 0.05) "significant" else "not significant"
    message("  ", m, ": log2FC = ", round(r$log2FC, 3), " (", dir, "), ", sig, " (adj.p = ", format(r$adj.pvalue, digits = 3), ")")
  } else {
    message("  ", m, ": not in results")
  }
}
message("\nDone. Outputs in ", OUT_DIR)
