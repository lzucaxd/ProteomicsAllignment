#!/usr/bin/env Rscript
# CCLE lineage DA — first pass (Breast vs Lung), limma on gene_matrix.
# Metadata: data/ccle_peptide/sample_info_ccle.csv (Tissue of Origin).
# Choice rationale: data/results/CCLE/ccle_lineage_contrast_choice.txt
#
# Run from repo root:
#   Rscript data/scripts/ccle_DA_lineage_firstpass.R

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Install ggrepel: install.packages(\"ggrepel\")")
  library(ggrepel)
})

root <- getwd()
gm_path <- file.path(root, "data", "results", "CCLE", "gene_matrix.csv")
ann_path <- file.path(root, "data", "ccle_peptide", "sample_info_ccle.csv")
if (!file.exists(gm_path)) {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  gm_path <- file.path(root, "data", "results", "CCLE", "gene_matrix.csv")
  ann_path <- file.path(root, "data", "ccle_peptide", "sample_info_ccle.csv")
}
if (!file.exists(gm_path)) stop("Cannot find gene_matrix.csv; run from repo root.")

out_dir <- file.path(dirname(gm_path), "DA_lineage_firstpass")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bucket <- function(tissue) {
  switch(tissue,
    "Breast" = "breast",
    "Kidney" = "kidney",
    "Lung" = "lung",
    "Pancreas" = "pancreas",
    "Stomach" = "stomach",
    "Central Nervous System" = "cns",
    NA_character_
  )
}

ann <- fread(ann_path, sep = ",", showProgress = FALSE, fill = TRUE)
# first row per cell line (tissue stable)
ann <- unique(ann, by = "Cell Line")
ann[, b := vapply(`Tissue of Origin`, function(x) bucket(as.character(x)), character(1))]

message("Reading gene matrix ...")
gm <- fread(gm_path, showProgress = FALSE)
cols_all <- names(gm)[-(1:2)]

breast <- intersect(cols_all, ann[b == "breast", `Cell Line`])
lung   <- intersect(cols_all, ann[b == "lung", `Cell Line`])
if (length(breast) < 5L || length(lung) < 5L) {
  stop("Too few columns for breast/lung: ", length(breast), " / ", length(lung))
}

samples <- c(breast, lung)
M <- as.matrix(gm[, ..samples])
storage.mode(M) <- "numeric"
rownames(M) <- gm[[2]]
gene_sym <- as.character(gm[[1]])
M[!is.finite(M)] <- NA

min_per_group <- 10L
ok <- rowSums(is.finite(M[, breast, drop = FALSE])) >= min_per_group &
  rowSums(is.finite(M[, lung, drop = FALSE])) >= min_per_group
M <- M[ok, , drop = FALSE]
gene_sym <- gene_sym[ok]
message("Genes after filter (>=", min_per_group, " observed per group): ", nrow(M))

group <- factor(
  c(rep("Breast", length(breast)), rep("Lung", length(lung))),
  levels = c("Breast", "Lung")
)
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Breast", "Lung")
fit <- lmFit(M, design)
contr <- makeContrasts(Breast_vs_Lung = Breast - Lung, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
tt <- as.data.table(topTable(fit2, coef = 1, number = Inf, sort.by = "P"), keep.rownames = "UniProtID")
tt[, GeneSymbol := gene_sym[match(UniProtID, rownames(M))]]
fwrite(tt, file.path(out_dir, "DA_Breast_vs_Lung_limma.csv"))

tt$sig <- tt$adj.P.Val < 0.1 & abs(tt$logFC) > 0.5
p_volc <- ggplot(tt, aes(logFC, -log10(P.Value))) +
  geom_point(aes(colour = sig), alpha = 0.25, size = 0.4) +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey50")) +
  theme_bw() +
  labs(
    title = "CCLE Breast vs Lung (lineage benchmark, first pass)",
    subtitle = paste0(
      "limma; n=", length(breast), " breast vs ", length(lung),
      " lung lines; positive logFC = higher in breast — exploratory, not definitive"
    ),
    x = "logFC (Breast - Lung)",
    y = "-log10 p"
  )
ggsave(file.path(out_dir, "volcano_Breast_vs_Lung.pdf"), p_volc, width = 8, height = 5.5)

Mp <- M
for (j in seq_len(ncol(Mp))) {
  v <- Mp[, j]
  med <- median(v[is.finite(v)], na.rm = TRUE)
  v[!is.finite(v)] <- med
  Mp[, j] <- v
}
pc <- prcomp(t(Mp), center = TRUE, scale. = TRUE)
pcs <- data.frame(
  sample = rownames(pc$x),
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  lineage = group
)
fwrite(pcs, file.path(out_dir, "pca_Breast_vs_Lung_scores.csv"))
p_pca <- ggplot(pcs, aes(PC1, PC2, colour = lineage, label = sample)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(size = 1.8, max.overlaps = 40, segment.size = 0.2) +
  theme_bw() +
  labs(
    title = "PCA — Breast vs Lung lines only",
    subtitle = "PCs on row-median-imputed matrix (PCA display only)"
  )
ggsave(file.path(out_dir, "pca_Breast_vs_Lung.pdf"), p_pca, width = 10, height = 7)

# Example proteins — match UniProt substrings (GeneSymbol column is sparse in this export)
marker_defs <- c(
  KRT8 = "P05787", KRT18 = "P05783", EPCAM = "P16422", CDH1 = "P12830",
  SFTPB = "P07988", SFTPC = "P11633", NKX2_1 = "NKX21_HUMAN", MUC1 = "P15941"
)
uid_all <- as.character(gm[[2]])
mk_rows <- lapply(names(marker_defs), function(nm) {
  pat <- marker_defs[[nm]]
  hi <- grep(pat, uid_all, fixed = TRUE)
  if (length(hi) == 0L) {
    return(data.table(gene = nm, pattern = pat, in_matrix = FALSE))
  }
  i <- hi[1L]
  u <- uid_all[i]
  if (!u %in% rownames(M)) {
    vec <- suppressWarnings(as.numeric(gm[i, ..samples]))
    mb <- mean(vec[seq_along(breast)], na.rm = TRUE)
    ml <- mean(vec[length(breast) + seq_along(lung)], na.rm = TRUE)
    return(data.table(
      gene = nm, pattern = pat, in_matrix = TRUE, in_limma_filtered = FALSE,
      UniProtID = u, mean_Breast = mb, mean_Lung = ml
    ))
  }
  ii <- match(u, rownames(M))
  tt_hit <- match(u, tt$UniProtID)
  mb <- mean(M[ii, breast], na.rm = TRUE)
  ml <- mean(M[ii, lung], na.rm = TRUE)
  if (is.na(tt_hit)) {
    return(data.table(
      gene = nm, pattern = pat, in_matrix = TRUE, in_limma_filtered = FALSE,
      UniProtID = u, mean_Breast = mb, mean_Lung = ml
    ))
  }
  row <- tt[tt_hit]
  data.table(
    gene = nm, pattern = pat, in_matrix = TRUE, in_limma_filtered = TRUE,
    UniProtID = u, mean_Breast = mb, mean_Lung = ml,
    logFC = row$logFC, adj.P.Val = row$adj.P.Val, P.Value = row$P.Value
  )
})
mk <- rbindlist(mk_rows, fill = TRUE)
fwrite(mk, file.path(out_dir, "lineage_example_markers.csv"))

n_sig_05 <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
n_sig_01 <- sum(tt$adj.P.Val < 0.1, na.rm = TRUE)

writeLines(c(
  "CCLE lineage DA — first pass (Breast vs Lung)",
  "=============================================",
  "",
  "Contrast: Breast - Lung (positive logFC = higher in breast cell lines).",
  paste("Samples:", length(breast), "breast lines x", length(lung), "lung lines."),
  paste("Genes after presence filter:", nrow(M)),
  paste("Features FDR < 0.05:", n_sig_05, "| FDR < 0.1:", n_sig_01),
  "",
  "Why this contrast (see ccle_lineage_contrast_choice.txt):",
  "  Stronger TMT mixture overlap between breast and lung than breast vs kidney,",
  "  reducing pure batch confounding for an exploratory lineage benchmark.",
  "",
  "Benchmark usability:",
  "  - Use for method sanity (effect sizes, volcano shape, PCA separation),",
  "    not for definitive tissue biology.",
  "  - n is unequal (30 vs 76); variance structure may differ by group size.",
  "  - No replication within lineage; interpret conservatively.",
  "",
  "Outputs:",
  "  DA_Breast_vs_Lung_limma.csv, volcano_Breast_vs_Lung.pdf,",
  "  pca_Breast_vs_Lung.pdf, pca_Breast_vs_Lung_scores.csv, lineage_example_markers.csv"
), file.path(out_dir, "README.txt"))

message("Done. Outputs in ", out_dir)
