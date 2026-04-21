#!/usr/bin/env Rscript
# Luminal vs Basal (v1: 4+4 breast lines) on CCLE Table S2 protein matrix
# (CCLE-code columns from table_s2_protein_matrix_cclecode_matched.csv.gz).
#
# Column names are CCLE codes (e.g. MCF7_BREAST), not display Cell Line names.
# Same contrasts as ccle_DA_luminal_basal_v1.R.
#
# Run from repo root:
#   Rscript --vanilla data/scripts/ccle_DA_luminal_basal_table_s2.R
#
# Requires: table built by ccle_sum_table_s2_to_matched_matrix.py

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Install ggrepel: install.packages(\"ggrepel\")")
  library(ggrepel)
})

root <- getwd()
gm_path <- file.path(root, "data", "results", "CCLE", "ccle_sum", "table_s2_protein_matrix_cclecode_matched.csv.gz")
out_dir <- file.path(root, "data", "results", "CCLE", "ccle_sum", "DA_luminal_vs_basal_table_s2")

if (!file.exists(gm_path)) {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  gm_path <- file.path(root, "data", "results", "CCLE", "ccle_sum", "table_s2_protein_matrix_cclecode_matched.csv.gz")
  out_dir <- file.path(root, "data", "results", "CCLE", "ccle_sum", "DA_luminal_vs_basal_table_s2")
}
if (!file.exists(gm_path)) {
  stop("Missing ", gm_path, " — run: .venv/bin/python data/scripts/ccle_sum_table_s2_to_matched_matrix.py")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# v1 lines: display name -> CCLE Code (from sample_info_ccle.csv)
luminal_ccle <- c(
  "MCF7_BREAST", "T47D_BREAST", "CAMA1_BREAST", "ZR751_BREAST"
)
basal_ccle <- c(
  "HCC1806_BREAST", "HCC1143_BREAST", "HCC70_BREAST", "MDAMB468_BREAST"
)
luminal <- luminal_ccle
basal <- basal_ccle
samples <- c(luminal, basal)

message("Reading Table S2 matrix ...")
gm <- fread(gm_path, showProgress = FALSE)
meta <- c("Protein_Id", "Gene_Symbol", "Description", "Group_ID", "Uniprot", "Uniprot_Acc")
miss <- setdiff(samples, names(gm))
if (length(miss)) {
  stop(
    "Missing CCLE columns (not in MSI-matched Table S2 subset): ",
    paste(miss, collapse = ", "),
    "\nRebuild matrix including these lines or use gene_matrix DA."
  )
}

uid <- gm[["Protein_Id"]]
M <- as.matrix(gm[, ..samples])
storage.mode(M) <- "numeric"
rownames(M) <- uid
M[!is.finite(M)] <- NA

ok <- rowSums(is.finite(M[, luminal, drop = FALSE])) >= 2L &
  rowSums(is.finite(M[, basal, drop = FALSE])) >= 2L
M <- M[ok, , drop = FALSE]
message("Proteins after filter (>=2/4 observed per group): ", nrow(M))

group <- factor(rep(c("Luminal", "Basal"), c(4L, 4L)), levels = c("Luminal", "Basal"))
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Luminal", "Basal")
fit <- lmFit(M, design)
contr <- makeContrasts(Luminal_vs_Basal = Luminal - Basal, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
tt <- as.data.table(topTable(fit2, coef = 1, number = Inf, sort.by = "P"), keep.rownames = "Protein_Id")
fwrite(tt, file.path(out_dir, "DA_luminal_vs_basal_limma.csv"))

tt$sig <- tt$adj.P.Val < 0.1 & abs(tt$logFC) > 0.5
p_volc <- ggplot(tt, aes(logFC, -log10(P.Value))) +
  geom_point(aes(colour = sig), alpha = 0.35, size = 0.6) +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey50")) +
  theme_bw() +
  labs(
    title = "CCLE Luminal vs Basal — Table S2 (4 vs 4 lines)",
    subtitle = "limma on paper Table S2 normalized protein matrix (CCLE-code columns)",
    x = "logFC (Luminal - Basal)",
    y = "-log10 p"
  )
ggsave(file.path(out_dir, "volcano_luminal_vs_basal.pdf"), p_volc, width = 7, height = 5)

Mp <- M
for (j in seq_len(ncol(Mp))) {
  v <- Mp[, j]
  med <- median(v[is.finite(v)], na.rm = TRUE)
  v[!is.finite(v)] <- med
  Mp[, j] <- v
}
pc <- prcomp(t(Mp), center = TRUE, scale. = TRUE)
pcs <- data.frame(
  sample_ccle_code = rownames(pc$x),
  display_name = c("MCF7", "T-47D", "CAMA-1", "ZR-75-1", "HCC 1806", "HCC1143", "HCC70", "MDA-MB-468"),
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  group = group
)
fwrite(pcs, file.path(out_dir, "ccle_pca_scores.csv"))
p_pca <- ggplot(pcs, aes(PC1, PC2, colour = group, label = display_name)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size = 2.5, max.overlaps = 20) +
  theme_bw() +
  labs(
    title = "PCA — Luminal vs Basal (Table S2)",
    subtitle = "PCs on imputed matrix (display only)"
  )
ggsave(file.path(out_dir, "ccle_luminal_basal_pca.pdf"), p_pca, width = 7, height = 5)

marker_patterns <- c(
  ESR1 = "P03372", GATA3 = "P23771", FOXA1 = "P55317", KRT18 = "P05783", PGR = "P06454",
  KRT5 = "P13647", KRT14 = "P02533", KRT17 = "Q04695", EGFR = "P00533", FOXC1 = "Q12948"
)
uid_full <- gm[["Protein_Id"]]
mk_rows <- lapply(names(marker_patterns), function(nm) {
  pat <- marker_patterns[[nm]]
  hi <- grep(pat, uid_full, fixed = FALSE)[1]
  if (is.na(hi)) {
    return(data.table(gene = nm, pattern = pat, in_matrix = FALSE))
  }
  u <- uid_full[hi]
  vec <- suppressWarnings(as.numeric(gm[hi, ..samples]))
  mL <- mean(vec[1:4], na.rm = TRUE)
  mB <- mean(vec[5:8], na.rm = TRUE)
  hit_tt <- match(u, tt$Protein_Id)
  if (is.na(hit_tt)) {
    return(data.table(
      gene = nm, pattern = pat, in_matrix = TRUE, Protein_Id = u,
      mean_Luminal = mL, mean_Basal = mB, delta_mean = mL - mB,
      in_limma_filtered = FALSE, logFC = NA_real_, adj.P.Val = NA_real_
    ))
  }
  row <- tt[hit_tt]
  data.table(
    gene = nm, pattern = pat, in_matrix = TRUE, Protein_Id = u,
    mean_Luminal = mL, mean_Basal = mB, delta_mean = mL - mB,
    in_limma_filtered = TRUE,
    logFC = row$logFC, adj.P.Val = row$adj.P.Val, P.Value = row$P.Value
  )
})
mk <- rbindlist(mk_rows, fill = TRUE)
fwrite(mk, file.path(out_dir, "canonical_markers_check.csv"))

writeLines(c(
  "CCLE Luminal vs Basal — Table S2 (exploratory v1)",
  "=================================================",
  "Matrix: data/results/CCLE/ccle_sum/table_s2_protein_matrix_cclecode_matched.csv.gz",
  "  (MSI-matched CCLE subset; same 8 breast lines as gene_matrix v1 if present).",
  "CCLE code columns: MCF7_BREAST, T47D_BREAST, CAMA1_BREAST, ZR751_BREAST vs",
  "  HCC1806_BREAST, HCC1143_BREAST, HCC70_BREAST, MDAMB468_BREAST.",
  "Design: limma Luminal - Basal; filter >=2/4 observed per group.",
  "Caveats: same as gene_matrix v1 (no replication; batch possible).",
  "",
  "Compare: data/results/CCLE/DA_luminal_vs_basal/ (gene_matrix)."
), file.path(out_dir, "README.txt"))

message("Done. Outputs in ", out_dir)
