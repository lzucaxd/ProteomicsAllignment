#!/usr/bin/env Rscript
# Same MSI vs MSS limma workflow as ccle_DA_msi_vs_mss.R, but expression from
# CCLE paper Table S2 (data/ccle_sum) collapsed to CCLE codes — see
# data/scripts/ccle_sum_table_s2_to_matched_matrix.py
#
# Join: mapping ccle_code_from_sample_info == matrix column names.
#
# Run (after Python matrix build + build_ccle_msi_labels.py):
#   Rscript --vanilla data/scripts/ccle_DA_msi_vs_mss_table_s2.R

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Install ggrepel: install.packages(\"ggrepel\")")
  library(ggrepel)
})

root <- getwd()
map_path <- file.path(root, "data", "results", "CCLE", "msi_vs_mss", "ccle_msi_label_mapping.csv")
gm_path <- file.path(root, "data", "results", "CCLE", "ccle_sum", "table_s2_protein_matrix_cclecode_matched.csv.gz")
out_dir <- file.path(root, "data", "results", "CCLE", "msi_vs_mss", "DA_MSI_vs_MSS_table_s2")

if (!file.exists(map_path)) {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  map_path <- file.path(root, "data", "results", "CCLE", "msi_vs_mss", "ccle_msi_label_mapping.csv")
  gm_path <- file.path(root, "data", "results", "CCLE", "ccle_sum", "table_s2_protein_matrix_cclecode_matched.csv.gz")
  out_dir <- file.path(root, "data", "results", "CCLE", "msi_vs_mss", "DA_MSI_vs_MSS_table_s2")
}
if (!file.exists(gm_path)) {
  stop("Missing Table S2 matrix. Run: .venv/bin/python data/scripts/ccle_sum_table_s2_to_matched_matrix.py")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

map <- fread(map_path, showProgress = FALSE)
lab <- map[match_status == "matched" & msi_label %in% c("MSI", "MSS")]
if (nrow(lab) < 10L) stop("Too few labeled samples: ", nrow(lab))

norm_sex <- function(s) {
  s <- trimws(as.character(s))
  s[s %in% c("", "NA")] <- "Unknown"
  s[!s %in% c("Male", "Female", "Unknown")] <- "Unknown"
  s
}

collapse_lineage <- function(lin, min_n = 5L) {
  lin <- as.character(lin)
  tab <- table(lin)
  rare <- names(tab)[tab < min_n]
  lin[lin %in% rare] <- "Other_rare"
  factor(lin)
}

lab[, sex_norm := norm_sex(depmap_sex)]
lab[, sex_f := factor(sex_norm, levels = c("Female", "Male", "Unknown"))]
lab[, msi_f := factor(msi_label, levels = c("MSS", "MSI"))]
lab[, lineage_f := collapse_lineage(depmap_oncotree_lineage, 5L)]
lab[, lineage_f2 := collapse_lineage(depmap_oncotree_lineage, 10L)]

message("Reading Table S2 matrix ...")
gm <- fread(gm_path, showProgress = FALSE)
meta <- c("Protein_Id", "Gene_Symbol", "Description", "Group_ID", "Uniprot", "Uniprot_Acc")
samples <- setdiff(names(gm), meta)
lab <- lab[ccle_code_from_sample_info %in% samples]
lab <- lab[order(match(ccle_code_from_sample_info, samples))]
if (nrow(lab) != length(samples)) {
  stop("Phenotype rows (", nrow(lab), ") != matrix columns (", length(samples), ") — duplicate CCLE?")
}

msi_cols <- lab[msi_label == "MSI", ccle_code_from_sample_info]
mss_cols <- lab[msi_label == "MSS", ccle_code_from_sample_info]

M <- as.matrix(gm[, ..samples])
storage.mode(M) <- "numeric"
rownames(M) <- gm[["Protein_Id"]]
M[!is.finite(M)] <- NA

sym_col <- as.character(gm[["Gene_Symbol"]])
uid_col <- as.character(gm[["Protein_Id"]])

min_msi <- max(3L, min(10L, length(msi_cols)))
min_mss <- max(3L, min(10L, length(mss_cols)))
ok <- rowSums(is.finite(M[, msi_cols, drop = FALSE])) >= min_msi &
  rowSums(is.finite(M[, mss_cols, drop = FALSE])) >= min_mss
M <- M[ok, , drop = FALSE]
message("Proteins after filter: ", nrow(M))

try_fit <- function(formula, lab_dt, label_txt) {
  design <- model.matrix(formula, data = lab_dt)
  ne <- nonEstimable(design)
  rk <- qr(design)$rank
  ok <- is.null(ne) && rk == ncol(design)
  list(design = design, ok = ok, label = label_txt, ne = ne, rank = rk, ncol = ncol(design))
}

t1 <- try_fit(~ msi_f + lineage_f + sex_f, lab, "MSI + lineage (min 5) + Sex")
t2 <- try_fit(~ msi_f + lineage_f2 + sex_f, lab, "MSI + lineage (min 10) + Sex")
t3 <- try_fit(~ msi_f + lineage_f2, lab, "MSI + lineage (min 10), no Sex")
t4 <- try_fit(~ msi_f + sex_f, lab, "MSI + Sex only")
t5 <- try_fit(~ msi_f, lab, "MSI vs MSS intercept-only (no covariates)")

design <- NULL
model_label <- NULL
coef_name <- "msi_fMSI"

if (t1$ok) {
  design <- t1$design
  model_label <- t1$label
} else if (t2$ok) {
  design <- t2$design
  model_label <- t2$label
} else if (t3$ok) {
  design <- t3$design
  model_label <- t3$label
} else if (t4$ok) {
  design <- t4$design
  model_label <- t4$label
} else if (t5$ok) {
  design <- t5$design
  model_label <- t5$label
}

fallback <- FALSE
if (is.null(design)) {
  fallback <- TRUE
  model_label <- "FALLBACK: two-group contrast MSI - MSS (covariate designs rank-deficient or non-estimable)"
  group <- factor(c(rep("MSI", length(msi_cols)), rep("MSS", length(mss_cols))), levels = c("MSI", "MSS"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("MSI", "MSS")
  fit <- lmFit(M, design)
  contr <- makeContrasts(MSI_vs_MSS = MSI - MSS, levels = design)
  fit <- contrasts.fit(fit, contr)
  fit <- eBayes(fit)
  tt <- as.data.table(topTable(fit, coef = 1, number = Inf, sort.by = "P"), keep.rownames = "Protein_Id")
  tt[, GeneSymbol := sym_col[match(Protein_Id, uid_col)]]
  coef_name <- "contrast_MSI_minus_MSS"
} else {
  fit <- lmFit(M, design)
  fit <- eBayes(fit)
  if (!coef_name %in% colnames(design)) {
    stop("msi_fMSI not in design columns: ", paste(colnames(design), collapse = ", "))
  }
  tt <- as.data.table(topTable(fit, coef = coef_name, number = Inf, sort.by = "P"), keep.rownames = "Protein_Id")
  tt[, GeneSymbol := sym_col[match(Protein_Id, uid_col)]]
}

fwrite(tt, file.path(out_dir, "DA_MSI_vs_MSS_limma.csv"))

des_out <- data.table(sample_ccle_code = samples, design)
fwrite(des_out, file.path(out_dir, "design_matrix_used.csv"))

pheno_out <- lab[, .(
  ccle_code = ccle_code_from_sample_info,
  cell_line_display = proteomics_column_name,
  msi_label, msi_f,
  depmap_oncotree_lineage,
  lineage_f, lineage_f2,
  sex_norm, sex_f,
  MSIScore, depmap_model_id
)]
fwrite(pheno_out, file.path(out_dir, "sample_phenotype_for_DA.csv"))

diag_lines <- c(
  "Design diagnostics — CCLE Table S2 MSI vs MSS",
  "=============================================",
  paste("Expression matrix: data/results/CCLE/ccle_sum/table_s2_protein_matrix_cclecode_matched.csv.gz"),
  paste("Selected model:", model_label),
  paste("Fallback two-group:", fallback),
  "",
  "Tried models (first OK wins):",
  paste("  [1] MSI+lineage(min5)+Sex: OK=", t1$ok, " rank=", t1$rank, "/", t1$ncol),
  paste("  [2] MSI+lineage(min10)+Sex: OK=", t2$ok, " rank=", t2$rank, "/", t2$ncol),
  paste("  [3] MSI+lineage(min10): OK=", t3$ok, " rank=", t3$rank, "/", t3$ncol),
  paste("  [4] MSI+Sex: OK=", t4$ok, " rank=", t4$rank, "/", t4$ncol),
  paste("  [5] MSI only: OK=", t5$ok, " rank=", t5$rank, "/", t5$ncol),
  "",
  "Table S2: CCLE supplementary normalized protein expression; columns collapsed from *_TenPxNN to CCLE code.",
  ""
)
writeLines(diag_lines, file.path(out_dir, "design_model_summary.txt"))

tt$sig <- tt$adj.P.Val < 0.1 & abs(tt$logFC) > 0.5
sub <- paste0(
  "Table S2 | n=", length(msi_cols), " MSI, ", length(mss_cols), " MSS | ",
  if (!fallback) paste("coef", coef_name) else "fallback MSI-MSS"
)
p_volc <- ggplot(tt, aes(logFC, -log10(P.Value))) +
  geom_point(aes(colour = sig), alpha = 0.2, size = 0.35) +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey50")) +
  theme_bw() +
  labs(
    title = "CCLE Table S2 — MSI vs MSS (limma)",
    subtitle = sub,
    x = if (!fallback) paste("log2 FC:", coef_name) else "logFC (MSI - MSS)",
    y = "-log10 p"
  )
ggsave(file.path(out_dir, "volcano_MSI_vs_MSS.pdf"), p_volc, width = 8, height = 5.5)

Mp <- M
for (j in seq_len(ncol(Mp))) {
  v <- Mp[, j]
  med <- median(v[is.finite(v)], na.rm = TRUE)
  v[!is.finite(v)] <- med
  Mp[, j] <- v
}
pc <- prcomp(t(Mp), center = TRUE, scale. = TRUE)
msi_fac <- factor(ifelse(samples %in% msi_cols, "MSI", "MSS"), levels = c("MSS", "MSI"))
pcs <- data.frame(
  sample_ccle_code = samples,
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  msi_group = msi_fac,
  lineage = lab$depmap_oncotree_lineage,
  sex = lab$sex_norm
)
fwrite(pcs, file.path(out_dir, "pca_MSI_vs_MSS_scores.csv"))
p_pca <- ggplot(pcs, aes(PC1, PC2, colour = msi_group, label = sample_ccle_code)) +
  geom_point(size = 1.8, alpha = 0.75) +
  geom_text_repel(size = 1.5, max.overlaps = 12, segment.size = 0.15) +
  theme_bw() +
  labs(
    title = "PCA — Table S2 (CCLE code columns)",
    subtitle = "Labels: DepMap MSI/MSS"
  )
ggsave(file.path(out_dir, "pca_MSI_vs_MSS.pdf"), p_pca, width = 11, height = 8)

n_fdr <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
writeLines(
  c(
    "CCLE MSI vs MSS — Table S2 DA README",
    "====================================",
    "",
    "Matrix: data/results/CCLE/ccle_sum/table_s2_protein_matrix_cclecode_matched.csv.gz",
    "  (from Table_S2 Protein Quant Normalized xlsx; TenPx columns collapsed to CCLE code; median if duplicate).",
    "Mapping: data/results/CCLE/msi_vs_mss/ccle_msi_label_mapping.csv (same DepMap MSI labels as gene_matrix run).",
    "",
    "Presence filter:",
    paste("  Per protein: >= ", min_msi, " observed in MSI lines AND >= ", min_mss, " in MSS lines.", sep = ""),
    paste("  Proteins tested:", nrow(M)),
    "",
    "Model:",
    paste(" ", model_label),
    paste("  Coefficient:", coef_name),
    "",
    paste("Samples: ", length(msi_cols), " MSI; ", length(mss_cols), " MSS (364 CCLE codes total)."),
    paste("Proteins with FDR < 0.05:", n_fdr),
    "",
    "Compare to data/results/CCLE/msi_vs_mss/DA_MSI_vs_MSS/ (gene_matrix from pipeline).",
    ""
  ),
  file.path(out_dir, "README.txt")
)

message("Selected: ", model_label)
message("Done. Outputs in ", out_dir)
