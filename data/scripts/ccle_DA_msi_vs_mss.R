#!/usr/bin/env Rscript
# MSI vs MSS on CCLE gene_matrix — limma with DepMap lineage + sex covariates when feasible.
#
# Priority:
#   1) ~ msi_f + lineage + sex_f  (rare lineages collapsed)
#   2) ~ msi_f + lineage2 + sex_f (stronger collapse, min 10 per lineage level)
#   3) ~ msi_f + lineage2          (drop sex if needed)
#   4) ~ msi_f + sex_f
#   5) ~ msi_f
#   6) Fallback: two-group MSI - MSS contrast (no covariates)
#
# Run: Rscript --vanilla data/scripts/ccle_DA_msi_vs_mss.R

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
gm_path <- file.path(root, "data", "results", "CCLE", "gene_matrix.csv")
out_dir <- file.path(root, "data", "results", "CCLE", "msi_vs_mss", "DA_MSI_vs_MSS")

if (!file.exists(map_path)) {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  map_path <- file.path(root, "data", "results", "CCLE", "msi_vs_mss", "ccle_msi_label_mapping.csv")
  gm_path <- file.path(root, "data", "results", "CCLE", "gene_matrix.csv")
  out_dir <- file.path(root, "data", "results", "CCLE", "msi_vs_mss", "DA_MSI_vs_MSS")
}
if (!file.exists(map_path)) stop("Run build_ccle_msi_labels.py first; missing ", map_path)
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

samples <- lab$proteomics_column_name
msi_cols <- lab[msi_label == "MSI", proteomics_column_name]
mss_cols <- lab[msi_label == "MSS", proteomics_column_name]

message("Reading gene matrix ...")
gm <- fread(gm_path, showProgress = FALSE)
miss <- setdiff(samples, names(gm))
if (length(miss)) stop("Columns missing from gene_matrix: ", paste(head(miss, 10), collapse = ", "))

M <- as.matrix(gm[, ..samples])
storage.mode(M) <- "numeric"
rownames(M) <- gm[[2]]
M[!is.finite(M)] <- NA

min_msi <- max(3L, min(10L, length(msi_cols)))
min_mss <- max(3L, min(10L, length(mss_cols)))
ok <- rowSums(is.finite(M[, msi_cols, drop = FALSE])) >= min_msi &
  rowSums(is.finite(M[, mss_cols, drop = FALSE])) >= min_mss
M <- M[ok, , drop = FALSE]
message("Genes after presence filter (>=", min_msi, " MSI and >=", min_mss, " MSS lines): ", nrow(M))

sym_col <- as.character(gm[[1]])
uid_col <- as.character(gm[[2]])

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
  tt <- as.data.table(topTable(fit, coef = 1, number = Inf, sort.by = "P"), keep.rownames = "UniProtID")
  tt[, GeneSymbol := sym_col[match(UniProtID, uid_col)]]
  coef_name <- "contrast_MSI_minus_MSS"
} else {
  fit <- lmFit(M, design)
  fit <- eBayes(fit)
  if (!coef_name %in% colnames(design)) {
    stop("msi_fMSI not in design columns: ", paste(colnames(design), collapse = ", "))
  }
  tt <- as.data.table(topTable(fit, coef = coef_name, number = Inf, sort.by = "P"), keep.rownames = "UniProtID")
  tt[, GeneSymbol := sym_col[match(UniProtID, uid_col)]]
}

fwrite(tt, file.path(out_dir, "DA_MSI_vs_MSS_limma.csv"))

# Save design with sample IDs
des_out <- data.table(sample = samples, design)
fwrite(des_out, file.path(out_dir, "design_matrix_used.csv"))

pheno_out <- lab[, .(
  sample = proteomics_column_name,
  msi_label, msi_f,
  depmap_oncotree_lineage,
  lineage_f, lineage_f2,
  sex_norm, sex_f,
  MSIScore, depmap_model_id
)]
fwrite(pheno_out, file.path(out_dir, "sample_phenotype_for_DA.csv"))

diag_lines <- c(
  "Design diagnostics — CCLE MSI vs MSS",
  "====================================",
  paste("Selected model:", model_label),
  paste("Fallback two-group:", fallback),
  "",
  "Tried models (first OK wins):",
  paste("  [1] MSI+lineage(min5)+Sex: OK=", t1$ok, " rank=", t1$rank, "/", t1$ncol, if (!is.null(t1$ne)) paste(" NE:", paste(t1$ne, collapse = ",")) else ""),
  paste("  [2] MSI+lineage(min10)+Sex: OK=", t2$ok, " rank=", t2$rank, "/", t2$ncol, if (!is.null(t2$ne)) paste(" NE:", paste(t2$ne, collapse = ",")) else ""),
  paste("  [3] MSI+lineage(min10): OK=", t3$ok, " rank=", t3$rank, "/", t3$ncol),
  paste("  [4] MSI+Sex: OK=", t4$ok, " rank=", t4$rank, "/", t4$ncol),
  paste("  [5] MSI only: OK=", t5$ok, " rank=", t5$rank, "/", t5$ncol),
  "",
  "Covariates source: DepMap Model.csv (OncotreeLineage, Sex).",
  "MSI label: DepMap OmicsGlobalSignatures MSIScore > 20.",
  "",
  "Comparison to CCLE proteomics paper:",
  "  - Paper: modeled MSI with tissue and sex covariates.",
  "  - This workflow: limma with MSI vs MSS coefficient + DepMap lineage (tissue proxy) + sex when identifiable.",
  "  - Tissue: we use OncotreeLineage (not CCLE 'Tissue of Origin' from sample_info) for alignment with DepMap.",
  "  - Sex: included when full-rank; dropped if needed for estimability.",
  "",
  "vs simple Breast-vs-Lung lineage contrast:",
  "  - MSI-vs-MSS (adjusted) targets microsatellite biology + covariates — closer to paper MSI analysis.",
  "  - Breast-vs-Lung is organ comparison without MSI; different scientific question.",
  "  - Neither is a perfect main benchmark; MSI here uses DepMap score, not clinical MSI.",
  ""
)
writeLines(diag_lines, file.path(out_dir, "design_model_summary.txt"))

# Volcano
tt$sig <- tt$adj.P.Val < 0.1 & abs(tt$logFC) > 0.5
sub <- paste0(
  "n=", length(msi_cols), " MSI, ", length(mss_cols), " MSS | ",
  if (!fallback) paste("coef", coef_name) else "fallback MSI-MSS"
)
p_volc <- ggplot(tt, aes(logFC, -log10(P.Value))) +
  geom_point(aes(colour = sig), alpha = 0.2, size = 0.35) +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey50")) +
  theme_bw() +
  labs(
    title = "CCLE MSI vs MSS (limma)",
    subtitle = sub,
    x = if (!fallback) paste("log2 FC:", coef_name, "(MSI vs MSS, adjusted)") else "logFC (MSI - MSS)",
    y = "-log10 p"
  )
ggsave(file.path(out_dir, "volcano_MSI_vs_MSS.pdf"), p_volc, width = 8, height = 5.5)

# PCA
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
  sample = samples,
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  msi_group = msi_fac,
  lineage = lab$depmap_oncotree_lineage,
  sex = lab$sex_norm
)
fwrite(pcs, file.path(out_dir, "pca_MSI_vs_MSS_scores.csv"))
p_pca <- ggplot(pcs, aes(PC1, PC2, colour = msi_group, label = sample)) +
  geom_point(size = 1.8, alpha = 0.75) +
  geom_text_repel(size = 1.5, max.overlaps = 12, segment.size = 0.15) +
  theme_bw() +
  labs(
    title = "PCA — MSI vs MSS (labeled lines)",
    subtitle = "Imputed matrix for display only"
  )
ggsave(file.path(out_dir, "pca_MSI_vs_MSS.pdf"), p_pca, width = 11, height = 8)

n_fdr <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
writeLines(
  c(
    "CCLE MSI vs MSS — DA README",
    "===========================",
    "",
    "Matrix: data/results/CCLE/gene_matrix.csv",
    "Mapping: data/results/CCLE/msi_vs_mss/ccle_msi_label_mapping.csv",
    "",
    "Presence filter:",
    paste("  Per gene: >= ", min_msi, " observed in MSI lines AND >= ", min_mss, " in MSS lines.", sep = ""),
    paste("  Genes tested:", nrow(M)),
    "",
    "Model:",
    paste(" ", model_label),
    paste("  Coefficient in topTable:", coef_name),
    "",
    paste("Samples: ", length(msi_cols), " MSI; ", length(mss_cols), " MSS."),
    paste("Genes with FDR < 0.05:", n_fdr),
    "",
    "Benchmark usability:",
    "  Exploratory / methods comparison with covariate adjustment (paper-like): appropriate.",
    "  Definitive MSI proteome biology: not without clinical MSI labels and/or Bowel restriction.",
    "",
    "See design_model_summary.txt for formula diagnostics and comparison to CCLE paper.",
    ""
  ),
  file.path(out_dir, "README.txt")
)

message("Selected: ", model_label)
message("Done. Outputs in ", out_dir)
