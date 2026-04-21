#!/usr/bin/env Rscript
# =============================================================================
# Subtype differential abundance — CPTAC breast (PDC000120)
# Uses DA_subtype_tumor_only.csv (one row per matrix column, Tumor + PAM50).
# Design: ~ 0 + PAM50; contrasts e.g. Basal vs LumA, LumA vs LumB.
# =============================================================================
# Prerequisite: run build_PDC000120_subtype_mapping.py to generate
#   results/PDC000120/DA_subtype_tumor_only.csv and DA_subtype_counts.csv
# Usage (from data/): Rscript scripts/DA_subtype_CPTAC_breast.R
# Outputs: results/PDC000120/subtype_DA_*.csv, subtype_volcano_*.pdf
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", update = FALSE, ask = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
})
library(limma)
library(data.table)

# Paths
DATA_DIR <- if (file.exists(file.path(getwd(), "results", "PDC000120", "gene_matrix.csv"))) getwd() else file.path(getwd(), "data")
setwd(DATA_DIR)
RESULTS_DIR <- file.path(DATA_DIR, "results", "PDC000120")
GENE_MATRIX_PATH <- file.path(RESULTS_DIR, "gene_matrix.csv")
SUBTYPE_ANNOT_PATH <- file.path(RESULTS_DIR, "DA_subtype_tumor_only.csv")
ANNOT_PATH <- file.path(RESULTS_DIR, "annotation_filled_corrected.csv")
OUT_DIR <- RESULTS_DIR
dir.create(OUT_DIR, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load gene matrix and subtype design
# -----------------------------------------------------------------------------
message("Step 1: Load matrix and subtype annotation")
gm <- fread(GENE_MATRIX_PATH)
rownames(gm) <- gm$GeneSymbol
gm$UniProtID <- NULL
all_cols <- setdiff(names(gm), "GeneSymbol")
mat_full <- as.matrix(gm[, ..all_cols])
rownames(mat_full) <- gm$GeneSymbol

if (!file.exists(SUBTYPE_ANNOT_PATH)) {
  stop("Run build_PDC000120_subtype_mapping.py first to create ", SUBTYPE_ANNOT_PATH)
}
design_dt <- fread(SUBTYPE_ANNOT_PATH)
setnames(design_dt, trimws(gsub("^\uFEFF", "", names(design_dt))))
id_col <- names(design_dt)[grepl("matrix_sample_id|bioreplicate", names(design_dt), ignore.case = TRUE)][1]
pam50_col <- names(design_dt)[grepl("pam50", names(design_dt), ignore.case = TRUE)][1]
if (is.na(id_col)) id_col <- "matrix_sample_id"
if (is.na(pam50_col)) pam50_col <- "pam50"

# Subset to samples in design (tumor + PAM50)
samples <- design_dt[[id_col]]
samples <- samples[!is.na(samples) & nzchar(trimws(samples))]
in_matrix <- samples %in% colnames(mat_full)
if (!all(in_matrix)) {
  message("  Dropping ", sum(!in_matrix), " design samples not in gene matrix")
  samples <- samples[in_matrix]
}
mat <- mat_full[, colnames(mat_full) %in% samples, drop = FALSE]
design_dt <- design_dt[get(id_col) %in% colnames(mat)]
# Align order to matrix columns
idx <- match(colnames(mat), design_dt[[id_col]])
design_dt <- design_dt[idx]
pam50 <- factor(trimws(design_dt[[pam50_col]]))
message("  Matrix: ", nrow(mat), " proteins x ", ncol(mat), " tumor samples")
message("  PAM50: ", paste(levels(pam50), collapse = ", "))

# -----------------------------------------------------------------------------
# 2. Design and contrasts
# -----------------------------------------------------------------------------
message("Step 2: Design ~ 0 + PAM50")
design <- model.matrix(~ 0 + pam50)
colnames(design) <- gsub("pam50", "", colnames(design))

# Contrasts (only if both groups have at least 3 samples)
# design colnames are like pam50Basal, pam50LumA (no space)
levs <- levels(pam50)
design_levs <- colnames(design)
n_per <- setNames(as.integer(table(pam50)), levs)
contrast_list <- list(
  "Basal_vs_LumA"    = c("Basal", "LumA"),
  "Basal_vs_LumB"    = c("Basal", "LumB"),
  "Basal_vs_Her2"    = c("Basal", "Her2"),
  "LumA_vs_LumB"     = c("LumA", "LumB"),
  "LumA_vs_Her2"     = c("LumA", "Her2"),
  "LumB_vs_Her2"     = c("LumB", "Her2"),
  "Normal-like_vs_LumA" = c("Normal-like", "LumA")
)
valid_contrasts <- list()
for (nm in names(contrast_list)) {
  pair <- contrast_list[[nm]]
  # Find design columns that end with or equal these level names
  c1 <- design_levs[design_levs == pair[1] | endsWith(design_levs, pair[1])][1]
  c2 <- design_levs[design_levs == pair[2] | endsWith(design_levs, pair[2])][1]
  if (is.na(c1) || is.na(c2)) next
  n1 <- n_per[pair[1]]
  n2 <- n_per[pair[2]]
  if (is.na(n1)) n1 <- 0
  if (is.na(n2)) n2 <- 0
  if (n1 >= 3 && n2 >= 3)
    valid_contrasts[[nm]] <- setNames(c(1, -1), c(c1, c2))
}
if (length(valid_contrasts) == 0 && length(levs) >= 2) {
  c1 <- design_levs[1]
  c2 <- design_levs[2]
  valid_contrasts[[paste0(levs[1], "_vs_", levs[2])]] <- setNames(c(1, -1), c(c1, c2))
}

# -----------------------------------------------------------------------------
# 3. limma fit and extract results per contrast
# -----------------------------------------------------------------------------
message("Step 3: limma fit and contrasts")
fit <- lmFit(mat, design)
out_tables <- list()
for (cn in names(valid_contrasts)) {
  vec <- valid_contrasts[[cn]]
  nms <- names(vec)
  if (!all(nms %in% colnames(design))) next
  cont <- makeContrasts(
    contrasts = paste0(nms[1], " - ", nms[2]),
    levels = design
  )
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
  res$GeneSymbol <- rownames(res)
  res$contrast <- cn
  out_tables[[cn]] <- as.data.frame(res)
}

# -----------------------------------------------------------------------------
# 4. Write results and simple volcano
# -----------------------------------------------------------------------------
message("Step 4: Write results")
all_res <- rbindlist(out_tables)
fwrite(all_res, file.path(OUT_DIR, "subtype_DA_limma_all.csv"))

# One wide file per contrast (GeneSymbol, logFC, P.Value, adj.P.Val)
for (cn in names(out_tables)) {
  r <- out_tables[[cn]]
  r <- r[, .(GeneSymbol, logFC, P.Value, adj.P.Val)]
  setnames(r, c("logFC", "P.Value", "adj.P.Val"), paste0(c("logFC_", "P_", "FDR_"), cn))
  fwrite(r, file.path(OUT_DIR, paste0("subtype_DA_", gsub("_vs_", "_vs_", cn, fixed = TRUE), ".csv")))
}

# Simple volcano (first contrast)
if (length(out_tables) > 0 && requireNamespace("graphics", quietly = TRUE)) {
  r <- out_tables[[1]]
  pdf(file.path(OUT_DIR, paste0("subtype_volcano_", names(out_tables)[1], ".pdf")), width = 6, height = 5)
  plot(r$logFC, -log10(r$P.Value + 1e-20), pch = 20, col = ifelse(r$adj.P.Val < 0.05, "red", "grey50"),
       xlab = "log2 FC", ylab = "-log10(P)", main = names(out_tables)[1])
  abline(h = -log10(0.05), lty = 2, col = "grey")
  dev.off()
  message("  Wrote ", file.path(OUT_DIR, paste0("subtype_volcano_", names(out_tables)[1], ".pdf")))
}

message("Done. Subtype DA results in ", OUT_DIR)
