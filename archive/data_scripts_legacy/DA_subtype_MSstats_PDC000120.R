#!/usr/bin/env Rscript
# =============================================================================
# Subtype-aware differential abundance — CPTAC breast (PDC000120)
# MSstats-style analysis from summarized gene/protein matrix.
# Tumor-only samples; PAM50 from DA_subtype_tumor_only.csv.
# =============================================================================
# Usage (from project root or data/):
#   Rscript data/scripts/DA_subtype_MSstats_PDC000120.R
#   # or from data/: Rscript scripts/DA_subtype_MSstats_PDC000120.R
#
# If you see "cannot open file 'renv/activate.R'", run without loading profiles:
#   Rscript --vanilla scripts/DA_subtype_MSstats_PDC000120.R
#
# Outputs: results/PDC000120/DA_MSstats_<contrast>.csv, volcano_*.pdf, *_summary.txt
# =============================================================================

# -----------------------------------------------------------------------------
# CONFIG — change these to run other contrasts later
# -----------------------------------------------------------------------------
# Primary contrast: Luminal - Basal (positive log2FC => higher in Luminal); matches CCLE scripts.
contrast_name <- "Luminal_vs_Basal"
group1         <- "Luminal"   # numerator
group2         <- "Basal"     # denominator (LumA + LumB pooled to Luminal)

# Coverage filter: require non-missing in at least this fraction of samples (before DA)
pct_overall_min <- 0.35   # 35% of all samples (range 30–40%)
pct_group_min   <- 0.25   # 25% of each subtype group (range 20–30%)

# Other contrasts (uncomment and set contrast_name/group1/group2):
# contrast_name <- "Basal_vs_LumA";    group1 <- "Basal"; group2 <- "LumA"
# contrast_name <- "Her2_vs_LumA";    group1 <- "Her2";  group2 <- "LumA"
# contrast_name <- "LumA_vs_LumB";     group1 <- "LumA";  group2 <- "LumB"

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
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

# Paths: run from project root or data/
DATA_DIR <- getwd()
if (!file.exists(file.path(DATA_DIR, "results", "PDC000120", "gene_matrix.csv")))
  DATA_DIR <- file.path(getwd(), "data")
if (!file.exists(file.path(DATA_DIR, "results", "PDC000120", "gene_matrix.csv")))
  stop("Cannot find results/PDC000120/gene_matrix.csv. Run from project root or data/.")

RESULTS_DIR <- file.path(DATA_DIR, "results", "PDC000120")
GENE_MATRIX_PATH   <- file.path(RESULTS_DIR, "gene_matrix.csv")
SUBTYPE_ANNOT_PATH <- file.path(RESULTS_DIR, "DA_subtype_tumor_only.csv")
ANNOT_PATH         <- file.path(RESULTS_DIR, "annotation_filled_corrected.csv")
OUT_DIR <- RESULTS_DIR
dir.create(OUT_DIR, showWarnings = FALSE)

# Trim BOM and whitespace from headers
trim_header <- function(dt) {
  setnames(dt, trimws(gsub("^\uFEFF", "", names(dt))))
  invisible(dt)
}

# =============================================================================
# STEP 1 — Data prep
# =============================================================================
message("========== STEP 1: Data prep ==========")

# 1.1 Load subtype design
design_dt <- fread(SUBTYPE_ANNOT_PATH)
trim_header(design_dt)
id_col  <- names(design_dt)[grepl("matrix_sample_id|bioreplicate", names(design_dt), ignore.case = TRUE)][1]
pam50_col <- names(design_dt)[grepl("pam50", names(design_dt), ignore.case = TRUE)][1]
mix_col   <- names(design_dt)[grepl("mixture", names(design_dt), ignore.case = TRUE)][1]
if (is.na(id_col))    id_col    <- "matrix_sample_id"
if (is.na(pam50_col)) pam50_col <- "pam50"
if (is.na(mix_col))   mix_col   <- "mixture"

design_dt[, (pam50_col) := trimws(as.character(get(pam50_col)))]
design_dt <- design_dt[get(pam50_col) != "" & !is.na(get(pam50_col))]

# Pool LumA + LumB to Luminal when Luminal is in the contrast
if (group1 == "Luminal" || group2 == "Luminal" || contrast_name == "Luminal_vs_Basal") {
  design_dt[get(pam50_col) %in% c("LumA", "LumB"), (pam50_col) := "Luminal"]
}
design_dt <- design_dt[get(pam50_col) %in% c(group1, group2)]
if (id_col != "matrix_sample_id") setnames(design_dt, id_col, "matrix_sample_id")
if (pam50_col != "pam50") setnames(design_dt, pam50_col, "pam50")
if (mix_col != "mixture" && !is.na(mix_col)) setnames(design_dt, mix_col, "mixture")
if (!"mixture" %in% names(design_dt)) design_dt[, mixture := NA_character_]

message("  Subtype design rows (tumor, non-missing PAM50, in contrast): ", nrow(design_dt))
message("  ", group1, ": ", sum(design_dt$pam50 == group1))
message("  ", group2, ": ", sum(design_dt$pam50 == group2))

# 1.2 Load gene matrix
gm <- fread(GENE_MATRIX_PATH)
trim_header(gm)
gene_col <- names(gm)[1]
if (!identical(tolower(gene_col), "genesymbol")) {
  if ("GeneSymbol" %in% names(gm)) gene_col <- "GeneSymbol"
  else if ("genesymbol" %in% names(gm)) gene_col <- "genesymbol"
}
non_sample_cols <- c(gene_col, "UniProtID", "uniprotid")
non_sample_cols <- intersect(non_sample_cols, names(gm))
all_matrix_cols <- setdiff(names(gm), non_sample_cols)

# Case-insensitive match: design matrix_sample_id (often lowercase) -> actual matrix column
design_ids <- design_dt$matrix_sample_id
design_lower <- tolower(trimws(design_ids))
col_lower <- tolower(trimws(all_matrix_cols))
match_idx <- match(design_lower, col_lower)
n_miss <- sum(is.na(match_idx))
if (n_miss > 0) {
  message("  WARNING: ", n_miss, " design sample(s) not found in gene matrix (case-insensitive). Dropping them.")
  design_dt <- design_dt[!is.na(match_idx)]
  match_idx <- match_idx[!is.na(match_idx)]
}
mat_cols <- all_matrix_cols[match_idx]
design_dt[, matrix_col := mat_cols]

# Build matrix: rows = genes, columns = design samples (order preserved)
gm_mat <- as.matrix(gm[, ..mat_cols])
rownames(gm_mat) <- gm[[gene_col]]
# Coerce to numeric; empty -> NA
storage.mode(gm_mat) <- "double"
n_genes <- nrow(gm_mat)
message("  Genes/proteins: ", n_genes)
message("  ", group1, " samples: ", sum(design_dt$pam50 == group1))
message("  ", group2, " samples: ", sum(design_dt$pam50 == group2))

# Mixture distribution per subtype
if ("mixture" %in% names(design_dt) && any(!is.na(design_dt$mixture))) {
  mix_tab <- design_dt[, .N, by = .(pam50, mixture)]
  n_mix_group1 <- mix_tab[pam50 == group1, uniqueN(mixture)]
  n_mix_group2 <- mix_tab[pam50 == group2, uniqueN(mixture)]
  message("  Unique mixtures in ", group1, ": ", n_mix_group1)
  message("  Unique mixtures in ", group2, ": ", n_mix_group2)
  mix_balance <- dcast(mix_tab, mixture ~ pam50, value.var = "N", fill = 0)
  if (nrow(mix_balance) > 1) {
    message("  Mixture balance (first 10 rows):")
    print(head(mix_balance, 10))
  }
} else {
  mix_balance <- NULL
  message("  Mixture: not available in design.")
}

# Warn if heavily confounded (e.g. one subtype only in one mixture)
n_mix_both <- NA_integer_
if (!is.null(mix_balance) && nrow(mix_balance) > 1) {
  g1_col <- names(mix_balance)[names(mix_balance) == group1]
  g2_col <- names(mix_balance)[names(mix_balance) == group2]
  if (length(g1_col) && length(g2_col)) {
    has_both <- (mix_balance[[g1_col]] > 0) & (mix_balance[[g2_col]] > 0)
    n_mix_both <- sum(has_both)
    if (n_mix_both < nrow(mix_balance) / 2)
      message("  WARNING: Contrast may be confounded by mixture — many mixtures have only one subtype. Consider including mixture in the model.")
  }
}

# -----------------------------------------------------------------------------
# STEP 1b — Group-wise coverage filter (before DA)
# -----------------------------------------------------------------------------
message("\n========== STEP 1b: Coverage filter (non-missing per group) ==========")
n_tot <- ncol(gm_mat)
n_grp1 <- max(1, sum(design_dt$pam50 == group1))
n_grp2 <- max(1, sum(design_dt$pam50 == group2))
is_grp1 <- design_dt$pam50 == group1
is_grp2 <- design_dt$pam50 == group2

keep_genes <- vapply(seq_len(n_genes), function(i) {
  x <- as.numeric(gm_mat[i, ])
  valid <- !is.na(x) & is.finite(x)
  pct_overall <- sum(valid) / n_tot
  pct_grp1   <- sum(valid[is_grp1]) / n_grp1
  pct_grp2   <- sum(valid[is_grp2]) / n_grp2
  (pct_overall >= pct_overall_min) && (pct_grp1 >= pct_group_min) && (pct_grp2 >= pct_group_min)
}, logical(1))

gm_mat <- gm_mat[keep_genes, , drop = FALSE]
n_genes <- nrow(gm_mat)
message("  Require non-missing in >= ", round(100 * pct_overall_min), "% overall and >= ", round(100 * pct_group_min), "% per subtype.")
message("  Genes retained after coverage filter: ", n_genes, " (dropped ", sum(!keep_genes), ")")

# =============================================================================
# STEP 2 — Build long format (MSstats-style, matrix-derived)
# =============================================================================
message("\n========== STEP 2: Build long format (matrix-derived, MSstats-style) ==========")
# This is a summarized gene matrix converted to long form for MSstats-like
# pipeline. We are NOT using raw feature-level peptide data; each protein is
# one summarized log2-abundance row. MSstats expects non-log intensity, so we
# pass 2^log2_abundance.

design_dt[, Run := matrix_col]
design_dt[, BioReplicate := matrix_col]
design_dt[, Condition := pam50]
setkey(design_dt, matrix_col)
cond_per_col <- setNames(design_dt$Condition, design_dt$matrix_col)

proteins <- rownames(gm_mat)
long_list <- lapply(seq_len(n_genes), function(i) {
  log2val <- as.numeric(gm_mat[i, ])
  data.table(
    ProteinName     = proteins[i],
    PeptideSequence = "PROTEIN_SUM",
    PrecursorCharge = 0L,
    FragmentIon     = "NA",
    ProductCharge   = 0L,
    IsotopeLabelType = "L",
    Condition       = cond_per_col[colnames(gm_mat)],
    BioReplicate    = colnames(gm_mat),
    Run             = colnames(gm_mat),
    Intensity       = 2^log2val
  )
})
msstats_long <- rbindlist(long_list)
# Drop rows with NA Intensity (missing in matrix)
msstats_long <- msstats_long[!is.na(Intensity) & is.finite(Intensity)]
message("  Long-format rows: ", nrow(msstats_long))

# =============================================================================
# STEP 3 — Run DA (MSstats then limma fallback)
# =============================================================================
message("\n========== STEP 3: Run subtype DA ==========")
method_used <- "limma"
res_dt <- NULL
msstats_ok <- FALSE

# Try MSstats
tryCatch({
  processed <- dataProcess(
    as.data.frame(msstats_long),
    normalization = FALSE,
    summaryMethod = "TMP",
    MBimpute = FALSE
  )
  contrast_mat <- matrix(0, nrow = 1, ncol = 2)
  rownames(contrast_mat) <- contrast_name
  colnames(contrast_mat) <- c(group1, group2)
  contrast_mat[1, group1] <- 1
  contrast_mat[1, group2] <- -1
  comparison <- groupComparison(contrast.matrix = contrast_mat, data = processed)
  res_ms <- as.data.frame(comparison$ComparisonResult)
  res_ms$ProteinName <- as.character(res_ms$ProteinName)
  if (!"log2FC" %in% names(res_ms) && "logFC" %in% names(res_ms))
    res_ms$log2FC <- res_ms$logFC
  if (!"adj.pvalue" %in% names(res_ms) && "Adjusted.Pvalue" %in% names(res_ms))
    res_ms$adj.pvalue <- res_ms$Adjusted.Pvalue
  res_dt <- as.data.table(res_ms)
  method_used <- "MSstats"
  msstats_ok <- TRUE
  message("  MSstats groupComparison succeeded.")
}, error = function(e) {
  message("  MSstats failed: ", conditionMessage(e))
  message("  Falling back to limma.")
})

# Fallback: limma on matrix
if (!msstats_ok) {
  condition <- factor(design_dt$pam50, levels = c(group2, group1))
  design_limma <- model.matrix(~ 0 + condition)
  colnames(design_limma) <- gsub("condition", "", colnames(design_limma))

  # Optional: include mixture as batch (use syntactically valid names for makeContrasts)
  use_mixture <- FALSE
  if ("mixture" %in% names(design_dt) && design_dt[, uniqueN(mixture) > 1] && design_dt[, all(!is.na(mixture))]) {
    mixture_f <- factor(design_dt$mixture)
    levels(mixture_f) <- make.names(levels(mixture_f))
    design_limma <- model.matrix(~ 0 + condition + mixture_f)
    colnames(design_limma) <- gsub("condition|mixture_f", "", colnames(design_limma))
    use_mixture <- TRUE
  }

  fit <- lmFit(gm_mat, design_limma)
  ctr <- makeContrasts(
    contrasts = paste0(group1, " - ", group2),
    levels = design_limma
  )
  fit2 <- contrasts.fit(fit, ctr)
  fit2 <- eBayes(fit2)
  res_limma <- topTable(fit2, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
  res_limma$GeneSymbol <- rownames(res_limma)
  res_limma$log2FC <- res_limma$logFC
  res_limma$adj.pvalue <- res_limma$adj.P.Val
  res_limma$pvalue <- res_limma$P.Value
  res_dt <- as.data.table(res_limma[, c("GeneSymbol", "log2FC", "pvalue", "adj.pvalue")])
  if (use_mixture) method_used <- "limma_with_mixture" else method_used <- "limma"
}

# Standardize output columns for consistent output
if ("ProteinName" %in% names(res_dt)) setnames(res_dt, "ProteinName", "GeneSymbol")
if (!"GeneSymbol" %in% names(res_dt)) res_dt[, GeneSymbol := res_dt[[names(res_dt)[1]]]]
if ("Adjusted.Pvalue" %in% names(res_dt) && !"adj.pvalue" %in% names(res_dt))
  setnames(res_dt, "Adjusted.Pvalue", "adj.pvalue")
if ("PValue" %in% names(res_dt) && !"pvalue" %in% names(res_dt))
  setnames(res_dt, "PValue", "pvalue")
res_dt[, contrast := contrast_name]
res_dt[, method_used := method_used]

# =============================================================================
# STEP 4 — Mixture reporting (already done above; ensure in summary)
# =============================================================================
mixture_in_model <- method_used %in% c("limma_with_mixture")

# =============================================================================
# STEP 5 — Save results and summary
# =============================================================================
message("\n========== STEP 5: Outputs ==========")

out_csv <- file.path(OUT_DIR, paste0("DA_MSstats_", contrast_name, ".csv"))
out_summary <- file.path(OUT_DIR, paste0("DA_MSstats_", contrast_name, "_summary.txt"))
fwrite(res_dt, out_csv)
message("  Wrote ", out_csv)

# Significance counts (FDR < 0.05, |log2FC| > 1)
sig <- res_dt[adj.pvalue < 0.05 & abs(log2FC) > 1]
n_sig <- nrow(sig)
# Positive log2FC = Luminal higher (Luminal - Basal); negative = Basal higher
group1_up <- res_dt[log2FC > 0][order(-log2FC)]
group2_up <- res_dt[log2FC < 0][order(log2FC)]

sink(out_summary)
cat("Subtype DA summary — ", contrast_name, "\n", sep = "")
cat("============================================\n\n")
cat("Group sizes:\n")
cat("  ", group1, ": ", sum(design_dt$pam50 == group1), "\n", sep = "")
cat("  ", group2, ": ", sum(design_dt$pam50 == group2), "\n\n", sep = "")
cat("Coverage filter: >= ", round(100 * pct_overall_min), "% non-missing overall, >= ", round(100 * pct_group_min), "% per subtype.\n", sep = "")
cat("Genes tested (after filter): ", n_genes, "\n", sep = "")
cat("Method used: ", method_used, "\n", sep = "")
cat("Mixture included in model: ", mixture_in_model, "\n\n", sep = "")
cat("Significant (FDR < 0.05, |log2FC| > 1): ", n_sig, "\n\n", sep = "")
cat("Top 10 ", group1, "-up genes:\n", sep = "")
print(group1_up[1:min(10, nrow(group1_up)), .(GeneSymbol, log2FC, adj.pvalue)])
cat("\nTop 10 ", group2, "-up genes:\n", sep = "")
print(group2_up[1:min(10, nrow(group2_up)), .(GeneSymbol, log2FC, adj.pvalue)])
if (!is.null(mix_balance)) {
  cat("\nMixture distribution:\n")
  print(mix_balance)
}
cat("\nInterpretation: positive log2FC = higher abundance in Luminal than Basal (Luminal - Basal).\n")
if (exists("n_miss") && n_miss > 0)
  cat("\nWARNING: ", n_miss, " design sample(s) were not found in gene matrix and were dropped.\n", sep = "")
if (!is.null(mix_balance) && nrow(mix_balance) > 1 && !is.na(n_mix_both) && n_mix_both < nrow(mix_balance) / 2)
  cat("\nWARNING: Contrast may be confounded by mixture (many mixtures have only one subtype).\n")
sink(NULL)
message("  Wrote ", out_summary)

# =============================================================================
# STEP 6 — Volcano plot
# =============================================================================
message("\n========== STEP 6: Volcano plot ==========")

# Optional breast subtype markers (only label if present)
markers_basal   <- c("EGFR", "KRT5", "KRT14", "KRT17")
markers_luminal <- c("ESR1", "GATA3", "KRT18", "FOXA1")
all_markers <- c(markers_basal, markers_luminal)

res_volcano <- as.data.frame(res_dt)
res_volcano$pval_plot <- res_volcano$adj.pvalue
res_volcano$pval_plot[res_volcano$pval_plot <= 0] <- min(res_volcano$pval_plot[res_volcano$pval_plot > 0], na.rm = TRUE)
y_volc <- -log10(res_volcano$pval_plot)
y_volc[is.infinite(y_volc)] <- max(y_volc[is.finite(y_volc)], na.rm = TRUE) + 0.5
sig_volc <- res_volcano$adj.pvalue < 0.05 & abs(res_volcano$log2FC) > 1
genes_volc <- res_volcano$GeneSymbol
top_n <- min(20, max(1, sum(sig_volc, na.rm = TRUE)))
ord <- order(res_volcano$adj.pvalue)
label_genes <- unique(genes_volc[ord[seq_len(top_n)]])
label_genes <- label_genes[!is.na(label_genes) & nzchar(label_genes)]

out_pdf <- file.path(OUT_DIR, paste0("volcano_MSstats_", contrast_name, ".pdf"))
pdf(out_pdf, width = 8, height = 7)
plot(res_volcano$log2FC, y_volc, pch = 20,
     col = ifelse(sig_volc, "red", "gray50"),
     xlab = "log2FC", ylab = "-log10(adj.pvalue)",
     main = paste0(contrast_name, " (", method_used, ")"))
abline(h = -log10(0.05), lty = 2, col = "gray40")
abline(v = c(-1, 1), lty = 2, col = "gray40")
if (length(label_genes) > 0) {
  idx <- match(label_genes, genes_volc)
  idx <- idx[!is.na(idx)]
  text(res_volcano$log2FC[idx], y_volc[idx], labels = genes_volc[idx], pos = 4, cex = 0.5)
}
# Highlight known markers if present
present_markers <- intersect(all_markers, genes_volc)
if (length(present_markers) > 0) {
  j <- match(present_markers, genes_volc)
  points(res_volcano$log2FC[j], y_volc[j], pch = 8, col = "blue", cex = 1.2)
  legend("topright", legend = c("FDR<0.05, |log2FC|>1", "Subtype markers"), col = c("red", "blue"), pch = c(20, 8), bty = "n")
} else {
  legend("topright", legend = "FDR<0.05, |log2FC|>1", col = "red", pch = 20, bty = "n")
}
dev.off()
message("  Wrote ", out_pdf)

# =============================================================================
# Final log
# =============================================================================
message("\n========== Done ==========")
message("Contrast: ", contrast_name, " (positive log2FC = Luminal vs Basal, Luminal higher)")
message("Method: ", method_used)
message("Significant (FDR<0.05, |log2FC|>1): ", n_sig)
message("Outputs: ", out_csv, ", ", out_pdf, ", ", out_summary)
