#!/usr/bin/env Rscript
# Diagnostics for PDC000120 Basal vs Luminal using limma on gene_matrix.csv
# (same spirit as data/scripts/DA_subtype_MSstats_PDC000120.R: mixture batch optional).
# Writes PDFs under data/results/PDC000120/diagnostics/
#
# Usage: cd data && Rscript --vanilla scripts/diagnostics_pdc_basal_luminal_limma.R

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", update = FALSE, ask = FALSE)
})
library(data.table)
library(limma)

DATA_DIR <- getwd()
if (!file.exists(file.path(DATA_DIR, "results", "PDC000120", "gene_matrix.csv")))
  DATA_DIR <- file.path(getwd(), "data")
if (!file.exists(file.path(DATA_DIR, "results", "PDC000120", "gene_matrix.csv")))
  stop("Run from project root or data/ with results/PDC000120/gene_matrix.csv present.")

RES <- file.path(DATA_DIR, "results", "PDC000120")
OUT <- file.path(RES, "diagnostics")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Contrast Luminal - Basal (positive coefficient => Luminal higher); matches DA_subtype_* scripts.
group1 <- "Luminal"
group2 <- "Basal"
contrast_name <- "Luminal_vs_Basal"
pct_overall_min <- 0.35
pct_group_min <- 0.25

MARKERS <- c("GATA3", "KRT18", "ESR1", "KRT5", "KRT17", "FOXC1")
PROFILE_GENES <- c("ESR1", "KRT18", "FOXC1")  # 3 panels; others used if missing

# --- Load design ---
design_dt <- fread(file.path(RES, "DA_subtype_tumor_only.csv"))
id_col <- names(design_dt)[grepl("matrix_sample_id|bioreplicate", names(design_dt), ignore.case = TRUE)][1]
pam50_col <- names(design_dt)[grepl("pam50", names(design_dt), ignore.case = TRUE)][1]
mix_col <- names(design_dt)[grepl("mixture", names(design_dt), ignore.case = TRUE)][1]
if (is.na(id_col)) id_col <- "matrix_sample_id"
if (is.na(pam50_col)) pam50_col <- "pam50"
design_dt[, (pam50_col) := trimws(as.character(get(pam50_col)))]
design_dt <- design_dt[get(pam50_col) != "" & !is.na(get(pam50_col))]
design_dt[get(pam50_col) %in% c("LumA", "LumB"), (pam50_col) := "Luminal"]
design_dt <- design_dt[get(pam50_col) %in% c("Luminal", "Basal")]
if (id_col != "matrix_sample_id") setnames(design_dt, id_col, "matrix_sample_id")
if (pam50_col != "pam50") setnames(design_dt, pam50_col, "pam50")
if (!is.na(mix_col) && mix_col != "mixture" && mix_col %in% names(design_dt))
  setnames(design_dt, mix_col, "mixture")
if (!"mixture" %in% names(design_dt)) design_dt[, mixture := NA_character_]

# --- Matrix ---
gm <- fread(file.path(RES, "gene_matrix.csv"))
gene_col <- names(gm)[1]
non_sample <- c(gene_col, "UniProtID", "uniprotid")
non_sample <- intersect(non_sample, names(gm))
cols <- setdiff(names(gm), non_sample)
design_lower <- tolower(trimws(design_dt$matrix_sample_id))
col_lower <- tolower(trimws(cols))
match_idx <- match(design_lower, col_lower)
design_dt <- design_dt[!is.na(match_idx)]
match_idx <- match_idx[!is.na(match_idx)]
mat_cols <- cols[match_idx]
gm_mat <- as.matrix(gm[, ..mat_cols])
storage.mode(gm_mat) <- "double"
rownames(gm_mat) <- gm[[gene_col]]
design_dt[, matrix_col := mat_cols]

n_tot <- ncol(gm_mat)
n_g1 <- max(1, sum(design_dt$pam50 == group1))
n_g2 <- max(1, sum(design_dt$pam50 == group2))
is_g1 <- design_dt$pam50 == group1
is_g2 <- design_dt$pam50 == group2

keep_genes <- vapply(seq_len(nrow(gm_mat)), function(i) {
  x <- gm_mat[i, ]
  valid <- !is.na(x) & is.finite(x)
  pct_overall <- sum(valid) / n_tot
  pct_1 <- sum(valid[is_g1]) / n_g1
  pct_2 <- sum(valid[is_g2]) / n_g2
  (pct_overall >= pct_overall_min) && (pct_1 >= pct_group_min) && (pct_2 >= pct_group_min)
}, logical(1))
gm_mat <- gm_mat[keep_genes, , drop = FALSE]

# --- limma with mixture batch (same as DA script when mixture available) ---
condition <- factor(design_dt$pam50, levels = c(group2, group1))
use_mixture <- FALSE
if ("mixture" %in% names(design_dt) && design_dt[, uniqueN(mixture) > 1] &&
    design_dt[, all(!is.na(mixture) & nzchar(as.character(mixture)))]) {
  mixture_f <- factor(design_dt$mixture)
  levels(mixture_f) <- make.names(levels(mixture_f))
  design_limma <- model.matrix(~ 0 + condition + mixture_f)
  colnames(design_limma) <- gsub("condition|mixture_f", "", colnames(design_limma))
  use_mixture <- TRUE
} else {
  design_limma <- model.matrix(~ 0 + condition)
  colnames(design_limma) <- gsub("condition", "", colnames(design_limma))
}

fit <- lmFit(gm_mat, design_limma)
ctr <- makeContrasts(contrasts = paste0(group1, " - ", group2), levels = design_limma)
fit <- contrasts.fit(fit, ctr)
fit <- eBayes(fit)

# --- 1) QQ plot of moderated t vs approximate reference t-distribution ---
tt <- fit$t[, 1]
n <- length(tt)
df_total <- fit$df.prior + fit$df.residual
df_total <- median(df_total, na.rm = TRUE)
df_total <- max(df_total, 2, na.rm = TRUE)
theo <- qt(ppoints(n), df = df_total)
obs <- sort(tt)

pdf(file.path(OUT, "qq_moderated_t_statistics.pdf"), width = 6, height = 6)
plot(theo, obs, xlab = sprintf("Expected t (df ~ %.1f)", df_total), ylab = "Observed moderated t",
     main = sprintf("QQ: moderated t (%s, limma)", contrast_name), pch = 16, col = grDevices::adjustcolor("black", 0.35))
abline(0, 1, col = "red", lty = 2)
dev.off()

# --- 2) Mean-variance (SA) plot — limma residual dispersion diagnostic ---
pdf(file.path(OUT, "mean_variance_SA_plot.pdf"), width = 6.5, height = 6)
plotSA(fit)
title(sub = sprintf("limma, mixture in model: %s", use_mixture))
dev.off()

# --- 3) Residual vs fitted for one marker (GATA3): per-sample linear model ---
gata3_row <- which(rownames(gm_mat) == "GATA3")[1]
if (is.na(gata3_row)) gata3_row <- 1L
y <- as.numeric(gm_mat[gata3_row, ])
X <- design_limma
ok <- is.finite(y)
fit_lm <- lm(y[ok] ~ 0 + X[ok, , drop = FALSE])
fv <- fitted(fit_lm)
rv <- residuals(fit_lm)
pdf(file.path(OUT, "residual_vs_fitted_GATA3.pdf"), width = 6, height = 6)
plot(fv, rv, xlab = "Fitted log2 abundance (GATA3)", ylab = "Residual",
     main = "Residual vs fitted — single protein (GATA3)", pch = 16,
     col = grDevices::adjustcolor("steelblue", 0.6))
abline(h = 0, lty = 2, col = "grey50")
dev.off()

# --- 4) Profile plots: log2 abundance by subtype for 3 genes ---
avail <- PROFILE_GENES[PROFILE_GENES %in% rownames(gm_mat)]
if (length(avail) < 3) {
  extra <- MARKERS[MARKERS %in% rownames(gm_mat)]
  avail <- unique(c(avail, extra))
}
avail <- head(avail, 3)
grp <- as.character(design_dt$pam50)
grp_f <- factor(grp, levels = c(group2, group1))

pdf(file.path(OUT, "profiles_three_markers_by_subtype.pdf"), width = 9, height = 3.5)
par(mfrow = c(1, length(avail)), mar = c(4, 4, 2, 1))
for (g in avail) {
  i <- which(rownames(gm_mat) == g)[1]
  yy <- as.numeric(gm_mat[i, ])
  boxplot(yy ~ grp_f, ylab = "log2 abundance", xlab = "PAM50 (tumor)", main = g, outline = FALSE)
  stripchart(yy ~ grp_f, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, col = grDevices::adjustcolor("steelblue", 0.55), cex = 0.85)
}
dev.off()

# Meta
sink(file.path(OUT, "run_meta.txt"))
cat("diagnostics_pdc_basal_luminal_limma.R\n")
cat("Contrast: ", contrast_name, "\n", sep = "")
cat("Samples: ", ncol(gm_mat), " (", group1, ": ", sum(grp == group1), ", ",
    group2, ": ", sum(grp == group2), ")\n", sep = "")
cat("Genes after coverage filter: ", nrow(gm_mat), "\n", sep = "")
cat("Mixture in design: ", use_mixture, "\n", sep = "")
cat("Profile genes plotted: ", paste(avail, collapse = ", "), "\n", sep = "")
sink(NULL)

message("Wrote diagnostics to ", OUT)
