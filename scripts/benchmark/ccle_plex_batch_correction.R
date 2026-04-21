#!/usr/bin/env Rscript
# =============================================================================
# CCLE Plex Batch Correction Diagnostic
# =============================================================================
# 1. Quantify how much variance Protein 10-Plex ID explains
# 2. Try: removeBatchEffect, ComBat, linear residualization
# 3. Re-evaluate: PCA, UMAP, tissue R², breast vs lung structure
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(sva)
})

.local_args <- commandArgs(trailingOnly = FALSE)
.local_file <- .local_args[startsWith(.local_args, "--file=")]
if (length(.local_file)) {
  .local_bench <- dirname(normalizePath(sub("^--file=", "", .local_file[1L])))
} else {
  .local_bench <- normalizePath(file.path(getwd(), "scripts", "benchmark"), mustWork = FALSE)
}
source(file.path(.local_bench, "harmonize_paths.R"))
REPO <- harmonize_repo_root()
CCLE_MATRIX <- file.path(REPO, "data/results/CCLE_corrected/gene_matrix.csv")
CCLE_SAMPLE <- file.path(REPO, "data/ccle_peptide/sample_info_ccle.csv")
OUTDIR <- file.path(REPO, "reports/benchmark_master/diagnostics/ccle_plex_correction")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ─── Load data ───────────────────────────────────────────────────────────
cat("Loading data...\n")
dt <- fread(CCLE_MATRIX, header = TRUE)
gene_col <- names(dt)[1]
id_cols <- intersect(c("GeneSymbol", "UniProtID", "Gene"), names(dt))
scols <- setdiff(names(dt), id_cols)
mat <- as.matrix(dt[, ..scols])
rownames(mat) <- dt[[gene_col]]
cat("  Raw matrix:", nrow(mat), "genes ×", ncol(mat), "samples\n")

# ─── Load sample info and map columns to plex ───────────────────────────
info <- read.csv(CCLE_SAMPLE, stringsAsFactors = FALSE, fill = TRUE)
info <- info[nchar(info$Cell.Line) > 0, ]

# Build cell line -> plex mapping
plex_map <- setNames(as.character(info$Protein.10.Plex.ID), info$Cell.Line)
tissue_map <- setNames(info$Tissue.of.Origin, info$Cell.Line)

# Also build CCLE code mapping
code_map_plex <- setNames(as.character(info$Protein.10.Plex.ID), info$CCLE.Code)
code_map_tissue <- setNames(info$Tissue.of.Origin, info$CCLE.Code)

# Match matrix columns to plex and tissue
sample_ids <- colnames(mat)
sample_plex <- rep(NA_character_, length(sample_ids))
sample_tissue <- rep(NA_character_, length(sample_ids))

normalize_name <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

info_norm <- data.frame(
  cell_line = info$Cell.Line,
  ccle_code = info$CCLE.Code,
  tissue = info$Tissue.of.Origin,
  plex = as.character(info$Protein.10.Plex.ID),
  norm_cl = normalize_name(info$Cell.Line),
  norm_cc = normalize_name(info$CCLE.Code),
  stringsAsFactors = FALSE
)
# Remove bridge plex (plex 0) entries to avoid double-matching
info_non_bridge <- info_norm[info_norm$plex != "0", ]

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  sid_norm <- normalize_name(sid)

  # Try exact match on cell line name first
  idx <- match(sid, info_non_bridge$cell_line)
  if (!is.na(idx)) {
    sample_plex[i] <- info_non_bridge$plex[idx]
    sample_tissue[i] <- info_non_bridge$tissue[idx]
    next
  }
  # Try normalized match on cell line
  idx <- match(sid_norm, info_non_bridge$norm_cl)
  if (!is.na(idx)) {
    sample_plex[i] <- info_non_bridge$plex[idx]
    sample_tissue[i] <- info_non_bridge$tissue[idx]
    next
  }
  # Try normalized match on CCLE code
  idx <- match(sid_norm, info_non_bridge$norm_cc)
  if (!is.na(idx)) {
    sample_plex[i] <- info_non_bridge$plex[idx]
    sample_tissue[i] <- info_non_bridge$tissue[idx]
    next
  }
}

cat("  Plex mapped:", sum(!is.na(sample_plex)), "of", length(sample_ids), "samples\n")
cat("  Tissue mapped:", sum(!is.na(sample_tissue)), "of", length(sample_ids), "samples\n")
cat("  Unique plexes:", length(unique(sample_plex[!is.na(sample_plex)])), "\n")

# Collapse tissues
tissue_remap <- c(
  "Breast" = "Breast", "Lung" = "Lung", "Ovary" = "Ovary",
  "Large Intestine" = "Colorectal",
  "Haematopoietic and Lymphoid Tissue" = "Blood/Lymphoid",
  "Acute Myeloid Leukemia" = "Blood/Lymphoid", "Lymphoma" = "Blood/Lymphoid",
  "Skin" = "Skin", "Central Nervous System" = "CNS",
  "Pancreas" = "Pancreas", "Stomach" = "Upper GI", "Oesophagus" = "Upper GI",
  "Endometrium" = "Endometrium", "Liver" = "Liver", "Kidney" = "Kidney",
  "Urinary Tract" = "Urinary Tract", "Upper Aerodigestive Tract" = "Head & Neck",
  "Prostate" = "Prostate", "Soft Tissue" = "Other", "Bone" = "Other",
  "Thyroid" = "Other", "Pleura" = "Other", "Autonomic Ganglia" = "Other",
  "Biliary Tract" = "Other"
)
sample_tissue_col <- ifelse(sample_tissue %in% names(tissue_remap),
                             tissue_remap[sample_tissue], "Other")

# ─── Filter genes: 70% prevalence + SD ──────────────────────────────────
prev <- rowMeans(!is.na(mat))
sd_vals <- apply(mat, 1, sd, na.rm = TRUE)
keep <- (prev >= 0.70) & (sd_vals >= 0.01)
mat_f <- mat[keep, ]
cat("  Filtered:", nrow(mat_f), "genes (70% prevalence + SD >= 0.01)\n")

# Impute remaining NAs with gene median
for (i in seq_len(nrow(mat_f))) {
  nas <- is.na(mat_f[i, ])
  if (any(nas)) mat_f[i, nas] <- median(mat_f[i, !nas], na.rm = TRUE)
}

# Keep only samples with plex + tissue info
has_info <- !is.na(sample_plex) & !is.na(sample_tissue)
mat_f <- mat_f[, has_info]
plex <- sample_plex[has_info]
tissue <- sample_tissue_col[has_info]
cat("  Final: ", ncol(mat_f), "samples with plex + tissue info\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 1: QUANTIFY PLEX VARIANCE
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  STEP 1: PLEX VARIANCE QUANTIFICATION\n")
cat(strrep("=", 70), "\n\n")

pca_base <- prcomp(t(mat_f), center = TRUE, scale. = FALSE, rank. = 10)
pve_base <- summary(pca_base)$importance[2, 1:10]

# R² for plex and tissue on each PC
r2_plex <- numeric(10)
r2_tissue <- numeric(10)
for (j in 1:10) {
  r2_plex[j] <- summary(lm(pca_base$x[, j] ~ factor(plex)))$r.squared
  r2_tissue[j] <- summary(lm(pca_base$x[, j] ~ factor(tissue)))$r.squared
}

cat("  PC  |  PVE%  | Plex R² | Tissue R²\n")
cat("  ----|--------|---------|----------\n")
for (j in 1:10) {
  cat(sprintf("  PC%d | %5.1f%% |  %5.3f  |   %5.3f\n",
              j, pve_base[j] * 100, r2_plex[j], r2_tissue[j]))
}

# Weighted summary
cat(sprintf("\n  Weighted Plex R²  (PC1-5): %.4f\n",
            sum(r2_plex[1:5] * pve_base[1:5]) / sum(pve_base[1:5])))
cat(sprintf("  Weighted Tissue R² (PC1-5): %.4f\n",
            sum(r2_tissue[1:5] * pve_base[1:5]) / sum(pve_base[1:5])))

# Joint model: tissue + plex together
r2_joint <- numeric(10)
for (j in 1:10) {
  r2_joint[j] <- summary(lm(pca_base$x[, j] ~ factor(plex) + factor(tissue)))$r.squared
}
cat(sprintf("  Weighted Joint R² (PC1-5): %.4f\n",
            sum(r2_joint[1:5] * pve_base[1:5]) / sum(pve_base[1:5])))

# ═══════════════════════════════════════════════════════════════════════════
# STEP 2: BATCH CORRECTIONS
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  STEP 2: APPLYING BATCH CORRECTIONS\n")
cat(strrep("=", 70), "\n")

corrections <- list()

# ── 2a: removeBatchEffect (limma) ────────────────────────────────────────
cat("\n  [2a] limma::removeBatchEffect...\n")
mat_rbe <- removeBatchEffect(mat_f, batch = factor(plex))
corrections[["removeBatchEffect"]] <- mat_rbe
cat("    Done.\n")

# ── 2b: ComBat (sva) ────────────────────────────────────────────────────
cat("\n  [2b] sva::ComBat...\n")
# ComBat needs at least 2 samples per batch — check
plex_ct <- table(plex)
ok_plexes <- names(plex_ct)[plex_ct >= 2]
combat_keep <- plex %in% ok_plexes
cat("    Plexes with >= 2 samples:", length(ok_plexes), "\n")
cat("    Samples used:", sum(combat_keep), "\n")

# ComBat with parametric adjustment, preserving tissue (biological variable)
mod_tissue <- model.matrix(~ factor(tissue[combat_keep]))
mat_combat <- ComBat(dat = mat_f[, combat_keep],
                     batch = factor(plex[combat_keep]),
                     mod = mod_tissue,
                     par.prior = TRUE, prior.plots = FALSE)
corrections[["ComBat"]] <- mat_combat
cat("    Done.\n")

# ── 2c: Linear regression residualization ────────────────────────────────
cat("\n  [2c] Linear regression residualization on plex...\n")
plex_design <- model.matrix(~ 0 + factor(plex))
mat_lm <- mat_f
fit_lm <- lmFit(mat_f, plex_design)
mat_lm <- mat_f - fitted(fit_lm)
# Add back the global gene mean so values are interpretable
mat_lm <- mat_lm + rowMeans(mat_f)
corrections[["lm_residual"]] <- mat_lm
cat("    Done.\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 3: EVALUATE EACH CORRECTION
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  STEP 3: EVALUATION\n")
cat(strrep("=", 70), "\n")

TISSUE_COLORS <- c(
  Breast = "#E91E63", Lung = "#4CAF50", Ovary = "#00BCD4",
  Colorectal = "#FF9800", "Blood/Lymphoid" = "#9C27B0",
  Skin = "#795548", CNS = "#3F51B5", Pancreas = "#CDDC39",
  "Upper GI" = "#FF5722", Endometrium = "#AB47BC",
  Liver = "#8D6E63", Kidney = "#5C6BC0", "Urinary Tract" = "#26A69A",
  "Head & Neck" = "#78909C", Prostate = "#AED581", Other = "#BDBDBD"
)

eval_and_plot <- function(mat_eval, tissue_eval, plex_eval, label, outdir) {
  n_samp <- ncol(mat_eval)
  n_gene <- nrow(mat_eval)

  # PCA
  pca <- prcomp(t(mat_eval), center = TRUE, scale. = FALSE, rank. = 10)
  pve <- summary(pca)$importance[2, 1:10]

  r2p <- numeric(5)
  r2t <- numeric(5)
  for (j in 1:5) {
    r2p[j] <- summary(lm(pca$x[, j] ~ factor(plex_eval)))$r.squared
    r2t[j] <- summary(lm(pca$x[, j] ~ factor(tissue_eval)))$r.squared
  }

  wr2_plex <- sum(r2p * pve[1:5]) / sum(pve[1:5])
  wr2_tissue <- sum(r2t * pve[1:5]) / sum(pve[1:5])

  cat(sprintf("\n  ── %s (%d samples × %d genes) ──\n", label, n_samp, n_gene))
  cat(sprintf("    Weighted Plex R²  (PC1-5): %.4f\n", wr2_plex))
  cat(sprintf("    Weighted Tissue R² (PC1-5): %.4f\n", wr2_tissue))
  for (j in 1:5) {
    cat(sprintf("    PC%d: PVE=%.1f%%, Plex R²=%.3f, Tissue R²=%.3f\n",
                j, pve[j]*100, r2p[j], r2t[j]))
  }

  # PCA plot colored by tissue
  tis_fac <- factor(tissue_eval)
  tis_levs <- levels(tis_fac)
  cols <- TISSUE_COLORS[tis_levs]
  cols[is.na(cols)] <- "#999999"

  png(file.path(outdir, paste0("pca_tissue_", gsub("[^a-zA-Z0-9]", "_", label), ".png")),
      width = 1100, height = 800, res = 130)
  par(mar = c(5, 5, 4, 10), xpd = TRUE)
  plot(pca$x[, 1], pca$x[, 2],
       col = cols[as.integer(tis_fac)], pch = 16, cex = 1.0,
       xlab = sprintf("PC1 (%.1f%%)", pve[1]*100),
       ylab = sprintf("PC2 (%.1f%%)", pve[2]*100),
       main = sprintf("%s — Tissue R²=%.3f, Plex R²=%.3f", label, wr2_tissue, wr2_plex))
  legend("topright", inset = c(-0.25, 0), legend = tis_levs,
         col = cols, pch = 16, cex = 0.6, title = "Tissue")
  dev.off()

  # PCA plot colored by plex
  plex_fac <- factor(plex_eval)
  n_plex <- nlevels(plex_fac)
  plex_cols <- rainbow(n_plex, s = 0.7, v = 0.8)

  png(file.path(outdir, paste0("pca_plex_", gsub("[^a-zA-Z0-9]", "_", label), ".png")),
      width = 1100, height = 800, res = 130)
  par(mar = c(5, 5, 4, 10), xpd = TRUE)
  plot(pca$x[, 1], pca$x[, 2],
       col = plex_cols[as.integer(plex_fac)], pch = 16, cex = 1.0,
       xlab = sprintf("PC1 (%.1f%%)", pve[1]*100),
       ylab = sprintf("PC2 (%.1f%%)", pve[2]*100),
       main = sprintf("%s — colored by Plex", label))
  dev.off()

  # Breast vs Lung focused PCA
  is_bl <- tissue_eval %in% c("Breast", "Lung")
  if (sum(is_bl) >= 10) {
    mat_bl <- mat_eval[, is_bl]
    tis_bl <- tissue_eval[is_bl]
    gvar <- apply(mat_bl, 1, var, na.rm = TRUE)
    top_g <- order(gvar, decreasing = TRUE)[1:min(2000, nrow(mat_bl))]
    pca_bl <- prcomp(t(mat_bl[top_g, ]), center = TRUE, scale. = FALSE)
    pve_bl <- summary(pca_bl)$importance[2, 1:2]
    r2_bl_pc1 <- summary(lm(pca_bl$x[, 1] ~ factor(tis_bl)))$r.squared

    blcols <- ifelse(tis_bl == "Breast", "#E91E63", "#4CAF50")
    png(file.path(outdir, paste0("pca_breast_vs_lung_", gsub("[^a-zA-Z0-9]", "_", label), ".png")),
        width = 900, height = 700, res = 130)
    par(mar = c(5, 5, 4, 8), xpd = TRUE)
    plot(pca_bl$x[, 1], pca_bl$x[, 2],
         col = blcols, pch = 16, cex = 1.2,
         xlab = sprintf("PC1 (%.1f%%)", pve_bl[1]*100),
         ylab = sprintf("PC2 (%.1f%%)", pve_bl[2]*100),
         main = sprintf("%s — Breast(%d) vs Lung(%d), PC1 R²=%.3f",
                        label, sum(tis_bl == "Breast"), sum(tis_bl == "Lung"), r2_bl_pc1))
    legend("topright", inset = c(-0.18, 0),
           legend = c(paste0("Breast (", sum(tis_bl == "Breast"), ")"),
                      paste0("Lung (", sum(tis_bl == "Lung"), ")")),
           col = c("#E91E63", "#4CAF50"), pch = 16, cex = 0.8)
    dev.off()

    cat(sprintf("    Breast vs Lung PC1 R²: %.4f (%d Breast, %d Lung)\n",
                r2_bl_pc1, sum(tis_bl == "Breast"), sum(tis_bl == "Lung")))
  }

  list(wr2_plex = wr2_plex, wr2_tissue = wr2_tissue,
       plex_r2_pcs = r2p, tissue_r2_pcs = r2t, pve = pve[1:5])
}

# ── Evaluate uncorrected ──
res_uncorrected <- eval_and_plot(mat_f, tissue, plex, "Uncorrected", OUTDIR)

# ── Evaluate each correction ──
results <- list(uncorrected = res_uncorrected)

for (nm in names(corrections)) {
  m <- corrections[[nm]]
  # Match tissue/plex for ComBat which may have dropped samples
  if (nm == "ComBat") {
    results[[nm]] <- eval_and_plot(m, tissue[combat_keep], plex[combat_keep], nm, OUTDIR)
  } else {
    results[[nm]] <- eval_and_plot(m, tissue, plex, nm, OUTDIR)
  }
}

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4: SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  SUMMARY TABLE\n")
cat(strrep("=", 70), "\n\n")

summary_dt <- rbindlist(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.table(
    method = nm,
    weighted_plex_r2 = round(r$wr2_plex, 4),
    weighted_tissue_r2 = round(r$wr2_tissue, 4),
    pc1_plex_r2 = round(r$plex_r2_pcs[1], 4),
    pc1_tissue_r2 = round(r$tissue_r2_pcs[1], 4),
    pc1_pve = round(r$pve[1] * 100, 1)
  )
}))
fwrite(summary_dt, file.path(OUTDIR, "plex_correction_summary.csv"))
print(summary_dt)

# ═══════════════════════════════════════════════════════════════════════════
# STEP 5: Save corrected matrices for downstream use
# ═══════════════════════════════════════════════════════════════════════════
cat("\nSaving corrected matrices...\n")
for (nm in names(corrections)) {
  m <- corrections[[nm]]
  out_dt <- data.table(GeneSymbol = rownames(m), as.data.table(m))
  fwrite(out_dt, file.path(OUTDIR, paste0("ccle_corrected_", nm, ".csv")))
  cat("  ", nm, ":", nrow(m), "×", ncol(m), "\n")
}

cat("\n", strrep("=", 70), "\n")
cat("  COMPLETE — outputs in:", OUTDIR, "\n")
cat(strrep("=", 70), "\n")
