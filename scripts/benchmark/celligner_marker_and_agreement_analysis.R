#!/usr/bin/env Rscript
# =============================================================================
# Celligner (all-data): Marker Gene Profiles + Up/Down Agreement Analysis
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
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
OUTDIR <- file.path(REPO, "reports/benchmark_master/diagnostics/celligner_all_markers")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ─── Load Celligner matrix ───────────────────────────────────────────────
cat("Loading Celligner all-data matrix...\n")
cell_dt <- fread(file.path(REPO, "reports/benchmark_master/celligner_all/celligner_aligned_matrix.csv"))
cell_samples <- cell_dt[[1]]
cell_genes <- names(cell_dt)[-1]
gm_cell <- t(as.matrix(cell_dt[, -1, with = FALSE]))
colnames(gm_cell) <- cell_samples
rownames(gm_cell) <- cell_genes
cat("  ", nrow(gm_cell), "genes ×", ncol(gm_cell), "samples\n")

# ─── Load metadata ──────────────────────────────────────────────────────
cell_meta <- fread(file.path(REPO, "reports/benchmark_master/celligner_all/sample_metadata.csv"))

# ─── Load subtype and tissue subsets ────────────────────────────────────
sm <- fread(file.path(REPO, "data/results/PDC000120/gene_matrix_subtype_mapping.csv"))
st_col <- if ("sample_type" %in% names(sm)) "sample_type" else "sample_type_if_available"
gm_breast_cols <- colnames(fread(file.path(REPO, "data/results/PDC000120/gene_matrix.csv"), nrows = 0))
gm_breast_cols <- setdiff(gm_breast_cols, c("GeneSymbol", "UniProtID"))

samples_pam50 <- sm[tolower(get(st_col)) == "sample" &
                      tolower(pam50) %in% c("basal", "luma", "lumb") &
                      exists_in_gene_matrix == TRUE]
samples_pam50[, subtype := ifelse(tolower(pam50) == "basal", "Basal", "Luminal")]
gm_cols_lower <- tolower(gm_breast_cols)
samples_pam50[, matched_col := {
  idx <- match(tolower(matrix_sample_id), gm_cols_lower)
  fifelse(is.na(idx), NA_character_, gm_breast_cols[idx])
}, by = seq_len(nrow(samples_pam50))]
tumors <- unique(samples_pam50[!is.na(matched_col)], by = "matched_col")

cptac_subtype <- tumors[, .(sample_id = matched_col, condition = subtype, domain = "CPTAC")]

ccle_basal <- c("HCC70", "HCC1806", "HCC1143", "MDA-MB-468")
ccle_luminal <- c("CAMA-1", "MCF7", "T-47D", "ZR-75-1")
match_ccle <- function(lines, cols) {
  out <- character(0)
  for (ln in lines) {
    pat <- gsub("-", ".", ln, fixed = TRUE)
    m <- grep(pat, cols, ignore.case = TRUE, value = TRUE)
    if (length(m) >= 1) out <- c(out, m[1])
  }
  out
}
ccle_b_ids <- match_ccle(ccle_basal, colnames(gm_cell))
ccle_l_ids <- match_ccle(ccle_luminal, colnames(gm_cell))
ccle_subtype <- data.table(
  sample_id = c(ccle_b_ids, ccle_l_ids),
  condition = c(rep("Basal", length(ccle_b_ids)), rep("Luminal", length(ccle_l_ids))),
  domain = "CCLE"
)

# Breast vs lung
ccle_info_df <- read.csv(file.path(REPO, "data/ccle_peptide/sample_info_ccle.csv"),
                          stringsAsFactors = FALSE, fill = TRUE)
ccle_info <- as.data.table(ccle_info_df)
ccle_info <- ccle_info[nchar(Cell.Line) > 0]
setnames(ccle_info, "Cell.Line", "Cell Line", skip_absent = TRUE)
setnames(ccle_info, "Tissue.of.Origin", "Tissue of Origin", skip_absent = TRUE)

tissue_map <- setNames(ccle_info[["Tissue of Origin"]], ccle_info[["Cell Line"]])
ccle_cols <- colnames(gm_cell)[colnames(gm_cell) %in% cell_meta[domain == "CCLE", sample_id]]

normalize_name <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
ccle_col_tissue <- rep(NA_character_, length(ccle_cols))
for (i in seq_along(ccle_cols)) {
  col <- ccle_cols[i]
  if (col %in% names(tissue_map)) { ccle_col_tissue[i] <- tissue_map[[col]]; next }
  col_norm <- normalize_name(col)
  for (nm in names(tissue_map)) {
    if (normalize_name(nm) == col_norm) { ccle_col_tissue[i] <- tissue_map[[nm]]; break }
  }
}
names(ccle_col_tissue) <- ccle_cols
ccle_breast_ids <- ccle_cols[!is.na(ccle_col_tissue) & tolower(ccle_col_tissue) == "breast"]
ccle_lung_ids <- ccle_cols[!is.na(ccle_col_tissue) & tolower(ccle_col_tissue) == "lung"]

cptac_breast_ids <- intersect(gm_breast_cols, colnames(gm_cell))
cptac_lung_ids <- intersect(
  colnames(fread(file.path(REPO, "data/results/PDC000153/gene_matrix.csv"), nrows = 0)),
  colnames(gm_cell)
)
cptac_lung_ids <- setdiff(cptac_lung_ids, c("GeneSymbol", "UniProtID", "prot"))

all_meta_subtype <- rbind(cptac_subtype[, .(sample_id, condition, domain)], ccle_subtype)
all_meta_bvl <- rbind(
  data.table(sample_id = cptac_breast_ids, condition = "Breast", domain = "CPTAC"),
  data.table(sample_id = cptac_lung_ids, condition = "Lung", domain = "CPTAC"),
  data.table(sample_id = ccle_breast_ids, condition = "Breast", domain = "CCLE"),
  data.table(sample_id = ccle_lung_ids, condition = "Lung", domain = "CCLE")
)

# ─── Marker panels ──────────────────────────────────────────────────────
subtype_markers <- c("FOXA1", "GATA3", "KRT5", "KRT14", "KRT17", "EGFR",
                      "ESR1", "PGR", "ERBB2", "CDH1")
# Expected direction: positive logFC = up in Luminal
subtype_expected <- c(FOXA1 = "up_Luminal", GATA3 = "up_Luminal",
                       KRT5 = "up_Basal", KRT14 = "up_Basal", KRT17 = "up_Basal",
                       EGFR = "up_Basal", ESR1 = "up_Luminal", PGR = "up_Luminal",
                       ERBB2 = "up_Luminal", CDH1 = "up_Luminal")

bvl_markers <- c("NKX2-1", "SFTPB", "SFTPC", "NAPSA", "GATA3", "FOXA1",
                  "ESR1", "KRT19", "EGFR", "ERBB2", "CDH1", "VIM")
bvl_expected <- c("NKX2-1" = "up_Lung", SFTPB = "up_Lung", SFTPC = "up_Lung",
                   NAPSA = "up_Lung", GATA3 = "up_Breast", FOXA1 = "up_Breast",
                   ESR1 = "up_Breast", KRT19 = "up_Breast",
                   EGFR = "up_Lung", ERBB2 = "up_Breast", CDH1 = "up_Breast",
                   VIM = "up_Lung")

# ═══════════════════════════════════════════════════════════════════════════
# PART 1: MARKER GENE PROFILE PLOTS
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  MARKER GENE PROFILE PLOTS\n")
cat(strrep("=", 70), "\n")

plot_marker <- function(mat, meta, gene, task_label, outpath) {
  if (!gene %in% rownames(mat)) { cat("  SKIP:", gene, "(not in matrix)\n"); return() }
  sids <- intersect(meta$sample_id, colnames(mat))
  vals <- mat[gene, sids]
  m <- meta[match(sids, sample_id)]
  grp <- paste(m$domain, m$condition, sep = "\n")
  grps <- sort(unique(grp))
  gvals <- lapply(grps, function(g) vals[grp == g])
  ns <- sapply(gvals, length)
  labels <- paste0(grps, "\n(n=", ns, ")")

  pal <- c("#2196F3", "#E91E63", "#64B5F6", "#F48FB1",
           "#4CAF50", "#FF9800", "#81C784", "#FFB74D")[seq_along(grps)]

  png(outpath, width = 700, height = 500, res = 120)
  par(mar = c(8, 5, 3, 1))
  boxplot(gvals, names = labels, las = 2, main = paste0(gene, " — ", task_label),
          ylab = "Celligner Abundance", col = pal, outline = FALSE, cex.axis = 0.75)
  stripchart(gvals, vertical = TRUE, method = "jitter", pch = 16,
             cex = 0.4, col = "gray30", add = TRUE)
  dev.off()
  cat("  ", gene, "\n")
}

# Subtype markers
cat("\n── Breast Subtype Markers ──\n")
for (g in subtype_markers) {
  plot_marker(gm_cell, all_meta_subtype, g, "Celligner Breast Subtype",
              file.path(OUTDIR, paste0(g, "_subtype_celligner.png")))
}

# BvL markers
cat("\n── Breast vs Lung Markers ──\n")
for (g in bvl_markers) {
  plot_marker(gm_cell, all_meta_bvl, g, "Celligner Breast vs Lung",
              file.path(OUTDIR, paste0(g, "_bvl_celligner.png")))
}

# ═══════════════════════════════════════════════════════════════════════════
# PART 2: LIMMA DA + DETAILED UP/DOWN AGREEMENT
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  UP/DOWN AGREEMENT ANALYSIS\n")
cat(strrep("=", 70), "\n")

run_limma <- function(mat, groups, contrast_name = "B_vs_A") {
  grp <- factor(groups)
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  fit <- lmFit(mat, design)
  lvls <- levels(grp)
  cstr <- paste0(lvls[2], "-", lvls[1])
  cm <- makeContrasts(contrasts = cstr, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  setDT(tt)
  tt
}

run_da_for_task <- function(mat, meta, task_name, markers, expected_dirs) {
  cat("\n── ", task_name, " ──\n")
  results <- list()

  for (dom in c("CPTAC", "CCLE")) {
    dm <- meta[toupper(domain) == dom]
    sids <- intersect(dm$sample_id, colnames(mat))
    if (length(sids) < 4) { cat("  skip ", dom, "\n"); next }
    dmat <- mat[, sids, drop = FALSE]
    dgrp <- dm[match(sids, sample_id), condition]

    keep <- rowSums(!is.na(dmat)) >= ncol(dmat) * 0.3
    dmat <- dmat[keep, , drop = FALSE]
    for (i in seq_len(nrow(dmat))) {
      nas <- is.na(dmat[i, ]); if (any(nas)) dmat[i, nas] <- median(dmat[i, !nas])
    }

    lvls <- sort(unique(dgrp[!is.na(dgrp)]))
    if (length(lvls) < 2) { cat("  skip ", dom, ": <2 groups\n"); next }

    cname <- paste0(lvls[2], "_vs_", lvls[1])
    cat("  ", dom, ":", length(sids), "samples →", cname, "\n")
    da <- run_limma(dmat, dgrp, cname)
    results[[dom]] <- da
  }

  if (is.null(results$CPTAC) || is.null(results$CCLE)) {
    cat("  Cannot compute agreement (missing domain)\n")
    return(NULL)
  }

  # Merge
  m <- merge(results$CPTAC[, .(gene, fc_cptac = logFC, p_cptac = adj.P.Val)],
             results$CCLE[, .(gene, fc_ccle = logFC, p_ccle = adj.P.Val)],
             by = "gene")
  m <- m[is.finite(fc_cptac) & is.finite(fc_ccle)]

  # Classify genes
  fdr_thresh <- 0.05
  fc_thresh <- 0.0  # any direction

  m[, sig_cptac := p_cptac < fdr_thresh]
  m[, sig_ccle := p_ccle < fdr_thresh]
  m[, dir_cptac := ifelse(fc_cptac > 0, "up", "down")]
  m[, dir_ccle := ifelse(fc_ccle > 0, "up", "down")]
  m[, dir_agree := dir_cptac == dir_ccle]

  cat("\n  === Agreement Summary (", task_name, ") ===\n")
  cat("  Total shared genes:", nrow(m), "\n")
  cat("  Direction agreement (all genes):", round(mean(m$dir_agree) * 100, 1), "%\n")
  cat("  Direction agreement (sig in CPTAC):",
      round(mean(m[sig_cptac == TRUE, dir_agree]) * 100, 1), "%",
      "(n=", sum(m$sig_cptac), ")\n")
  cat("  Direction agreement (sig in both):",
      round(mean(m[sig_cptac == TRUE & sig_ccle == TRUE, dir_agree]) * 100, 1), "%",
      "(n=", sum(m$sig_cptac & m$sig_ccle), ")\n")

  # Breakdown by direction
  cat("\n  --- Directional Breakdown ---\n")
  cat("  CPTAC up-regulated (logFC > 0):", sum(m$fc_cptac > 0), "\n")
  cat("    → CCLE agrees (also up):", sum(m$fc_cptac > 0 & m$fc_ccle > 0), "\n")
  cat("    → CCLE disagrees (down):", sum(m$fc_cptac > 0 & m$fc_ccle <= 0), "\n")
  cat("  CPTAC down-regulated (logFC < 0):", sum(m$fc_cptac < 0), "\n")
  cat("    → CCLE agrees (also down):", sum(m$fc_cptac < 0 & m$fc_ccle < 0), "\n")
  cat("    → CCLE disagrees (up):", sum(m$fc_cptac < 0 & m$fc_ccle >= 0), "\n")

  # Significant genes
  sig_both <- m[sig_cptac == TRUE & sig_ccle == TRUE]
  if (nrow(sig_both) > 0) {
    cat("\n  --- Significant in Both Domains (FDR < 0.05) ---\n")
    cat("  Total:", nrow(sig_both), "\n")
    cat("  Same direction:", sum(sig_both$dir_agree), "(",
        round(mean(sig_both$dir_agree) * 100, 1), "%)\n")
    cat("  Opposite direction:", sum(!sig_both$dir_agree), "(",
        round(mean(!sig_both$dir_agree) * 100, 1), "%)\n")

    # Top concordant
    concordant <- sig_both[dir_agree == TRUE][order(-abs(fc_cptac))][1:min(15, .N)]
    if (nrow(concordant) > 0) {
      cat("\n  Top 15 concordant significant genes:\n")
      cat(sprintf("  %-12s  CPTAC_FC=%6.3f  CCLE_FC=%6.3f\n",
                  concordant$gene, concordant$fc_cptac, concordant$fc_ccle), sep = "")
    }

    # Top discordant
    discordant <- sig_both[dir_agree == FALSE][order(-abs(fc_cptac))][1:min(15, .N)]
    if (nrow(discordant) > 0) {
      cat("\n  Top discordant significant genes:\n")
      cat(sprintf("  %-12s  CPTAC_FC=%6.3f  CCLE_FC=%6.3f\n",
                  discordant$gene, discordant$fc_cptac, discordant$fc_ccle), sep = "")
    }
  }

  # Marker direction check
  cat("\n  --- Marker Gene Check ---\n")
  for (g in names(expected_dirs)) {
    if (!g %in% m$gene) { cat("  ", g, ": NOT IN DATA\n"); next }
    row <- m[gene == g]
    exp <- expected_dirs[g]
    cptac_ok <- (grepl("up", exp) & row$fc_cptac > 0) | (grepl("down", exp) & row$fc_cptac < 0)
    # For markers: "up_Luminal" means logFC > 0 (Luminal vs Basal), "up_Basal" means logFC < 0
    if (grepl("up_Basal|up_Lung", exp)) {
      cptac_ok <- row$fc_cptac < 0  # second level vs first; if Basal expected up, logFC should be negative (Luminal-Basal)
      ccle_ok <- row$fc_ccle < 0
    } else {
      cptac_ok <- row$fc_cptac > 0
      ccle_ok <- row$fc_ccle > 0
    }
    status_c <- ifelse(cptac_ok, "✓", "✗")
    status_e <- ifelse(ccle_ok, "✓", "✗")
    cat(sprintf("  %-10s  expected: %-12s  CPTAC: %s (FC=%6.3f, p=%.1e)  CCLE: %s (FC=%6.3f, p=%.1e)\n",
                g, exp, status_c, row$fc_cptac, row$p_cptac, status_e, row$fc_ccle, row$p_ccle))
  }

  # FC agreement plot
  r_val <- cor(m$fc_cptac, m$fc_ccle, method = "pearson")
  rho_val <- cor(m$fc_cptac, m$fc_ccle, method = "spearman")

  png(file.path(OUTDIR, paste0("fc_agreement_", gsub(" ", "_", task_name), ".png")),
      width = 900, height = 800, res = 130)
  par(mar = c(5, 5, 4, 1))
  plot(m$fc_cptac, m$fc_ccle, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.15),
       xlab = "CPTAC logFC", ylab = "CCLE logFC",
       main = paste0("Celligner — ", task_name))
  abline(0, 1, col = "red", lwd = 1.5, lty = 2)
  abline(h = 0, v = 0, col = "gray70")

  # Color significant concordant/discordant
  sig_conc <- m[sig_cptac & sig_ccle & dir_agree == TRUE]
  sig_disc <- m[sig_cptac & sig_ccle & dir_agree == FALSE]
  points(sig_conc$fc_cptac, sig_conc$fc_ccle, pch = 16, cex = 0.5, col = rgb(0, 0.6, 0, 0.5))
  points(sig_disc$fc_cptac, sig_disc$fc_ccle, pch = 16, cex = 0.5, col = rgb(0.8, 0, 0, 0.5))

  # Highlight markers
  for (g in intersect(names(expected_dirs), m$gene)) {
    row <- m[gene == g]
    points(row$fc_cptac, row$fc_ccle, pch = 18, cex = 1.5, col = "blue")
    text(row$fc_cptac, row$fc_ccle, g, pos = 4, cex = 0.55, col = "blue")
  }

  legend("topleft", legend = c(
    paste0("r = ", round(r_val, 3)),
    paste0("ρ = ", round(rho_val, 3)),
    paste0("dir agree = ", round(mean(m$dir_agree) * 100, 1), "%"),
    paste0("n = ", nrow(m)),
    paste0("sig concordant = ", nrow(sig_conc)),
    paste0("sig discordant = ", nrow(sig_disc))
  ), bty = "n", cex = 0.75,
  text.col = c("black", "black", "black", "black", "darkgreen", "red"))
  dev.off()
  cat("  Saved FC agreement plot\n")

  # Concordance by FC magnitude bins
  cat("\n  --- Concordance by FC Magnitude ---\n")
  bins <- c(0, 0.25, 0.5, 1.0, 2.0, Inf)
  for (i in seq_len(length(bins) - 1)) {
    lo <- bins[i]; hi <- bins[i + 1]
    sub <- m[abs(fc_cptac) >= lo & abs(fc_cptac) < hi]
    if (nrow(sub) > 0) {
      cat(sprintf("  |FC_CPTAC| [%.2f, %.2f): n=%d, dir_agree=%.1f%%\n",
                  lo, hi, nrow(sub), mean(sub$dir_agree) * 100))
    }
  }

  # Save full table
  fwrite(m, file.path(OUTDIR, paste0("agreement_table_", gsub(" ", "_", task_name), ".csv")))

  m
}

# Run both tasks
ag_subtype <- run_da_for_task(gm_cell, all_meta_subtype, "Breast Subtype",
                               subtype_markers, subtype_expected)
ag_bvl <- run_da_for_task(gm_cell, all_meta_bvl, "Breast vs Lung",
                           bvl_markers, bvl_expected)

cat("\n", strrep("=", 70), "\n")
cat("  DONE — outputs in:", OUTDIR, "\n")
cat(strrep("=", 70), "\n")
