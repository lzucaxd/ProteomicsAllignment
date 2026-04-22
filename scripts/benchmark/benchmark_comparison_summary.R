#!/usr/bin/env Rscript
# =============================================================================
# Benchmark Comparison Summary — Raw vs Celligner vs Bridge
# =============================================================================
# Runs additional metrics: pooled cross-domain DA, marker concordance, domain R²

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
OUTDIR <- file.path(REPO, "reports/benchmark_master/diagnostics/benchmark_comparison")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

load_gm <- function(path) {
  dt <- fread(path, header = TRUE)
  id_cols <- intersect(c("GeneSymbol", "UniProtID", "Gene"), names(dt))
  scols <- setdiff(names(dt), id_cols)
  mat <- as.matrix(dt[, ..scols])
  rownames(mat) <- dt[[names(dt)[1]]]
  mat
}

run_limma <- function(mat, groups) {
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

# ─── Load data ───────────────────────────────────────────────────────────
cat("Loading matrices...\n")
gm_breast <- load_gm(file.path(REPO, "data/results/PDC000120/gene_matrix.csv"))
gm_ccle <- load_gm(file.path(REPO, "data/results/CCLE_corrected/gene_matrix.csv"))

cell_dt <- fread(file.path(REPO, "reports/benchmark_master/celligner_all/celligner_aligned_matrix.csv"))
gm_cell <- t(as.matrix(cell_dt[, -1, with = FALSE]))
colnames(gm_cell) <- cell_dt[[1]]
rownames(gm_cell) <- names(cell_dt)[-1]

gm_bshift <- load_gm(file.path(REPO, "reports/benchmark_master/methods/bridge_aware/bridge_aware_shift_only_matrix.csv"))
gm_bscale <- load_gm(file.path(REPO, "reports/benchmark_master/methods/bridge_aware/bridge_aware_shift_scale_matrix.csv"))

# ─── Metadata ────────────────────────────────────────────────────────────
sm <- fread(file.path(REPO, "data/annotations/cptac/PDC000120/gene_matrix_subtype_mapping.csv"))
st_col <- if ("sample_type" %in% names(sm)) "sample_type" else "sample_type_if_available"

gm_cols <- colnames(gm_breast)
gm_cols_lower <- tolower(gm_cols)
samples_pam50 <- sm[tolower(get(st_col)) == "sample" &
                      tolower(pam50) %in% c("basal", "luma", "lumb") &
                      exists_in_gene_matrix == TRUE]
samples_pam50[, subtype := ifelse(tolower(pam50) == "basal", "Basal", "Luminal")]
samples_pam50[, matched_col := {
  idx <- match(tolower(matrix_sample_id), gm_cols_lower)
  fifelse(is.na(idx), NA_character_, gm_cols[idx])
}, by = seq_len(nrow(samples_pam50))]
tumors <- unique(samples_pam50[!is.na(matched_col)], by = "matched_col")

cptac_sub <- tumors[, .(sample_id = matched_col, condition = subtype, domain = "CPTAC")]

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
ccle_b_ids <- match_ccle(ccle_basal, colnames(gm_ccle))
ccle_l_ids <- match_ccle(ccle_luminal, colnames(gm_ccle))
ccle_sub <- data.table(
  sample_id = c(ccle_b_ids, ccle_l_ids),
  condition = c(rep("Basal", length(ccle_b_ids)), rep("Luminal", length(ccle_l_ids))),
  domain = "CCLE"
)
all_meta <- rbind(cptac_sub[, .(sample_id, condition, domain)], ccle_sub)

subtype_markers <- c("FOXA1","GATA3","KRT5","KRT14","KRT17","EGFR","ESR1","PGR","ERBB2","CDH1")

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: WITHIN-DOMAIN DA + CROSS-DOMAIN AGREEMENT (existing approach)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  ANALYSIS 1: Within-domain DA → cross-domain FC agreement\n")
cat(strrep("=", 70), "\n")

reps <- list(
  raw = NULL,
  celligner = gm_cell,
  bridge_shift = gm_bshift,
  bridge_scale = gm_bscale
)

shared_raw <- intersect(rownames(gm_breast), rownames(gm_ccle))
raw_mat <- cbind(gm_breast[shared_raw, ], gm_ccle[shared_raw, ])

results_within <- list()
for (rn in names(reps)) {
  mat <- if (rn == "raw") raw_mat else reps[[rn]]
  da_list <- list()
  for (dom in c("CPTAC", "CCLE")) {
    dm <- all_meta[toupper(domain) == dom]
    sids <- intersect(dm$sample_id, colnames(mat))
    if (length(sids) < 4) next
    dmat <- mat[, sids, drop = FALSE]
    dgrp <- dm[match(sids, sample_id), condition]
    keep <- rowSums(!is.na(dmat)) >= ncol(dmat) * 0.3
    dmat <- dmat[keep, , drop = FALSE]
    for (i in seq_len(nrow(dmat))) {
      nas <- is.na(dmat[i, ]); if (any(nas)) dmat[i, nas] <- median(dmat[i, !nas])
    }
    lvls <- sort(unique(dgrp))
    if (length(lvls) < 2) next
    da_list[[dom]] <- run_limma(dmat, dgrp)
  }
  if (length(da_list) == 2) {
    m <- merge(da_list$CPTAC[, .(gene, fc_cptac = logFC, p_cptac = adj.P.Val)],
               da_list$CCLE[, .(gene, fc_ccle = logFC, p_ccle = adj.P.Val)], by = "gene")
    m <- m[is.finite(fc_cptac) & is.finite(fc_ccle)]
    results_within[[rn]] <- data.table(
      representation = rn,
      n_genes = nrow(m),
      pearson_r = cor(m$fc_cptac, m$fc_ccle),
      spearman_rho = cor(m$fc_cptac, m$fc_ccle, method = "spearman"),
      dir_agree = mean(sign(m$fc_cptac) == sign(m$fc_ccle)),
      sig_both = sum(m$p_cptac < 0.05 & m$p_ccle < 0.05),
      sig_concordant = sum(m$p_cptac < 0.05 & m$p_ccle < 0.05 &
                             sign(m$fc_cptac) == sign(m$fc_ccle)),
      rmse = sqrt(mean((m$fc_cptac - m$fc_ccle)^2))
    )
    cat("  ", rn, ": r =", round(cor(m$fc_cptac, m$fc_ccle), 3),
        "dir_agree =", round(mean(sign(m$fc_cptac) == sign(m$fc_ccle))*100, 1), "%\n")
  }
}
within_dt <- rbindlist(results_within)

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: POOLED CROSS-DOMAIN DA (where bridge shines)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  ANALYSIS 2: Pooled cross-domain DA (subtype, adjusting for domain)\n")
cat(strrep("=", 70), "\n")

results_pooled <- list()
for (rn in names(reps)) {
  mat <- if (rn == "raw") raw_mat else reps[[rn]]
  sids <- intersect(all_meta$sample_id, colnames(mat))
  if (length(sids) < 4) next
  dmat <- mat[, sids, drop = FALSE]
  m_meta <- all_meta[match(sids, sample_id)]

  keep <- rowSums(!is.na(dmat)) >= ncol(dmat) * 0.3
  dmat <- dmat[keep, , drop = FALSE]
  for (i in seq_len(nrow(dmat))) {
    nas <- is.na(dmat[i, ]); if (any(nas)) dmat[i, nas] <- median(dmat[i, !nas])
  }

  grp <- factor(m_meta$condition)
  dom <- factor(m_meta$domain)
  design <- model.matrix(~ 0 + grp + dom)
  colnames(design) <- gsub("^grp|^dom", "", colnames(design))

  fit <- lmFit(dmat, design)
  cm <- makeContrasts(Luminal - Basal, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  setDT(tt)

  n_sig <- sum(tt$adj.P.Val < 0.05)
  n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
  n_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)

  # Marker check
  mk <- tt[gene %in% subtype_markers]
  n_markers_present <- nrow(mk)
  n_markers_sig <- sum(mk$adj.P.Val < 0.05)

  results_pooled[[rn]] <- data.table(
    representation = rn,
    n_samples = length(sids),
    n_genes = nrow(dmat),
    n_sig_005 = n_sig,
    n_up = n_up,
    n_down = n_down,
    markers_present = n_markers_present,
    markers_sig = n_markers_sig,
    median_abs_fc = median(abs(tt$logFC))
  )
  cat("  ", rn, ": ", length(sids), "samples,", n_sig, "sig genes,",
      n_markers_sig, "/", n_markers_present, "markers sig\n")

  # Save marker table
  fwrite(mk[order(adj.P.Val)], file.path(OUTDIR, paste0("pooled_markers_", rn, ".csv")))
  fwrite(tt[order(adj.P.Val)], file.path(OUTDIR, paste0("pooled_da_", rn, ".csv")))
}
pooled_dt <- rbindlist(results_pooled)

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: DOMAIN R² (how much does PC1 capture domain vs biology)
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  ANALYSIS 3: Domain R² and condition R² from PCA\n")
cat(strrep("=", 70), "\n")

results_pca <- list()
for (rn in names(reps)) {
  mat <- if (rn == "raw") raw_mat else reps[[rn]]
  sids <- intersect(all_meta$sample_id, colnames(mat))
  if (length(sids) < 4) next
  dmat <- mat[, sids, drop = FALSE]
  m_meta <- all_meta[match(sids, sample_id)]

  keep <- rowSums(!is.na(dmat)) >= ncol(dmat) * 0.5
  dmat <- dmat[keep, , drop = FALSE]
  dmat[is.na(dmat)] <- 0

  gv <- apply(dmat, 1, var)
  top <- order(gv, decreasing = TRUE)[1:min(2000, nrow(dmat))]
  pca <- prcomp(t(dmat[top, ]), center = TRUE, scale. = FALSE)
  pve <- summary(pca)$importance[2, 1:5]

  pc1 <- pca$x[, 1]
  pc2 <- pca$x[, 2]
  domain_r2_pc1 <- summary(lm(pc1 ~ factor(m_meta$domain)))$r.squared
  domain_r2_pc2 <- summary(lm(pc2 ~ factor(m_meta$domain)))$r.squared
  cond_r2_pc1 <- summary(lm(pc1 ~ factor(m_meta$condition)))$r.squared
  cond_r2_pc2 <- summary(lm(pc2 ~ factor(m_meta$condition)))$r.squared

  results_pca[[rn]] <- data.table(
    representation = rn,
    domain_r2_pc1 = round(domain_r2_pc1, 4),
    domain_r2_pc2 = round(domain_r2_pc2, 4),
    condition_r2_pc1 = round(cond_r2_pc1, 4),
    condition_r2_pc2 = round(cond_r2_pc2, 4),
    pve_pc1 = round(pve[1] * 100, 1),
    pve_pc2 = round(pve[2] * 100, 1)
  )
  cat("  ", rn, ": domain_R2_PC1=", round(domain_r2_pc1, 3),
      "cond_R2_PC1=", round(cond_r2_pc1, 3),
      "PVE1=", round(pve[1]*100,1), "%\n")
}
pca_dt <- rbindlist(results_pca)

# ═══════════════════════════════════════════════════════════════════════════
# ANALYSIS 4: MARKER GENE DIRECTION + SIGNIFICANCE TABLE
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  ANALYSIS 4: Marker gene concordance (per-domain DA)\n")
cat(strrep("=", 70), "\n")

expected <- c(FOXA1 = 1, GATA3 = 1, KRT5 = -1, KRT14 = -1, KRT17 = -1,
              EGFR = -1, ESR1 = 1, PGR = 1, ERBB2 = 1, CDH1 = 1)

marker_rows <- list()
for (rn in names(reps)) {
  mat <- if (rn == "raw") raw_mat else reps[[rn]]
  for (dom in c("CPTAC", "CCLE")) {
    dm <- all_meta[toupper(domain) == dom]
    sids <- intersect(dm$sample_id, colnames(mat))
    if (length(sids) < 4) next
    dmat <- mat[, sids, drop = FALSE]
    dgrp <- dm[match(sids, sample_id), condition]
    keep <- rowSums(!is.na(dmat)) >= ncol(dmat) * 0.3
    dmat <- dmat[keep, , drop = FALSE]
    for (i in seq_len(nrow(dmat))) {
      nas <- is.na(dmat[i, ]); if (any(nas)) dmat[i, nas] <- median(dmat[i, !nas])
    }
    lvls <- sort(unique(dgrp))
    if (length(lvls) < 2) next
    da <- run_limma(dmat, dgrp)
    for (g in names(expected)) {
      if (!g %in% da$gene) {
        marker_rows[[length(marker_rows) + 1]] <- data.table(
          representation = rn, domain = dom, gene = g,
          logFC = NA, adj_p = NA, expected_sign = expected[g],
          correct_dir = NA, significant = NA, present = FALSE
        )
      } else {
        row <- da[gene == g]
        marker_rows[[length(marker_rows) + 1]] <- data.table(
          representation = rn, domain = dom, gene = g,
          logFC = round(row$logFC, 4), adj_p = row$adj.P.Val,
          expected_sign = expected[g],
          correct_dir = sign(row$logFC) == expected[g],
          significant = row$adj.P.Val < 0.05,
          present = TRUE
        )
      }
    }
  }
}
marker_dt <- rbindlist(marker_rows)

# Summarize by representation
marker_summary <- marker_dt[, .(
  markers_present = sum(present),
  correct_direction = sum(correct_dir, na.rm = TRUE),
  sig_correct = sum(correct_dir & significant, na.rm = TRUE),
  sig_wrong = sum(!correct_dir & significant, na.rm = TRUE)
), by = .(representation, domain)]

cat("\n  Marker summary:\n")
print(marker_summary)

# ═══════════════════════════════════════════════════════════════════════════
# FINAL COMPARISON TABLE
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  FINAL COMPARISON TABLE\n")
cat(strrep("=", 70), "\n\n")

final <- merge(within_dt, pca_dt, by = "representation", all = TRUE)
final <- merge(final, pooled_dt[, .(representation, pooled_n_sig = n_sig_005,
                                     pooled_markers_sig = markers_sig)],
               by = "representation", all = TRUE)

cat("═══ CROSS-DOMAIN FC AGREEMENT (subtype, within-domain DA) ═══\n")
print(within_dt[order(-pearson_r)])

cat("\n═══ PCA STRUCTURE ═══\n")
print(pca_dt)

cat("\n═══ POOLED DA (domain-adjusted, both domains combined) ═══\n")
print(pooled_dt[order(-n_sig_005)])

cat("\n═══ MARKER GENE CONCORDANCE (per domain) ═══\n")
print(marker_summary[order(representation, domain)])

# Save all
fwrite(final, file.path(OUTDIR, "full_comparison.csv"))
fwrite(within_dt, file.path(OUTDIR, "within_domain_agreement.csv"))
fwrite(pca_dt, file.path(OUTDIR, "pca_structure.csv"))
fwrite(pooled_dt, file.path(OUTDIR, "pooled_da_summary.csv"))
fwrite(marker_dt, file.path(OUTDIR, "marker_detail.csv"))
fwrite(marker_summary, file.path(OUTDIR, "marker_summary.csv"))

cat("\n  All outputs saved to:", OUTDIR, "\n")
cat("  Done.\n")
