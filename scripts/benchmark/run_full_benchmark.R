#!/usr/bin/env Rscript
# =============================================================================
# Full Benchmark — Raw + Celligner + Bridge × (subtype + breast_vs_lung) × (CPTAC+CCLE)
# =============================================================================
# Sources harmonize_paths.R + subset_strategies.R (union CCLE subtype panel).
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
source(file.path(.local_bench, "subset_strategies.R"))
REPO <- harmonize_repo_root()

CPTAC_BREAST  <- file.path(REPO, "data/results/PDC000120/gene_matrix.csv")
CPTAC_LUNG    <- file.path(REPO, "data/results/PDC000153/gene_matrix.csv")
CPTAC_OVARIAN <- file.path(REPO, "data/results/PDC000127/gene_matrix.csv")
CPTAC_UTERINE <- file.path(REPO, "data/results/PDC000204/gene_matrix.csv")
CCLE_MATRIX   <- file.path(REPO, "data/results/CCLE_corrected/gene_matrix.csv")
SUBTYPE_MAP   <- file.path(REPO, "data/annotations/cptac/PDC000120/gene_matrix_subtype_mapping.csv")
CCLE_SAMPLE   <- file.path(REPO, "data/ccle_peptide/sample_info_ccle.csv")
CELLIGNER_MAT <- file.path(REPO, "reports/benchmark_master/celligner_all/celligner_aligned_matrix.csv")
CELLIGNER_META<- file.path(REPO, "reports/benchmark_master/celligner_all/sample_metadata.csv")
BRIDGE_SHIFT  <- file.path(REPO, "reports/benchmark_master/methods/bridge_aware/bridge_aware_shift_only_matrix.csv")
BRIDGE_SCALE  <- file.path(REPO, "reports/benchmark_master/methods/bridge_aware/bridge_aware_shift_scale_matrix.csv")
OUTDIR        <- file.path(REPO, "reports/benchmark_master")

# ─── HELPER: limma DA ─────────────────────────────────────────────────────
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
  setcolorder(tt, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
  tt[, contrast := contrast_name]
  tt
}

# ─── HELPER: FC agreement ─────────────────────────────────────────────────
fc_agreement <- function(da1, da2) {
  m <- merge(da1[, .(gene, fc1 = logFC)], da2[, .(gene, fc2 = logFC)], by = "gene")
  m <- m[is.finite(fc1) & is.finite(fc2)]
  if (nrow(m) < 10) return(data.table(n = nrow(m)))
  data.table(
    n = nrow(m),
    pearson_r = cor(m$fc1, m$fc2, method = "pearson"),
    spearman_rho = cor(m$fc1, m$fc2, method = "spearman"),
    direction_agree = mean(sign(m$fc1) == sign(m$fc2)),
    median_abs_diff = median(abs(m$fc1 - m$fc2)),
    rmse = sqrt(mean((m$fc1 - m$fc2)^2))
  )
}

# ─── HELPER: PCA domain R² ────────────────────────────────────────────────
domain_r2 <- function(mat, domain_labels) {
  ok <- complete.cases(t(mat))
  if (sum(ok) < 10) return(NA_real_)
  pca <- prcomp(t(mat[, ok, drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 2)
  summary(lm(pca$x[, 1] ~ factor(domain_labels[ok])))$r.squared
}

# ─── HELPER: PCA plot ──────────────────────────────────────────────────────
pca_plot <- function(mat, meta, color_col, shape_col = NULL, title, outpath) {
  ok <- complete.cases(t(mat))
  mat_ok <- mat[, ok, drop = FALSE]
  meta_ok <- meta[ok, ]
  gv <- apply(mat_ok, 1, var, na.rm = TRUE)
  top <- order(gv, decreasing = TRUE)[1:min(2000, nrow(mat_ok))]
  m <- mat_ok[top, ]
  m[is.na(m)] <- 0
  pca <- prcomp(t(m), center = TRUE, scale. = FALSE)
  pve <- round(summary(pca)$importance[2, 1:2] * 100, 1)
  colfac <- factor(meta_ok[[color_col]])
  palette <- c("#E91E63", "#4CAF50", "#2196F3", "#FF9800", "#9C27B0",
               "#795548", "#00BCD4", "#607D8B")
  cols <- palette[as.integer(colfac)]
  pchs <- if (!is.null(shape_col)) c(16, 17)[as.integer(factor(meta_ok[[shape_col]]))] else 16
  png(outpath, width = 1000, height = 750, res = 130)
  par(mar = c(5, 5, 3, 9), xpd = TRUE)
  plot(pca$x[, 1], pca$x[, 2], col = cols, pch = pchs,
       xlab = paste0("PC1 (", pve[1], "%)"), ylab = paste0("PC2 (", pve[2], "%)"),
       main = title, cex = 1.1)
  legend("topright", inset = c(-0.22, 0), legend = levels(colfac),
         col = palette[seq_along(levels(colfac))], pch = 16, cex = 0.7, title = color_col)
  if (!is.null(shape_col)) {
    sfac <- levels(factor(meta_ok[[shape_col]]))
    legend("bottomright", inset = c(-0.22, 0), legend = sfac,
           pch = c(16, 17)[seq_along(sfac)], cex = 0.7, title = shape_col)
  }
  dev.off()
}

# ─── HELPER: marker profile plot ──────────────────────────────────────────
marker_profile <- function(mat, meta, gene, outpath) {
  if (!gene %in% rownames(mat)) return(invisible(NULL))
  vals <- mat[gene, ]
  grp <- paste(meta$domain, meta$condition, sep = "\n")
  grps <- sort(unique(grp))
  gvals <- lapply(grps, function(g) vals[grp == g])
  pal <- c("#2196F3", "#E91E63", "#64B5F6", "#F48FB1")[seq_along(grps)]
  png(outpath, width = 600, height = 450, res = 120)
  par(mar = c(7, 5, 3, 1))
  boxplot(gvals, names = grps, las = 2, main = gene, ylab = "Abundance",
          col = pal, outline = FALSE)
  stripchart(gvals, vertical = TRUE, method = "jitter", pch = 16,
             cex = 0.5, col = "gray30", add = TRUE)
  dev.off()
}

# ─── HELPER: load gene matrix ─────────────────────────────────────────────
load_gm <- function(path) {
  dt <- fread(path, header = TRUE)
  gene_col <- names(dt)[1]
  id_cols <- intersect(c("GeneSymbol", "UniProtID", "Gene"), names(dt))
  scols <- setdiff(names(dt), id_cols)
  mat <- as.matrix(dt[, ..scols])
  rownames(mat) <- dt[[gene_col]]
  mat
}

# ─── HELPER: run one task on one representation ───────────────────────────
run_task <- function(mat, meta, task_name, rep_name, markers, outbase) {
  outdir <- file.path(outbase, rep_name, task_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  results <- list()
  for (dom in c("CPTAC", "CCLE")) {
    dm <- meta[toupper(domain) == dom]
    sids <- intersect(dm$sample_id, colnames(mat))
    if (length(sids) < 4) { message("  skip ", dom, " (", length(sids), " samples)"); next }

    dmat <- mat[, sids, drop = FALSE]
    dgrp <- dm[match(sids, sample_id), condition]

    # Filter genes: >=30% non-NA
    keep <- rowSums(!is.na(dmat)) >= ncol(dmat) * 0.3
    dmat <- dmat[keep, , drop = FALSE]
    # Impute remaining NAs with row median for limma
    for (i in seq_len(nrow(dmat))) {
      nas <- is.na(dmat[i, ])
      if (any(nas)) dmat[i, nas] <- median(dmat[i, !nas], na.rm = TRUE)
    }

    lvls <- sort(unique(dgrp[!is.na(dgrp)]))
    if (length(lvls) < 2) {
      message("  skip ", dom, ": only group(s) ", paste(lvls, collapse = ","), " present")
      next
    }
    cname <- paste0(lvls[2], "_vs_", lvls[1], "_", dom)
    message("  ", dom, ": ", length(sids), " samples × ", nrow(dmat), " genes → ", cname)

    da <- run_limma(dmat, dgrp, cname)
    da[, inference_type := "representation_level_limma"]
    da[, representation := rep_name]
    da[, domain := dom]

    ddir <- file.path(outdir, tolower(dom))
    dir.create(ddir, showWarnings = FALSE)
    fwrite(da, file.path(ddir, "da_limma_result.csv"))

    mk <- da[gene %in% markers]
    fwrite(mk, file.path(ddir, "marker_summary.csv"))

    note <- data.table(domain = dom, task = task_name, representation = rep_name,
                       n_group1 = sum(dgrp == lvls[1]), n_group2 = sum(dgrp == lvls[2]),
                       n_genes = nrow(dmat), inference_type = "representation_level_limma")
    fwrite(note, file.path(ddir, "sample_counts.csv"))

    results[[dom]] <- da
  }

  # Cross-domain agreement
  if (!is.null(results[["CPTAC"]]) && !is.null(results[["CCLE"]])) {
    ag <- fc_agreement(results$CPTAC, results$CCLE)
    ag[, representation := rep_name]
    ag[, task := task_name]
    fwrite(ag, file.path(outdir, "cross_domain_agreement.csv"))
    fwrite(merge(results$CPTAC[, .(gene, fc_cptac = logFC)],
                 results$CCLE[, .(gene, fc_ccle = logFC)], by = "gene"),
           file.path(outdir, "fc_scatter_data.csv"))
    results$agreement <- ag
  }

  results
}

# ==========================================================================
cat("\n", strrep("=", 70), "\n")
cat("  FULL BENCHMARK RUN\n")
cat(strrep("=", 70), "\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 1: Load all data
# ═══════════════════════════════════════════════════════════════════════════
cat("STEP 1: Loading data...\n")
gm_breast  <- load_gm(CPTAC_BREAST)
gm_lung    <- load_gm(CPTAC_LUNG)
gm_ovarian <- load_gm(CPTAC_OVARIAN)
gm_uterine <- load_gm(CPTAC_UTERINE)
gm_ccle    <- load_gm(CCLE_MATRIX)

# Celligner all-data matrix (samples×genes); transpose to genes×samples
cell_dt <- fread(CELLIGNER_MAT, header = TRUE)
cell_samples <- cell_dt[[1]]
cell_genes <- names(cell_dt)[-1]
gm_cell <- t(as.matrix(cell_dt[, -1, with = FALSE]))
colnames(gm_cell) <- cell_samples
rownames(gm_cell) <- cell_genes

# Celligner metadata (domain, primary_site per sample)
cell_meta <- fread(CELLIGNER_META)

# Bridge-aware matrices (genes × samples, first col = GeneSymbol)
gm_bridge_shift <- load_gm(BRIDGE_SHIFT)
gm_bridge_scale <- load_gm(BRIDGE_SCALE)

sm <- fread(SUBTYPE_MAP)
ccle_info_df <- read.csv(CCLE_SAMPLE, stringsAsFactors = FALSE, fill = TRUE)
ccle_info <- as.data.table(ccle_info_df)
ccle_info <- ccle_info[nchar(Cell.Line) > 0]
setnames(ccle_info, "Cell.Line", "Cell Line", skip_absent = TRUE)
setnames(ccle_info, "Tissue.of.Origin", "Tissue of Origin", skip_absent = TRUE)
setnames(ccle_info, "CCLE.Code", "CCLE Code", skip_absent = TRUE)

cat("  CPTAC breast: ", ncol(gm_breast), "samples ×", nrow(gm_breast), "genes\n")
cat("  CPTAC lung:   ", ncol(gm_lung), "samples ×", nrow(gm_lung), "genes\n")
cat("  CPTAC ovarian:", ncol(gm_ovarian), "samples ×", nrow(gm_ovarian), "genes\n")
cat("  CPTAC uterine:", ncol(gm_uterine), "samples ×", nrow(gm_uterine), "genes\n")
cat("  CCLE:         ", ncol(gm_ccle), "samples ×", nrow(gm_ccle), "genes\n")
cat("  Celligner:    ", ncol(gm_cell), "samples ×", nrow(gm_cell), "genes\n")
cat("  Bridge shift: ", ncol(gm_bridge_shift), "samples ×", nrow(gm_bridge_shift), "genes\n")
cat("  Bridge scale: ", ncol(gm_bridge_scale), "samples ×", nrow(gm_bridge_scale), "genes\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 2: Build subsets
# ═══════════════════════════════════════════════════════════════════════════
cat("\nSTEP 2: Building subsets...\n")

# ── Subtype (Task A) ──────────────────────────────────────────────────────
st_col <- if ("sample_type" %in% names(sm)) "sample_type" else "sample_type_if_available"
# "sample" = non-pool (i.e. tumor or NAT); filter to those with PAM50
samples_with_pam50 <- sm[tolower(get(st_col)) == "sample" &
                          tolower(pam50) %in% c("basal", "luma", "lumb") &
                          exists_in_gene_matrix == TRUE]
samples_with_pam50[, subtype := ifelse(tolower(pam50) == "basal", "Basal", "Luminal")]

# Case-insensitive match of matrix_sample_id to gene matrix columns
gm_cols <- colnames(gm_breast)
gm_cols_lower <- tolower(gm_cols)
samples_with_pam50[, matched_col := {
  idx <- match(tolower(matrix_sample_id), gm_cols_lower)
  fifelse(is.na(idx), NA_character_, gm_cols[idx])
}, by = seq_len(nrow(samples_with_pam50))]
# Deduplicate to one row per matched column
tumors <- unique(samples_with_pam50[!is.na(matched_col)], by = "matched_col")

mix_ct <- tumors[, .(nB = sum(subtype == "Basal"), nL = sum(subtype == "Luminal"),
                     tot = .N), by = mixture]
mix_ct[, keep := (nB >= 1) & (nL >= 1) & !(pmin(nB, nL) == 1 & tot >= 6)]
kept_mix <- mix_ct[keep == TRUE, mixture]
cptac_subtype <- tumors[mixture %in% kept_mix,
                        .(sample_id = matched_col, condition = subtype,
                          domain = "CPTAC", mixture)]
cat("  CPTAC subtype: ", cptac_subtype[condition == "Basal", .N], "Basal,",
    cptac_subtype[condition == "Luminal", .N], "Luminal\n")

# CCLE: same Basal/Luminal panel as data/processed/union/sample_meta_breast_subtype.csv
# (sample_id mapped to colnames(gm_ccle); lines absent from matrix are dropped).
UNION_SUBTYPE_META <- file.path(REPO, "data/processed/union/sample_meta_breast_subtype.csv")
if (!file.exists(UNION_SUBTYPE_META))
  UNION_SUBTYPE_META <- file.path(REPO, "data/processed/sample_meta_breast_subtype.csv")
if (!file.exists(UNION_SUBTYPE_META))
  stop("Missing sample_meta_breast_subtype.csv for CCLE union panel.")
ccle_subtype <- build_subtype_subset_ccle(
  CCLE_MATRIX,
  union_meta_path = UNION_SUBTYPE_META
)[, .(sample_id, condition = subtype, domain)]
cat("  CCLE subtype (union meta): ", ccle_subtype[condition == "Basal", .N], "Basal,",
    ccle_subtype[condition == "Luminal", .N], "Luminal\n")

# ── Breast vs Lung (Task B) ──────────────────────────────────────────────
breast_ids <- colnames(gm_breast)
lung_ids <- colnames(gm_lung)
cptac_bvl <- data.table(
  sample_id = c(breast_ids, lung_ids),
  condition = c(rep("Breast", length(breast_ids)), rep("Lung", length(lung_ids))),
  domain = "CPTAC"
)

# Build tissue lookup: exact match Cell Line name → matrix column
tissue_map <- setNames(ccle_info[["Tissue of Origin"]], ccle_info[["Cell Line"]])
ccle_cols <- colnames(gm_ccle)
ccle_col_tissue <- rep(NA_character_, length(ccle_cols))
names(ccle_col_tissue) <- ccle_cols

for (i in seq_along(ccle_cols)) {
  col <- ccle_cols[i]
  if (col %in% names(tissue_map)) { ccle_col_tissue[i] <- tissue_map[[col]]; next }
  # Normalize: remove hyphens, underscores, uppercase
  col_norm <- toupper(gsub("[^A-Za-z0-9]", "", col))
  for (nm in names(tissue_map)) {
    nm_norm <- toupper(gsub("[^A-Za-z0-9]", "", nm))
    if (col_norm == nm_norm) { ccle_col_tissue[i] <- tissue_map[[nm]]; break }
  }
}
cat("  CCLE tissue mapped:", sum(!is.na(ccle_col_tissue)), "of", length(ccle_cols), "\n")
ccle_breast_ids <- ccle_cols[!is.na(ccle_col_tissue) & tolower(ccle_col_tissue) == "breast"]
ccle_lung_ids   <- ccle_cols[!is.na(ccle_col_tissue) & tolower(ccle_col_tissue) == "lung"]
ccle_bvl <- data.table(
  sample_id = c(ccle_breast_ids, ccle_lung_ids),
  condition = c(rep("Breast", length(ccle_breast_ids)), rep("Lung", length(ccle_lung_ids))),
  domain = "CCLE"
)
cat("  CPTAC breast vs lung:", length(breast_ids), "vs", length(lung_ids), "\n")
cat("  CCLE breast vs lung:", length(ccle_breast_ids), "vs", length(ccle_lung_ids), "\n")

# Save subset tables
sub_dir <- file.path(OUTDIR, "tasks")
dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(rbind(cptac_subtype[, .(sample_id, condition, domain)], ccle_subtype),
       file.path(sub_dir, "subtype_subset_all.csv"))
fwrite(rbind(cptac_bvl, ccle_bvl),
       file.path(sub_dir, "breast_vs_lung_subset_all.csv"))

# ═══════════════════════════════════════════════════════════════════════════
# STEP 3: Build representation matrices
# ═══════════════════════════════════════════════════════════════════════════
cat("\nSTEP 3: Building representation matrices...\n")

# ── RAW: inner join of CPTAC breast + CCLE ──
shared_bc <- intersect(rownames(gm_breast), rownames(gm_ccle))
raw_subtype_mat <- cbind(gm_breast[shared_bc, ], gm_ccle[shared_bc, ])
cat("  Raw (subtype): ", ncol(raw_subtype_mat), "samples ×", nrow(raw_subtype_mat), "genes\n")

# ── RAW: inner join of CPTAC breast + CPTAC lung + CCLE for breast_vs_lung
shared_blc <- Reduce(intersect, list(rownames(gm_breast), rownames(gm_lung), rownames(gm_ccle)))
raw_bvl_mat <- cbind(gm_breast[shared_blc, ], gm_lung[shared_blc, ], gm_ccle[shared_blc, ])
cat("  Raw (bvl):     ", ncol(raw_bvl_mat), "samples ×", nrow(raw_bvl_mat), "genes\n")

# ── CELLIGNER: already combined, subset to task-relevant samples
cat("  Celligner:     ", ncol(gm_cell), "samples ×", nrow(gm_cell), "genes\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4: Native-domain DA
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n  STEP 4: Native-domain DA\n", strrep("=", 70), "\n")

native_dir <- file.path(OUTDIR, "native_domain_da")

# ── Subtype: CPTAC native ──
cat("\n  [Native] CPTAC subtype DA...\n")
nat_cptac_sub_dir <- file.path(native_dir, "breast_subtype", "cptac")
dir.create(nat_cptac_sub_dir, recursive = TRUE, showWarnings = FALSE)
nat_sids <- intersect(cptac_subtype$sample_id, colnames(gm_breast))
nat_mat <- gm_breast[, nat_sids, drop = FALSE]
nat_grp <- cptac_subtype[match(nat_sids, sample_id), condition]
keep <- rowSums(!is.na(nat_mat)) >= ncol(nat_mat) * 0.5
nat_mat <- nat_mat[keep, , drop = FALSE]
for (i in seq_len(nrow(nat_mat))) {
  nas <- is.na(nat_mat[i, ]); if (any(nas)) nat_mat[i, nas] <- median(nat_mat[i, !nas])
}
cat("    ", sum(nat_grp == "Basal"), "Basal,", sum(nat_grp == "Luminal"), "Luminal,",
    nrow(nat_mat), "genes\n")
da_nat_cptac_sub <- run_limma(nat_mat, nat_grp, "Luminal_vs_Basal_CPTAC_native")
da_nat_cptac_sub[, inference_type := "native_domain_limma"]
fwrite(da_nat_cptac_sub, file.path(nat_cptac_sub_dir, "da_result.csv"))
fwrite(data.table(domain="CPTAC", n_basal=sum(nat_grp=="Basal"),
                  n_luminal=sum(nat_grp=="Luminal"), n_genes=nrow(nat_mat),
                  inference="native_domain_limma",
                  note="gene-level limma on PDC000120 gene_matrix"),
       file.path(nat_cptac_sub_dir, "sample_counts.csv"))

# ── Subtype: CCLE native ──
cat("  [Native] CCLE subtype DA...\n")
nat_ccle_sub_dir <- file.path(native_dir, "breast_subtype", "ccle")
dir.create(nat_ccle_sub_dir, recursive = TRUE, showWarnings = FALSE)
nat_ccle_sids <- intersect(ccle_subtype$sample_id, colnames(gm_ccle))
nat_ccle_mat <- gm_ccle[, nat_ccle_sids, drop = FALSE]
nat_ccle_grp <- ccle_subtype[match(nat_ccle_sids, sample_id), condition]
keep2 <- rowSums(!is.na(nat_ccle_mat)) >= ncol(nat_ccle_mat) * 0.5
nat_ccle_mat <- nat_ccle_mat[keep2, , drop = FALSE]
for (i in seq_len(nrow(nat_ccle_mat))) {
  nas <- is.na(nat_ccle_mat[i, ]); if (any(nas)) nat_ccle_mat[i, nas] <- median(nat_ccle_mat[i, !nas])
}
cat("    ", sum(nat_ccle_grp == "Basal"), "Basal,", sum(nat_ccle_grp == "Luminal"), "Luminal,",
    nrow(nat_ccle_mat), "genes\n")
da_nat_ccle_sub <- run_limma(nat_ccle_mat, nat_ccle_grp, "Luminal_vs_Basal_CCLE_native")
da_nat_ccle_sub[, inference_type := "native_domain_limma"]
fwrite(da_nat_ccle_sub, file.path(nat_ccle_sub_dir, "da_result.csv"))

# ── Breast vs Lung: CCLE native ──
cat("  [Native] CCLE breast vs lung DA...\n")
nat_ccle_bvl_dir <- file.path(native_dir, "breast_vs_lung", "ccle")
dir.create(nat_ccle_bvl_dir, recursive = TRUE, showWarnings = FALSE)
bvl_ccle_sids <- intersect(ccle_bvl$sample_id, colnames(gm_ccle))
if (length(bvl_ccle_sids) < 4) {
  cat("    skip CCLE breast vs lung: only", length(bvl_ccle_sids), "samples matched\n")
  da_nat_ccle_bvl <- data.table(gene = character(), logFC = numeric())
} else {
bvl_ccle_mat <- gm_ccle[, bvl_ccle_sids, drop = FALSE]
bvl_ccle_grp <- ccle_bvl[match(bvl_ccle_sids, sample_id), condition]
keep3 <- rowSums(!is.na(bvl_ccle_mat)) >= ncol(bvl_ccle_mat) * 0.3
bvl_ccle_mat <- bvl_ccle_mat[keep3, , drop = FALSE]
for (i in seq_len(nrow(bvl_ccle_mat))) {
  nas <- is.na(bvl_ccle_mat[i, ]); if (any(nas)) bvl_ccle_mat[i, nas] <- median(bvl_ccle_mat[i, !nas])
}
cat("    ", sum(bvl_ccle_grp == "Breast"), "Breast,", sum(bvl_ccle_grp == "Lung"), "Lung,",
    nrow(bvl_ccle_mat), "genes\n")
da_nat_ccle_bvl <- run_limma(bvl_ccle_mat, bvl_ccle_grp, "Lung_vs_Breast_CCLE_native")
da_nat_ccle_bvl[, inference_type := "native_domain_limma"]
fwrite(da_nat_ccle_bvl, file.path(nat_ccle_bvl_dir, "da_result.csv"))
}

# ── Breast vs Lung: CPTAC native (merged breast + lung) ──
cat("  [Native] CPTAC breast vs lung DA...\n")
nat_cptac_bvl_dir <- file.path(native_dir, "breast_vs_lung", "cptac")
dir.create(nat_cptac_bvl_dir, recursive = TRUE, showWarnings = FALSE)
bvl_cptac_sids <- intersect(cptac_bvl$sample_id, c(colnames(gm_breast), colnames(gm_lung)))
bvl_cptac_mat <- raw_bvl_mat[shared_blc, bvl_cptac_sids, drop = FALSE]
bvl_cptac_grp <- cptac_bvl[match(bvl_cptac_sids, sample_id), condition]
keep4 <- rowSums(!is.na(bvl_cptac_mat)) >= ncol(bvl_cptac_mat) * 0.3
bvl_cptac_mat <- bvl_cptac_mat[keep4, , drop = FALSE]
for (i in seq_len(nrow(bvl_cptac_mat))) {
  nas <- is.na(bvl_cptac_mat[i, ]); if (any(nas)) bvl_cptac_mat[i, nas] <- median(bvl_cptac_mat[i, !nas])
}
cat("    ", sum(bvl_cptac_grp == "Breast"), "Breast,", sum(bvl_cptac_grp == "Lung"), "Lung,",
    nrow(bvl_cptac_mat), "genes\n")
da_nat_cptac_bvl <- run_limma(bvl_cptac_mat, bvl_cptac_grp, "Lung_vs_Breast_CPTAC_native")
da_nat_cptac_bvl[, inference_type := "native_domain_limma"]
fwrite(da_nat_cptac_bvl, file.path(nat_cptac_bvl_dir, "da_result.csv"))
fwrite(data.table(domain="CPTAC", n_breast=sum(bvl_cptac_grp=="Breast"),
                  n_lung=sum(bvl_cptac_grp=="Lung"), n_genes=nrow(bvl_cptac_mat),
                  inference="native_domain_limma",
                  note="cross-study: cancer_type confounded with study"),
       file.path(nat_cptac_bvl_dir, "sample_counts.csv"))

# Native cross-domain agreement
cat("  [Native] Cross-domain agreement...\n")
ag_nat_sub <- fc_agreement(da_nat_cptac_sub, da_nat_ccle_sub)
ag_nat_sub[, c("representation", "task") := .("native", "breast_subtype")]
ag_nat_bvl <- fc_agreement(da_nat_cptac_bvl, da_nat_ccle_bvl)
ag_nat_bvl[, c("representation", "task") := .("native", "breast_vs_lung")]

# ═══════════════════════════════════════════════════════════════════════════
# STEP 5: Representation-level DA
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n  STEP 5: Representation-level DA\n", strrep("=", 70), "\n")

rep_dir <- file.path(OUTDIR, "representation_level_da")
subtype_markers <- c("FOXA1","GATA3","KRT5","KRT14","KRT17","EGFR","ESR1","PGR","ERBB2","CDH1")
bvl_markers <- c("NKX2-1","SFTPB","SFTPC","NAPSA","GATA3","FOXA1","ESR1","KRT19","EGFR","ERBB2","CDH1","VIM")

all_meta_subtype <- rbind(cptac_subtype[, .(sample_id, condition, domain)], ccle_subtype)
all_meta_bvl <- rbind(cptac_bvl, ccle_bvl)

# ── RAW ──────────────────────────────────────────────────────────────────
cat("\n── Raw representation ──\n")
cat("  Task A: breast subtype\n")
res_raw_sub <- run_task(raw_subtype_mat, all_meta_subtype, "breast_subtype", "raw",
                         subtype_markers, rep_dir)
cat("  Task B: breast vs lung\n")
res_raw_bvl <- run_task(raw_bvl_mat, all_meta_bvl, "breast_vs_lung", "raw",
                         bvl_markers, rep_dir)

# ── CELLIGNER ────────────────────────────────────────────────────────────
cat("\n── Celligner representation ──\n")
cat("  Task A: breast subtype\n")
res_cell_sub <- run_task(gm_cell, all_meta_subtype, "breast_subtype", "celligner",
                          subtype_markers, rep_dir)
cat("  Task B: breast vs lung\n")
res_cell_bvl <- run_task(gm_cell, all_meta_bvl, "breast_vs_lung", "celligner",
                          bvl_markers, rep_dir)

# ── BRIDGE SHIFT-ONLY ───────────────────────────────────────────────────
cat("\n── Bridge shift-only representation ──\n")
cat("  Task A: breast subtype\n")
res_bsh_sub <- run_task(gm_bridge_shift, all_meta_subtype, "breast_subtype", "bridge_shift",
                         subtype_markers, rep_dir)
cat("  Task B: breast vs lung\n")
res_bsh_bvl <- run_task(gm_bridge_shift, all_meta_bvl, "breast_vs_lung", "bridge_shift",
                         bvl_markers, rep_dir)

# ── BRIDGE SHIFT+SCALE ──────────────────────────────────────────────────
cat("\n── Bridge shift+scale representation ──\n")
cat("  Task A: breast subtype\n")
res_bsc_sub <- run_task(gm_bridge_scale, all_meta_subtype, "breast_subtype", "bridge_scale",
                         subtype_markers, rep_dir)
cat("  Task B: breast vs lung\n")
res_bsc_bvl <- run_task(gm_bridge_scale, all_meta_bvl, "breast_vs_lung", "bridge_scale",
                         bvl_markers, rep_dir)

# ═══════════════════════════════════════════════════════════════════════════
# STEP 6: Diagnostics — PCA structure plots
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n  STEP 6: Diagnostics\n", strrep("=", 70), "\n")

diag_dir <- file.path(OUTDIR, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

rep_mats <- list(raw_sub = raw_subtype_mat, celligner = gm_cell,
                  bridge_shift = gm_bridge_shift, bridge_scale = gm_bridge_scale)

# Subtype PCA
for (rep_name in c("raw", "celligner", "bridge_shift", "bridge_scale")) {
  mat <- if (rep_name == "raw") raw_subtype_mat else rep_mats[[rep_name]]
  meta_all <- all_meta_subtype
  sids <- intersect(meta_all$sample_id, colnames(mat))
  if (length(sids) < 4) { cat("  PCA skip: ", rep_name, " subtype (too few)\n"); next }
  meta_use <- meta_all[match(sids, sample_id)]
  m <- mat[, sids, drop = FALSE]
  keep <- rowSums(!is.na(m)) >= ncol(m) * 0.3
  m <- m[keep, ]; m[is.na(m)] <- 0
  tryCatch({
    pca_plot(m, meta_use, "condition", "domain",
             paste0(rep_name, " — Breast Subtype PCA"),
             file.path(diag_dir, paste0("pca_subtype_", rep_name, ".png")))
    cat("  PCA: ", rep_name, " subtype\n")
  }, error = function(e) warning(e$message))
}

# BvL PCA
for (rep_name in c("raw", "celligner", "bridge_shift", "bridge_scale")) {
  mat <- if (rep_name == "raw") raw_bvl_mat else rep_mats[[rep_name]]
  meta_all <- all_meta_bvl
  sids <- intersect(meta_all$sample_id, colnames(mat))
  if (length(sids) < 4) { cat("  PCA skip: ", rep_name, " bvl (too few)\n"); next }
  meta_use <- meta_all[match(sids, sample_id)]
  m <- mat[, sids, drop = FALSE]
  keep <- rowSums(!is.na(m)) >= ncol(m) * 0.3
  m <- m[keep, ]; m[is.na(m)] <- 0
  tryCatch({
    pca_plot(m, meta_use, "condition", "domain",
             paste0(rep_name, " — Breast vs Lung PCA"),
             file.path(diag_dir, paste0("pca_bvl_", rep_name, ".png")))
    cat("  PCA: ", rep_name, " bvl\n")
  }, error = function(e) warning(e$message))
}

# Domain R²
cat("\n  Domain R² summary:\n")
dr2_rows <- list()
all_reps <- c("raw", "celligner", "bridge_shift", "bridge_scale")
for (rep_name in all_reps) {
  for (task in c("subtype", "bvl")) {
    meta <- if(task=="subtype") all_meta_subtype else all_meta_bvl
    if (rep_name == "raw" && task == "subtype") mat <- raw_subtype_mat
    else if (rep_name == "raw" && task == "bvl") mat <- raw_bvl_mat
    else mat <- rep_mats[[rep_name]]
    sids <- intersect(meta$sample_id, colnames(mat))
    if (length(sids) < 4) { cat("    ", rep_name, "/", task, ": skip (too few samples)\n"); next }
    m <- mat[, sids, drop = FALSE]; m[is.na(m)] <- 0
    d <- meta[match(sids, sample_id), domain]
    r2 <- domain_r2(m, d)
    cat("    ", rep_name, "/", task, ": PC1 domain R² =", round(r2, 4), "\n")
    dr2_rows[[length(dr2_rows) + 1]] <- data.table(
      representation = rep_name,
      task = ifelse(task == "subtype", "breast_subtype", "breast_vs_lung"),
      pc1_domain_r2 = r2
    )
  }
}
fwrite(rbindlist(dr2_rows), file.path(diag_dir, "domain_r2_summary.csv"))

# Marker profile plots
prof_dir <- file.path(diag_dir, "marker_profiles")
dir.create(prof_dir, showWarnings = FALSE)
for (g in subtype_markers) {
  for (rn in all_reps) {
    mat <- if (rn == "raw") raw_subtype_mat else rep_mats[[rn]]
    sids <- intersect(all_meta_subtype$sample_id, colnames(mat))
    if (length(sids) < 4) next
    m <- mat[, sids, drop = FALSE]; m[is.na(m)] <- 0
    meta_use <- all_meta_subtype[match(sids, sample_id)]
    tryCatch(marker_profile(m, meta_use, g,
                            file.path(prof_dir, paste0(g, "_subtype_", rn, ".png"))),
             error = function(e) NULL)
  }
}
for (g in bvl_markers) {
  for (rn in all_reps) {
    mat <- if (rn == "raw") raw_bvl_mat else rep_mats[[rn]]
    sids <- intersect(all_meta_bvl$sample_id, colnames(mat))
    if (length(sids) < 4) next
    m <- mat[, sids, drop = FALSE]; m[is.na(m)] <- 0
    meta_use <- all_meta_bvl[match(sids, sample_id)]
    tryCatch(marker_profile(m, meta_use, g,
                            file.path(prof_dir, paste0(g, "_bvl_", rn, ".png"))),
             error = function(e) NULL)
  }
}

# ═══════════════════════════════════════════════════════════════════════════
# STEP 7: Comparison table
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n  STEP 7: Comparison table\n", strrep("=", 70), "\n")

ag_files <- Sys.glob(file.path(rep_dir, "*", "*", "cross_domain_agreement.csv"))
all_ag <- rbindlist(c(
  lapply(ag_files, fread),
  list(ag_nat_sub, ag_nat_bvl)
), fill = TRUE)
fwrite(all_ag, file.path(OUTDIR, "benchmark_method_comparison.tsv"), sep = "\t")
cat("\nBenchmark comparison table:\n")
print(all_ag[, .(representation, task, n, pearson_r, spearman_rho, direction_agree, rmse)])

# ═══════════════════════════════════════════════════════════════════════════
# STEP 8: FC scatter plots
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n  STEP 8: FC scatter plots\n", strrep("=", 70), "\n")

scatter_dir <- file.path(OUTDIR, "diagnostics", "fc_scatters")
dir.create(scatter_dir, showWarnings = FALSE)

for (rn in c("raw", "celligner", "bridge_shift", "bridge_scale")) {
  for (task in c("breast_subtype", "breast_vs_lung")) {
    fpath <- file.path(rep_dir, rn, task, "fc_scatter_data.csv")
    if (!file.exists(fpath)) next
    fcd <- fread(fpath)
    ag <- fread(file.path(rep_dir, rn, task, "cross_domain_agreement.csv"))

    png(file.path(scatter_dir, paste0("fc_scatter_", rn, "_", task, ".png")),
        width = 800, height = 700, res = 130)
    par(mar = c(5, 5, 4, 1))
    plot(fcd$fc_cptac, fcd$fc_ccle, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.15),
         xlab = "CPTAC log2FC", ylab = "CCLE log2FC",
         main = paste0(rn, " — ", gsub("_", " ", task)))
    abline(0, 1, col = "red", lwd = 1.5, lty = 2)
    abline(h = 0, v = 0, col = "gray70")
    legend("topleft", legend = c(
      paste0("r = ", round(ag$pearson_r, 3)),
      paste0("ρ = ", round(ag$spearman_rho, 3)),
      paste0("dir = ", round(ag$direction_agree * 100, 1), "%"),
      paste0("n = ", ag$n)
    ), bty = "n", cex = 0.9)
    dev.off()
    cat("  ", rn, "/", task, "\n")
  }
}

# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 70), "\n")
cat("  BENCHMARK COMPLETE\n")
cat(strrep("=", 70), "\n")
cat("  Outputs:", OUTDIR, "\n")
cat("  Key file: benchmark_method_comparison.tsv\n\n")
