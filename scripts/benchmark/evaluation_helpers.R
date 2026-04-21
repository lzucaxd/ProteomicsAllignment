#!/usr/bin/env Rscript
# =============================================================================
# Evaluation helpers — limma DA, marker summaries, metrics
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ---------------------------------------------------------------------------
# Representation-level limma DA
# ---------------------------------------------------------------------------
run_limma_da <- function(matrix, groups, contrast_name = "GroupB_vs_GroupA",
                          covariate_formula = NULL) {
  if (length(groups) != ncol(matrix))
    stop("groups length (", length(groups), ") != matrix columns (", ncol(matrix), ")")

  group_factor <- factor(groups)
  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- levels(group_factor)

  fit <- lmFit(matrix, design)

  lvls <- levels(group_factor)
  if (length(lvls) != 2) stop("Exactly 2 group levels required, got: ", paste(lvls, collapse = ", "))
  contrast_str <- paste0(lvls[2], "-", lvls[1])
  cm <- makeContrasts(contrasts = contrast_str, levels = design)

  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)

  tt <- topTable(fit2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$contrast <- contrast_name
  tt$inference_type <- "representation_level_limma"
  setDT(tt)
  setcolorder(tt, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "contrast", "inference_type"))
  tt
}

# ---------------------------------------------------------------------------
# Marker summary: extract rows for canonical genes
# ---------------------------------------------------------------------------
extract_marker_summary <- function(da_result, marker_genes) {
  da_result[gene %in% marker_genes]
}

# ---------------------------------------------------------------------------
# FC agreement metrics between two DA results (same gene set)
# ---------------------------------------------------------------------------
compute_fc_agreement <- function(da_cptac, da_ccle, fc_col = "logFC") {
  shared <- merge(da_cptac[, .(gene, fc_cptac = get(fc_col))],
                  da_ccle[, .(gene, fc_ccle = get(fc_col))],
                  by = "gene")
  shared <- shared[is.finite(fc_cptac) & is.finite(fc_ccle)]
  if (nrow(shared) < 3) return(list(n = nrow(shared)))

  list(
    n = nrow(shared),
    pearson_r = cor(shared$fc_cptac, shared$fc_ccle, method = "pearson"),
    spearman_rho = cor(shared$fc_cptac, shared$fc_ccle, method = "spearman"),
    direction_agree_frac = mean(sign(shared$fc_cptac) == sign(shared$fc_ccle)),
    median_abs_fc_diff = median(abs(shared$fc_cptac - shared$fc_ccle)),
    rmse = sqrt(mean((shared$fc_cptac - shared$fc_ccle)^2)),
    fc_data = shared
  )
}

# ---------------------------------------------------------------------------
# Spread / variance summary per domain
# ---------------------------------------------------------------------------
compute_spread_summary <- function(matrix, domain_labels) {
  domains <- unique(domain_labels)
  out <- list()
  for (d in domains) {
    cols <- which(domain_labels == d)
    gene_var <- apply(matrix[, cols, drop = FALSE], 1, var, na.rm = TRUE)
    out[[d]] <- data.table(gene = rownames(matrix), domain = d, variance = gene_var)
  }
  rbindlist(out)
}

# ---------------------------------------------------------------------------
# Domain effect: fraction of variance on PC1 explained by domain
# ---------------------------------------------------------------------------
compute_domain_effect <- function(matrix, domain_labels) {
  complete <- complete.cases(t(matrix))
  if (sum(complete) < 10) return(list(pc1_domain_r2 = NA_real_))
  mat_clean <- matrix[, complete]
  domains_clean <- domain_labels[complete]

  pca <- prcomp(t(mat_clean), center = TRUE, scale. = FALSE, rank. = 5)
  pc1 <- pca$x[, 1]
  r2 <- summary(lm(pc1 ~ factor(domains_clean)))$r.squared
  list(pc1_domain_r2 = r2, pca = pca, pve = summary(pca)$importance[2, 1:5])
}

# ---------------------------------------------------------------------------
# Marker direction check against known panel
# ---------------------------------------------------------------------------
check_marker_directions <- function(da_result, expected_directions) {
  # expected_directions: data.frame with gene, expected_sign (+1 or -1)
  merged <- merge(da_result[, .(gene, logFC)],
                  as.data.table(expected_directions),
                  by = "gene")
  merged[, observed_sign := sign(logFC)]
  merged[, direction_correct := (observed_sign == expected_sign)]
  merged
}
