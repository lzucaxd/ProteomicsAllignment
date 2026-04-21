#!/usr/bin/env Rscript
# =============================================================================
# Task B — Breast vs Lung
# =============================================================================
# Runs representation-level limma separately within CPTAC and CCLE,
# then compares cross-domain agreement.
# =============================================================================

LINEAGE_MARKERS <- c("NKX2-1", "SFTPB", "SFTPC", "NAPSA",  # lung
                      "GATA3", "FOXA1", "ESR1", "KRT19",     # breast
                      "EGFR", "ERBB2", "CDH1", "VIM")        # shared / EMT

run_task_breast_vs_lung <- function(matrix, sample_meta, feature_meta,
                                     representation_name, outdir,
                                     marker_genes = NULL, ...) {
  if (is.null(marker_genes)) marker_genes <- LINEAGE_MARKERS
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  sm <- as.data.table(sample_meta)
  if (!"condition" %in% names(sm) || !"domain" %in% names(sm))
    stop("sample_meta must have 'condition' and 'domain' columns")

  lineage_samples <- sm[tolower(condition) %in% c("breast", "lung")]
  if (nrow(lineage_samples) == 0) stop("No Breast/Lung samples found in sample_meta")

  results <- list()

  for (dom in c("CPTAC", "CCLE")) {
    dom_samples <- lineage_samples[toupper(domain) == dom]
    if (nrow(dom_samples) < 4) {
      message("  Skipping ", dom, ": only ", nrow(dom_samples), " lineage samples")
      next
    }

    sample_ids <- intersect(dom_samples$sample_id, colnames(matrix))
    if (length(sample_ids) < 4) next

    dom_mat <- matrix[, sample_ids, drop = FALSE]
    dom_groups <- dom_samples[match(sample_ids, sample_id), condition]

    keep_genes <- rowSums(!is.na(dom_mat)) >= ncol(dom_mat) * 0.3
    dom_mat <- dom_mat[keep_genes, , drop = FALSE]

    message("  ", dom, ": ", ncol(dom_mat), " samples, ", nrow(dom_mat), " genes")

    da <- run_limma_da(dom_mat, dom_groups,
                        contrast_name = paste0("Breast_vs_Lung_", dom))

    markers <- extract_marker_summary(da, marker_genes)

    dom_dir <- file.path(outdir, tolower(dom))
    dir.create(dom_dir, showWarnings = FALSE)
    fwrite(da, file.path(dom_dir, "da_limma_result.csv"))
    fwrite(markers, file.path(dom_dir, "marker_summary.csv"))

    sample_note <- data.table(
      domain = dom,
      n_breast = sum(tolower(dom_groups) == "breast"),
      n_lung = sum(tolower(dom_groups) == "lung"),
      n_genes = nrow(dom_mat),
      inference_type = "representation_level_limma",
      representation = representation_name,
      caveat = if (dom == "CPTAC") "cancer_type confounded with study" else "unbalanced design"
    )
    fwrite(sample_note, file.path(dom_dir, "sample_counts.csv"))

    results[[dom]] <- list(da = da, markers = markers, sample_note = sample_note)
  }

  # Cross-domain agreement
  if (!is.null(results[["CPTAC"]]) && !is.null(results[["CCLE"]])) {
    agreement <- compute_fc_agreement(results$CPTAC$da, results$CCLE$da)
    agreement_summary <- data.table(
      representation = representation_name,
      task = "breast_vs_lung",
      n_shared_genes = agreement$n,
      pearson_r = agreement$pearson_r,
      spearman_rho = agreement$spearman_rho,
      direction_agree_frac = agreement$direction_agree_frac,
      median_abs_fc_diff = agreement$median_abs_fc_diff,
      rmse = agreement$rmse
    )
    fwrite(agreement_summary, file.path(outdir, "cross_domain_agreement.csv"))
    if (!is.null(agreement$fc_data))
      fwrite(agreement$fc_data, file.path(outdir, "fc_scatter_data.csv"))
    results$agreement <- agreement_summary
  }

  message("Task breast_vs_lung complete for ", representation_name, " → ", outdir)
  invisible(results)
}
