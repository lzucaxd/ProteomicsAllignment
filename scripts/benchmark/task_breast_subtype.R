#!/usr/bin/env Rscript
# =============================================================================
# Task A — Breast Subtype (Luminal vs Basal)
# =============================================================================
# Runs representation-level limma separately within CPTAC and CCLE,
# then compares cross-domain agreement.
# =============================================================================

SUBTYPE_MARKERS <- c("FOXA1", "GATA3", "KRT5", "KRT14", "KRT17",
                      "EGFR", "ESR1", "PGR", "ERBB2", "CDH1")

run_task_breast_subtype <- function(matrix, sample_meta, feature_meta,
                                     representation_name, outdir,
                                     marker_genes = NULL, ...) {
  if (is.null(marker_genes)) marker_genes <- SUBTYPE_MARKERS
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  sm <- as.data.table(sample_meta)
  if (!"condition" %in% names(sm) || !"domain" %in% names(sm))
    stop("sample_meta must have 'condition' and 'domain' columns")

  # Subset to breast subtype samples
  subtype_samples <- sm[tolower(condition) %in% c("basal", "luminal")]
  if (nrow(subtype_samples) == 0) stop("No Basal/Luminal samples found in sample_meta")

  results <- list()

  for (dom in c("CPTAC", "CCLE")) {
    dom_samples <- subtype_samples[toupper(domain) == dom]
    if (nrow(dom_samples) < 4) {
      message("  Skipping ", dom, ": only ", nrow(dom_samples), " subtype samples")
      next
    }

    sample_ids <- intersect(dom_samples$sample_id, colnames(matrix))
    if (length(sample_ids) < 4) {
      message("  Skipping ", dom, ": only ", length(sample_ids), " samples found in matrix")
      next
    }

    dom_mat <- matrix[, sample_ids, drop = FALSE]
    dom_groups <- dom_samples[match(sample_ids, sample_id), condition]

    # Remove genes with >50% NA in either group
    for (g in unique(dom_groups)) {
      g_cols <- which(dom_groups == g)
      na_frac <- rowMeans(is.na(dom_mat[, g_cols, drop = FALSE]))
      dom_mat[na_frac > 0.5, g_cols] <- NA
    }
    keep_genes <- rowSums(!is.na(dom_mat)) >= ncol(dom_mat) * 0.5
    dom_mat <- dom_mat[keep_genes, , drop = FALSE]

    message("  ", dom, ": ", ncol(dom_mat), " samples, ", nrow(dom_mat), " genes")

    # Run limma (contrast: Luminal - Basal)
    da <- run_limma_da(dom_mat, dom_groups,
                        contrast_name = paste0("Luminal_vs_Basal_", dom))

    # Marker summary
    markers <- extract_marker_summary(da, marker_genes)

    # Save
    dom_dir <- file.path(outdir, tolower(dom))
    dir.create(dom_dir, showWarnings = FALSE)
    fwrite(da, file.path(dom_dir, "da_limma_result.csv"))
    fwrite(markers, file.path(dom_dir, "marker_summary.csv"))

    sample_note <- data.table(
      domain = dom,
      n_basal = sum(dom_groups == "Basal"),
      n_luminal = sum(dom_groups == "Luminal"),
      n_genes = nrow(dom_mat),
      inference_type = "representation_level_limma",
      representation = representation_name
    )
    fwrite(sample_note, file.path(dom_dir, "sample_counts.csv"))

    results[[dom]] <- list(da = da, markers = markers, sample_note = sample_note)
  }

  # Cross-domain agreement
  if (!is.null(results[["CPTAC"]]) && !is.null(results[["CCLE"]])) {
    agreement <- compute_fc_agreement(results$CPTAC$da, results$CCLE$da)
    agreement_summary <- data.table(
      representation = representation_name,
      task = "breast_subtype",
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

  # Marker direction check
  marker_directions <- data.table(
    gene = c("FOXA1", "GATA3", "ESR1", "PGR", "ERBB2", "CDH1",
             "KRT5", "KRT14", "KRT17", "EGFR"),
    expected_sign = c(1, 1, 1, 1, 1, 1, -1, -1, -1, -1)
  )
  for (dom in names(results)) {
    if (dom == "agreement") next
    dc <- check_marker_directions(results[[dom]]$da, marker_directions)
    fwrite(dc, file.path(outdir, tolower(dom), "marker_direction_check.csv"))
  }

  message("Task breast_subtype complete for ", representation_name, " → ", outdir)
  invisible(results)
}
