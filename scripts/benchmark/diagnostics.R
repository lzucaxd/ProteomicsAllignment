#!/usr/bin/env Rscript
# =============================================================================
# Diagnostics — assumption checks, structure plots, profile plots
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# PCA structure plot
# ---------------------------------------------------------------------------
plot_pca_structure <- function(matrix, sample_meta,
                                color_by = "condition", shape_by = "domain",
                                title = "PCA Structure", outpath = NULL) {
  complete_cols <- complete.cases(t(matrix))
  mat_clean <- matrix[, complete_cols, drop = FALSE]
  meta_clean <- sample_meta[complete_cols, ]

  gene_var <- apply(mat_clean, 1, var, na.rm = TRUE)
  top_genes <- order(gene_var, decreasing = TRUE)[1:min(2000, nrow(mat_clean))]
  mat_top <- mat_clean[top_genes, ]
  mat_top[is.na(mat_top)] <- 0

  pca <- prcomp(t(mat_top), center = TRUE, scale. = FALSE)
  pve <- round(summary(pca)$importance[2, 1:2] * 100, 1)

  pc_df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    color = meta_clean[[color_by]],
    shape = if (!is.null(shape_by) && shape_by %in% names(meta_clean)) meta_clean[[shape_by]] else "all"
  )

  if (!is.null(outpath)) {
    png(outpath, width = 900, height = 700, res = 120)
    par(mar = c(5, 5, 3, 8), xpd = TRUE)
    cols <- as.integer(factor(pc_df$color))
    palette <- c("steelblue", "tomato", "forestgreen", "purple", "orange", "brown")
    pchs <- c(16, 17, 15, 18)[as.integer(factor(pc_df$shape))]
    plot(pc_df$PC1, pc_df$PC2,
         col = palette[cols], pch = pchs,
         xlab = paste0("PC1 (", pve[1], "%)"),
         ylab = paste0("PC2 (", pve[2], "%)"),
         main = title, cex = 1.2)
    legend("topright", inset = c(-0.25, 0),
           legend = levels(factor(pc_df$color)),
           col = palette[seq_along(levels(factor(pc_df$color)))],
           pch = 16, cex = 0.8, title = color_by)
    if (!is.null(shape_by)) {
      legend("bottomright", inset = c(-0.25, 0),
             legend = levels(factor(pc_df$shape)),
             pch = c(16, 17, 15, 18)[seq_along(levels(factor(pc_df$shape)))],
             cex = 0.8, title = shape_by)
    }
    dev.off()
    message("PCA plot saved: ", outpath)
  }
  invisible(list(pca = pca, pve = pve, pc_df = pc_df))
}

# ---------------------------------------------------------------------------
# Marker profile plot (one gene, split by domain × condition)
# ---------------------------------------------------------------------------
plot_marker_profile <- function(matrix, sample_meta, gene,
                                 condition_col = "condition", domain_col = "domain",
                                 outpath = NULL) {
  if (!gene %in% rownames(matrix)) {
    message("  Gene ", gene, " not in matrix, skipping.")
    return(invisible(NULL))
  }
  vals <- matrix[gene, ]
  sm <- as.data.table(sample_meta)
  sm[, value := vals]
  sm[, group := paste(get(domain_col), get(condition_col), sep = "\n")]

  if (!is.null(outpath)) {
    png(outpath, width = 700, height = 500, res = 120)
    par(mar = c(7, 5, 3, 1))
    groups <- sort(unique(sm$group))
    group_vals <- lapply(groups, function(g) sm[group == g, value])
    boxplot(group_vals, names = groups, las = 2,
            main = gene, ylab = "Abundance (log2)",
            col = c("steelblue", "tomato", "lightblue", "salmon")[seq_along(groups)],
            outline = TRUE)
    stripchart(group_vals, vertical = TRUE, method = "jitter",
               pch = 16, cex = 0.6, col = "gray30", add = TRUE)
    dev.off()
    message("  Profile plot: ", outpath)
  }
  invisible(sm)
}

# ---------------------------------------------------------------------------
# Generate all diagnostics for a task result
# ---------------------------------------------------------------------------
generate_diagnostics <- function(matrix, sample_meta, representation_name,
                                  task_name, outdir, marker_genes = NULL) {
  diag_dir <- file.path(outdir, "diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  sm <- as.data.table(sample_meta)

  # PCA by condition
  tryCatch(
    plot_pca_structure(matrix, sm, color_by = "condition", shape_by = "domain",
                       title = paste(representation_name, task_name, "— PCA"),
                       outpath = file.path(diag_dir, "pca_condition_domain.png")),
    error = function(e) warning("PCA plot failed: ", e$message)
  )

  # PCA by domain only
  tryCatch(
    plot_pca_structure(matrix, sm, color_by = "domain", shape_by = NULL,
                       title = paste(representation_name, "— Domain PCA"),
                       outpath = file.path(diag_dir, "pca_domain.png")),
    error = function(e) warning("Domain PCA failed: ", e$message)
  )

  # Domain effect metric
  de <- compute_domain_effect(matrix, sm$domain)
  writeLines(c(
    paste("PC1 domain R²:", round(de$pc1_domain_r2, 4)),
    paste("PC1 variance explained:", round(de$pve[1] * 100, 2), "%")
  ), file.path(diag_dir, "domain_effect.txt"))

  # Spread summary
  spread <- compute_spread_summary(matrix, sm$domain)
  fwrite(spread, file.path(diag_dir, "spread_by_domain.csv"))

  # Marker profiles
  if (!is.null(marker_genes)) {
    prof_dir <- file.path(diag_dir, "marker_profiles")
    dir.create(prof_dir, showWarnings = FALSE)
    for (g in marker_genes) {
      tryCatch(
        plot_marker_profile(matrix, sm, g,
                            outpath = file.path(prof_dir, paste0("profile_", g, ".png"))),
        error = function(e) warning("Marker profile failed for ", g, ": ", e$message)
      )
    }
  }

  message("Diagnostics saved: ", diag_dir)
  invisible(diag_dir)
}
