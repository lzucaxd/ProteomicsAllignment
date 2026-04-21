#!/usr/bin/env Rscript
# Generate report figures into presentation_materials/figures/report/
#
# Run from repo root:
#   Rscript --vanilla scripts/presentation/generate_report_figures.R
#
# Reads:
# - reports/benchmark_master/benchmark_results/{method}/{task}/representation_da/{cptac,ccle}/da_limma_result.csv
# - data/processed/methods/{method}/transformed_{task}.csv
# - data/processed/union/sample_meta_{task}.csv
# - reports/benchmark_master/benchmark_results/comparison_summary.csv
# - reports/benchmark_master/benchmark_results/disconnect_scores.csv
# - reports/benchmark_master/diagnostics/fc_stratified_{task}.csv
# - reports/benchmark_master/benchmark_results/raw/{task}/calibration/null_distribution.csv
# - reports/benchmark_master/benchmark_results/raw/{task}/calibration/observed_metrics.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

HAS_GGREP <- requireNamespace("ggrepel", quietly = TRUE)

ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) {
  dirname(normalizePath(sub("^--file=", "", ff[1L])))
} else {
  normalizePath(file.path(getwd(), "scripts", "presentation"), mustWork = FALSE)
}
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

OUT <- file.path(PRES_OUT, "figures", "report")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

save_plot_both <- function(plot, base_name, width, height, dpi = 300) {
  pdf_path <- file.path(OUT, paste0(base_name, ".pdf"))
  png_path <- file.path(OUT, paste0(base_name, ".png"))
  ggplot2::ggsave(filename = pdf_path, plot = plot, width = width, height = height)
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

BENCH <- file.path(REPO, "reports", "benchmark_master", "benchmark_results")
METHODS <- file.path(REPO, "data", "processed", "methods")
UNION <- file.path(REPO, "data", "processed", "union")

METHOD_IDS <- c("raw", "bridge_shift", "bridge_scale", "celligner")
METHOD_LABELS <- c(
  raw = "Raw",
  bridge_shift = "Bridge shift",
  bridge_scale = "Bridge shift+scale",
  celligner = "Celligner"
)
METHOD_COLORS <- c(
  Raw = "#AAAAAA",
  `Bridge shift` = "#4878CF",
  `Bridge shift+scale` = "#6A9BD1",
  Celligner = "#D65F5F"
)

gene_col_from_dt <- function(dt) {
  cand <- intersect(c("gene", "Gene", "GeneSymbol", "Protein", "protein"), names(dt))
  if (length(cand) == 0L) stop("No gene column found in: ", paste(names(dt), collapse = ", "))
  cand[[1L]]
}

se_from_da <- function(dt) {
  if ("SE" %in% names(dt)) return(as.numeric(dt$SE))
  if (all(c("logFC", "t") %in% names(dt))) return(as.numeric(abs(dt$logFC / dt$t)))
  rep(NA_real_, nrow(dt))
}

load_da <- function(method, task, domain) {
  path <- file.path(
    BENCH, method, task, "representation_da", domain, "da_limma_result.csv"
  )
  if (!file.exists(path)) {
    alt <- file.path(BENCH, method, task, "representation_da", sprintf("da_%s.csv", domain))
    if (file.exists(alt)) path <- alt else return(NULL)
  }
  fread(path)
}

load_matrix <- function(method, task) {
  path <- file.path(METHODS, method, sprintf("transformed_%s.csv", task))
  if (!file.exists(path)) return(NULL)
  fread(path)
}

load_meta <- function(task) {
  fread(file.path(UNION, sprintf("sample_meta_%s.csv", task)))
}

method_display <- function(m) {
  m <- as.character(m)
  out <- unname(METHOD_LABELS[m])
  out[is.na(out)] <- m[is.na(out)]
  out
}

stratum_label <- function(x) {
  x <- as.character(x)
  y <- fifelse(
    x == "sig_both", "Sig both",
    fifelse(
      x == "sig_cptac_only", "Sig CPTAC only",
      fifelse(x == "sig_ccle_only", "Sig CCLE only", "Sig neither")
    )
  )
  factor(y, levels = c("Sig neither", "Sig CCLE only", "Sig CPTAC only", "Sig both"))
}

figure_se_comparison <- function(task) {
  cptac_da <- load_da("raw", task, "cptac")
  ccle_da <- load_da("raw", task, "ccle")
  if (is.null(cptac_da) || is.null(ccle_da)) return(invisible(NULL))

  gc <- gene_col_from_dt(cptac_da)
  ge <- gene_col_from_dt(ccle_da)

  merged <- merge(
    cptac_da[, .(gene = as.character(get(gc)), SE_cptac = se_from_da(.SD))],
    ccle_da[, .(gene = as.character(get(ge)), SE_ccle = se_from_da(.SD))],
    by = "gene"
  )
  merged <- merged[is.finite(SE_cptac) & is.finite(SE_ccle)]

  p <- ggplot(merged, aes(x = SE_cptac, y = SE_ccle)) +
    geom_point(alpha = 0.15, size = 0.6, color = "grey30") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    labs(
      x = "CPTAC SE (raw limma; or |logFC/t| if SE missing)",
      y = "CCLE SE (raw limma; or |logFC/t| if SE missing)",
      title = sprintf("SE comparison (%s)", task)
    ) +
    coord_equal() +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  save_plot_both(p, sprintf("se_comparison_%s", task), width = 5, height = 5)
}

figure_fc_scatter_panel <- function(task) {
  strat_colors <- c(
    "Sig both" = "#D62728",
    "Sig CPTAC only" = "#4878CF",
    "Sig CCLE only" = "#6ACC65",
    "Sig neither" = "#CCCCCC"
  )

  panels <- list()
  for (m in METHOD_IDS) {
    cptac_da <- load_da(m, task, "cptac")
    ccle_da <- load_da(m, task, "ccle")
    if (is.null(cptac_da) || is.null(ccle_da)) next

    gc <- gene_col_from_dt(cptac_da)
    ge <- gene_col_from_dt(ccle_da)
    stopifnot("logFC" %in% names(cptac_da), "logFC" %in% names(ccle_da))
    stopifnot("adj.P.Val" %in% names(cptac_da), "adj.P.Val" %in% names(ccle_da))

    merged <- merge(
      cptac_da[, .(gene = as.character(get(gc)), logFC_cptac = logFC, pval_cptac = adj.P.Val)],
      ccle_da[, .(gene = as.character(get(ge)), logFC_ccle = logFC, pval_ccle = adj.P.Val)],
      by = "gene"
    )

    merged[, stratum := fifelse(
      pval_cptac < 0.05 & pval_ccle < 0.05, "Sig both",
      fifelse(
        pval_cptac < 0.05, "Sig CPTAC only",
        fifelse(pval_ccle < 0.05, "Sig CCLE only", "Sig neither")
      )
    )]
    merged[, stratum := stratum_label(stratum)]
    setorder(merged, stratum)

    r_val <- cor(merged$logFC_cptac, merged$logFC_ccle, use = "complete.obs")

    p <- ggplot(merged, aes(x = logFC_cptac, y = logFC_ccle, color = stratum)) +
      geom_point(alpha = 0.3, size = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
      scale_color_manual(values = strat_colors, name = "Stratum") +
      annotate(
        "text", x = Inf, y = -Inf,
        label = sprintf("r = %.3f", r_val),
        hjust = 1.1, vjust = -0.5, size = 3.5, fontface = "bold"
      ) +
      labs(
        x = "CPTAC logFC",
        y = "CCLE logFC",
        title = method_display(m)
      ) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11), legend.position = "none")

    panels[[m]] <- p
  }

  if (length(panels) == 0L) return(invisible(NULL))
  ord <- METHOD_IDS[METHOD_IDS %in% names(panels)]
  panels <- panels[ord]

  p_legend <- panels[[1L]] +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 1L, override.aes = list(size = 3, alpha = 1)))

  legend <- get_legend(p_legend)
  grid <- plot_grid(plotlist = panels, nrow = 1L, align = "hv")
  final <- plot_grid(grid, legend, ncol = 1L, rel_heights = c(1, 0.10))

  save_plot_both(
    final,
    sprintf("fc_scatter_panel_%s", task),
    width = 4 * length(panels),
    height = 4.8
  )
}

figure_pca_panel <- function(task) {
  meta <- load_meta(task)
  meta[, domain := toupper(as.character(domain))]

  raw_mat <- load_matrix("raw", task)
  if (is.null(raw_mat)) return(invisible(NULL))

  gene_col <- gene_col_from_dt(raw_mat)
  sample_cols <- intersect(setdiff(names(raw_mat), gene_col), meta$sample_id)
  if (length(sample_cols) < 5L) return(invisible(NULL))

  mat_numeric <- as.matrix(raw_mat[, ..sample_cols])
  rownames(mat_numeric) <- as.character(raw_mat[[gene_col]])

  complete_rows <- complete.cases(mat_numeric)
  if (sum(complete_rows) < 10L) return(invisible(NULL))
  mat_complete <- mat_numeric[complete_rows, , drop = FALSE]

  # prcomp expects samples as rows, genes as cols
  pca_fit <- prcomp(t(mat_complete), center = TRUE, scale. = FALSE)
  ve <- summary(pca_fit)$importance[2, 1:2] * 100

  panels <- list()
  for (m in METHOD_IDS) {
    m_mat <- load_matrix(m, task)
    if (is.null(m_mat)) next

    gene_col_m <- gene_col_from_dt(m_mat)
    genes <- intersect(rownames(mat_complete), as.character(m_mat[[gene_col_m]]))
    if (length(genes) < 10L) next

    idx <- match(genes, as.character(m_mat[[gene_col_m]]))
    m_numeric <- as.matrix(m_mat[idx, ..sample_cols])
    rownames(m_numeric) <- genes

    projected <- predict(pca_fit, newdata = t(m_numeric))

    pca_df <- data.table(
      PC1 = projected[, 1L],
      PC2 = projected[, 2L],
      sample_id = sample_cols
    )
    pca_df <- merge(pca_df, meta, by = "sample_id", all.x = TRUE)

    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = domain, shape = condition)) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_manual(values = c(CPTAC = "#4878CF", CCLE = "#D65F5F")) +
      labs(
        x = sprintf("PC1 (%.1f%%)", ve[[1L]]),
        y = sprintf("PC2 (%.1f%%)", ve[[2L]]),
        title = method_display(m),
        subtitle = "PC basis fit on Raw; other methods projected"
      ) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold"), legend.position = "none")

    panels[[m]] <- p
  }

  if (length(panels) == 0L) return(invisible(NULL))
  ord <- METHOD_IDS[METHOD_IDS %in% names(panels)]
  panels <- panels[ord]

  p_legend <- panels[[1L]] +
    theme(legend.position = "bottom") +
    guides(
      color = guide_legend(title = "Domain"),
      shape = guide_legend(title = "Condition")
    )
  legend <- get_legend(p_legend)
  grid <- plot_grid(plotlist = panels, nrow = 1L, align = "hv")
  final <- plot_grid(grid, legend, ncol = 1L, rel_heights = c(1, 0.12))

  save_plot_both(
    final,
    sprintf("pca_panel_%s", task),
    width = 4 * length(panels),
    height = 5.2
  )
}

figure_permutation_null <- function(task) {
  null_path <- file.path(BENCH, "raw", task, "calibration", "null_distribution.csv")
  obs_path <- file.path(BENCH, "raw", task, "calibration", "observed_metrics.csv")
  if (!file.exists(null_path) || !file.exists(obs_path)) return(invisible(NULL))

  null_dt <- fread(null_path)
  obs_dt <- fread(obs_path)
  stopifnot("fc_correlation" %in% names(null_dt))
  stopifnot("fc_correlation" %in% names(obs_dt))

  obs_r <- as.numeric(obs_dt$fc_correlation[1L])
  null_vals <- as.numeric(null_dt$fc_correlation)
  z_val <- (obs_r - mean(null_vals, na.rm = TRUE)) / stats::sd(null_vals, na.rm = TRUE)

  p <- ggplot(data.frame(r = null_vals), aes(x = r)) +
    geom_histogram(bins = 50, fill = "grey70", color = "grey50", linewidth = 0.2) +
    geom_vline(xintercept = obs_r, color = "red", linewidth = 1) +
    annotate(
      "text",
      x = obs_r, y = Inf,
      label = sprintf("observed r = %.3f\nz = %.2f", obs_r, z_val),
      hjust = -0.05, vjust = 1.2, size = 3.5, color = "red", fontface = "bold"
    ) +
    labs(
      x = "Permuted FC correlation (intersection genes)",
      y = "Count",
      title = sprintf("Permutation null (%s)", task)
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  save_plot_both(p, sprintf("permutation_null_%s", task), width = 5.2, height = 4.2)
}

plot_stratified_fc <- function(task_name, bar_color = "#4878CF") {
  path <- file.path(REPO, "reports", "benchmark_master", "diagnostics", sprintf("fc_stratified_%s.csv", task_name))
  if (!file.exists(path)) return(NULL)

  dt <- fread(path)[method == "raw" & task == task_name]
  if (nrow(dt) == 0L) return(NULL)

  dt <- dt[is.finite(fc_correlation)]
  if (nrow(dt) == 0L) return(NULL)

  dt[, stratum := stratum_label(stratum)]
  dt[, label := sprintf("%s\n(n=%s)", as.character(stratum), formatC(n_genes, format = "f", digits = 0, big.mark = ","))]
  dt[, label := factor(label, levels = rev(unique(label)))]

  comp <- fread(file.path(BENCH, "comparison_summary.csv"))
  ceil <- comp[method == "raw" & task == task_name][["concordance_ceiling_fc_corr"]][1L]
  if (!is.finite(ceil)) ceil <- NA_real_

  ggplot(dt, aes(x = label, y = fc_correlation)) +
    geom_col(fill = bar_color, width = 0.65) +
    { if (is.finite(ceil)) geom_hline(yintercept = ceil, linetype = "dashed", color = "red", linewidth = 0.5) else NULL } +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    { if (is.finite(ceil)) annotate(
      "text",
      x = Inf, y = Inf,
      label = sprintf("ceiling = %.3f", ceil),
      hjust = 1.05, vjust = 1.2, size = 3, color = "red"
    ) else NULL } +
    coord_flip() +
    labs(x = NULL, y = "FC correlation (CPTAC vs CCLE)", title = sprintf("Stratified FC (%s)", task_name)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

figure_stratified_fc_combined <- function() {
  p_sub <- plot_stratified_fc("breast_subtype", bar_color = "#4878CF")
  p_bvl <- plot_stratified_fc("breast_vs_lung", bar_color = "#D65F5F")
  if (is.null(p_sub) && is.null(p_bvl)) return(invisible(NULL))
  if (is.null(p_sub)) {
    save_plot_both(p_bvl, "stratified_fc_combined", width = 6.8, height = 3.9)
    return(invisible(NULL))
  }
  if (is.null(p_bvl)) {
    save_plot_both(p_sub, "stratified_fc_combined", width = 6.8, height = 3.9)
    return(invisible(NULL))
  }

  final <- plot_grid(p_sub, p_bvl, nrow = 1L, rel_widths = c(1.25, 1.0), align = "h")
  save_plot_both(final, "stratified_fc_combined", width = 10.2, height = 3.9)
}

figure_destruction_overlap_subtype <- function() {
  task <- "breast_subtype"
  raw_c <- load_da("raw", task, "cptac")
  cell_c <- load_da("celligner", task, "cptac")
  if (is.null(raw_c) || is.null(cell_c)) return(invisible(NULL))

  gc <- gene_col_from_dt(raw_c)
  ge <- gene_col_from_dt(cell_c)
  stopifnot("adj.P.Val" %in% names(raw_c), "adj.P.Val" %in% names(cell_c))

  r <- unique(raw_c[, .(gene = as.character(get(gc)), sig = adj.P.Val < 0.05)])
  c <- unique(cell_c[, .(gene = as.character(get(ge)), sig = adj.P.Val < 0.05)])
  m <- merge(r, c, by = "gene", suffixes = c("_raw", "_cell"))

  retained <- sum(m$sig_raw & m$sig_cell)
  lost_gt <- sum(m$sig_raw & !m$sig_cell)
  new_cell <- sum(!m$sig_raw & m$sig_cell)

  df <- data.table(
    category = c(
      sprintf("Retained\n(%d)", retained),
      sprintf("Lost from\nraw CPTAC sig.\n(%d)", lost_gt),
      sprintf("New in\nCelligner sig.\n(%d)", new_cell)
    ),
    count = c(retained, lost_gt, new_cell)
  )
  df[, category := factor(category, levels = category)]

  p <- ggplot(df, aes(x = category, y = count)) +
    geom_col(fill = "#4878CF", width = 0.65) +
    geom_text(aes(label = count), vjust = -0.25, size = 4, fontface = "bold") +
    labs(
      x = NULL,
      y = "Gene count",
      title = "Subtype: overlap of significant CPTAC genes (raw vs Celligner)",
      subtitle = sprintf(
        "Raw CPTAC sig: %d | Celligner CPTAC sig: %d",
        sum(r$sig), sum(c$sig)
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  save_plot_both(p, "destruction_overlap", width = 7.2, height = 4.6)
}

figure_disconnect_scatter <- function() {
  disc_path <- file.path(BENCH, "disconnect_scores.csv")
  if (!file.exists(disc_path)) return(invisible(NULL))

  d <- fread(disc_path)
  # expected columns from benchmark_master/disconnect_scores.csv
  need <- c("method", "task", "geom_improvement", "da_improvement", "biology_cost")
  stopifnot(all(need %in% names(d)))

  d <- d[method %in% c("bridge_shift", "bridge_scale", "celligner")]
  d[, method_lab := method_display(method)]
  d[, task_lab := fifelse(task == "breast_subtype", "Subtype", "BvL")]

  p <- ggplot(d, aes(x = geom_improvement, y = da_improvement)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
    geom_point(aes(color = task_lab, shape = method_lab, size = pmax(biology_cost, 1e-6)), alpha = 0.85)

  if (isTRUE(HAS_GGREP)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = paste0(method_lab, " / ", task_lab)),
      size = 3.0,
      box.padding = 0.35,
      max.overlaps = 50,
      show.legend = FALSE
    )
  } else {
    p <- p + geom_text(
      aes(label = paste0(method_lab, " / ", task_lab)),
      size = 2.8,
      vjust = -1.0,
      show.legend = FALSE
    )
  }

  p <- p +
    scale_color_manual(values = c(Subtype = "#4878CF", BvL = "#D62728"), name = "Task") +
    scale_shape_manual(
      values = c(`Bridge shift` = 16, `Bridge shift+scale` = 15, Celligner = 17),
      breaks = c("Bridge shift", "Bridge shift+scale", "Celligner"),
      name = "Method"
    ) +
    scale_size_continuous(range = c(2.5, 10), name = "Biology cost") +
    labs(
      x = "Geometry improvement",
      y = "DA improvement",
      title = "Disconnect: geometry vs DA (non-raw methods)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  save_plot_both(p, "disconnect_scatter", width = 7.4, height = 5.2)
}

figure_geometry_bars_bvl <- function() {
  comp <- fread(file.path(BENCH, "comparison_summary.csv"))
  d <- comp[task == "breast_vs_lung" & method %in% c("raw", "bridge_shift", "bridge_scale", "celligner")]
  if (nrow(d) == 0L) return(invisible(NULL))

  rows <- rbindlist(list(
    d[, .(metric = "Domain R²", value = struct_domain_r2_pc1, method)],
    d[, .(metric = "Domain silh.", value = struct_silhouette_domain, method)],
    d[, .(metric = "Condition R²", value = struct_condition_r2_pc1, method)],
    d[, .(metric = "Condition silh.", value = struct_silhouette_condition, method)]
  ))
  rows[, method := method_display(method)]
  rows[, method := factor(method, levels = unname(METHOD_LABELS[METHOD_IDS]))]
  rows[, metric := factor(metric, levels = c("Domain R²", "Domain silh.", "Condition R²", "Condition silh."))]

  p <- ggplot(rows, aes(x = method, y = value, fill = method)) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    facet_wrap(~metric, scales = "free_y", nrow = 1L) +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(x = NULL, y = "Value", title = "Geometric metrics: breast vs lung") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  save_plot_both(p, "geometry_bars_bvl", width = 10.5, height = 4.2)
}

figure_marker_lollipop_subtype <- function() {
  task <- "breast_subtype"
  markers <- c("FOXA1", "EGFR", "KRT5", "KRT17")
  raw_c <- load_da("raw", task, "cptac")
  raw_e <- load_da("raw", task, "ccle")
  cell_c <- load_da("celligner", task, "cptac")
  cell_e <- load_da("celligner", task, "ccle")
  if (is.null(raw_c) || is.null(raw_e) || is.null(cell_c) || is.null(cell_e)) return(invisible(NULL))

  gc <- gene_col_from_dt(raw_c)
  ge <- gene_col_from_dt(raw_e)

  pick_fc <- function(dt, gcol, genes) {
    dt <- dt[as.character(get(gcol)) %in% genes, .(gene = as.character(get(gcol)), logFC)]
    unique(dt, by = "gene")
  }

  rc <- pick_fc(raw_c, gc, markers)
  re <- pick_fc(raw_e, ge, markers)
  cc <- pick_fc(cell_c, gc, markers)
  ce <- pick_fc(cell_e, ge, markers)

  cpt <- merge(rc, cc, by = "gene", suffixes = c("_raw", "_cell"))
  ccle <- merge(re, ce, by = "gene", suffixes = c("_raw", "_cell"))

  df <- rbind(
    cpt[, .(gene, domain = "CPTAC", raw_fc = logFC_raw, cell_fc = logFC_cell)],
    ccle[, .(gene, domain = "CCLE", raw_fc = logFC_raw, cell_fc = logFC_cell)],
    fill = TRUE
  )

  p <- ggplot(df, aes(y = gene)) +
    geom_segment(aes(x = raw_fc, xend = cell_fc, yend = gene), color = "grey55", linewidth = 0.8) +
    geom_point(aes(x = raw_fc), color = "#4878CF", size = 3) +
    geom_point(aes(x = cell_fc), color = "#D62728", size = 3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    facet_wrap(~domain, scales = "free_x", nrow = 1L) +
    labs(
      x = "logFC (limma; subtype contrast as in DA tables)",
      y = NULL,
      title = "Marker FC: raw vs Celligner",
      subtitle = "Blue = raw, Red = Celligner; segment shows change"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  save_plot_both(p, "marker_lollipop", width = 8.6, height = 4.2)
}

for (task in c("breast_subtype", "breast_vs_lung")) {
  figure_pca_panel(task)
  figure_fc_scatter_panel(task)
  figure_se_comparison(task)
  figure_permutation_null(task)
}

figure_stratified_fc_combined()
figure_destruction_overlap_subtype()
figure_disconnect_scatter()
figure_geometry_bars_bvl()
figure_marker_lollipop_subtype()

message("Wrote report figures to: ", OUT)
