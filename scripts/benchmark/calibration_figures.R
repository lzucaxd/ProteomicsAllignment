#!/usr/bin/env Rscript
# =============================================================================
# Calibration Analysis Figures (ggplot2)
# =============================================================================
# Generates diagnostic and analysis figures from calibration outputs:
#   1. FC scatterplot with calibration context (ceiling + null envelope)
#   2. Permutation null histograms (observed vs null per metric)
#   3. Biology destruction bar chart (retention + shrinkage across methods)
#   4. Residual dependence heatmap (sample x sample correlation)
#
# Usage: Rscript calibration_figures.R --repo-root /path/to/repo
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ── Parse arguments ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

results_dir <- file.path(repo_root, "reports/benchmark_master/benchmark_results")
fig_outdir  <- file.path(repo_root, "reports/benchmark_master/calibration/figures")
dir.create(fig_outdir, recursive = TRUE, showWarnings = FALSE)

theme_bench <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# =============================================================================
# 1. FC Scatterplot with Calibration Context
# =============================================================================
plot_fc_scatter_calibrated <- function(results_dir, fig_outdir) {
  method_dirs <- list.dirs(results_dir, recursive = FALSE)

  for (mdir in method_dirs) {
    method <- basename(mdir)
    task_dirs <- list.dirs(mdir, recursive = FALSE)

    for (tdir in task_dirs) {
      task <- basename(tdir)
      fc_file <- file.path(tdir, "representation_da", "fc_agreement.csv")
      if (!file.exists(fc_file)) next

      fc <- fread(fc_file)
      if (nrow(fc) < 5 || !all(c("logFC_cptac", "logFC_ccle") %in% names(fc))) next

      # Read calibration context if available
      ceil_file <- file.path(tdir, "calibration", "ceiling_summary.csv")
      null_file <- file.path(tdir, "calibration", "observed_vs_null_summary.csv")

      ceil_r <- NA_real_
      null_mean <- NA_real_
      null_sd <- NA_real_

      if (file.exists(ceil_file)) {
        ceil <- fread(ceil_file)
        if ("ceiling_fc_correlation" %in% names(ceil))
          ceil_r <- mean(ceil$ceiling_fc_correlation, na.rm = TRUE)
      }
      if (file.exists(null_file)) {
        null_dt <- fread(null_file)
        fc_row <- null_dt[metric == "fc_correlation"]
        if (nrow(fc_row) > 0) {
          null_mean <- fc_row$null_mean[1]
          null_sd <- fc_row$null_sd[1]
        }
      }

      obs_r <- cor(fc$logFC_cptac, fc$logFC_ccle, use = "pairwise.complete.obs")
      lim <- max(abs(c(fc$logFC_cptac, fc$logFC_ccle)), na.rm = TRUE) * 1.05

      p <- ggplot(fc, aes(x = logFC_cptac, y = logFC_ccle)) +
        geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "gray70", linewidth = 0.3) +
        geom_abline(slope = 1, intercept = 0, color = "gray50",
                     linetype = "dashed", linewidth = 0.4) +
        geom_point(alpha = 0.15, size = 0.8, color = "steelblue") +
        coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
        labs(
          title = paste0(method, " — ", task),
          subtitle = paste0("r = ", round(obs_r, 3),
                            if (!is.na(ceil_r)) paste0("  |  ceiling = ", round(ceil_r, 3)) else "",
                            if (!is.na(null_mean)) paste0("  |  null = ", round(null_mean, 3),
                                                           " \u00b1 ", round(null_sd, 3)) else ""),
          x = "logFC (CPTAC)", y = "logFC (CCLE)"
        ) +
        theme_bench

      # Annotate ceiling line if available
      if (!is.na(ceil_r)) {
        p <- p + geom_abline(slope = ceil_r, intercept = 0,
                              color = "darkgreen", linetype = "dotted", linewidth = 0.6)
      }

      fname <- paste0("fc_scatter_calibrated_", method, "_", task, ".png")
      ggsave(file.path(fig_outdir, fname), p, width = 5.5, height = 5.5, dpi = 200)
      cat("  Saved:", fname, "\n")
    }
  }
}

# =============================================================================
# 2. Permutation Null Histograms
# =============================================================================
plot_permutation_histograms <- function(results_dir, fig_outdir) {
  method_dirs <- list.dirs(results_dir, recursive = FALSE)
  all_nulls <- list()

  for (mdir in method_dirs) {
    method <- basename(mdir)
    task_dirs <- list.dirs(mdir, recursive = FALSE)

    for (tdir in task_dirs) {
      task <- basename(tdir)
      null_file <- file.path(tdir, "calibration", "null_distribution.csv")
      obs_file  <- file.path(tdir, "calibration", "observed_vs_null_summary.csv")
      if (!file.exists(null_file) || !file.exists(obs_file)) next

      null_dt <- fread(null_file)
      obs_dt  <- fread(obs_file)

      for (metric in c("fc_correlation", "same_direction_frac", "marker_concordance")) {
        if (!metric %in% names(null_dt)) next
        null_vals <- null_dt[[metric]]
        null_vals <- null_vals[is.finite(null_vals)]
        if (length(null_vals) < 10) next

        obs_row <- obs_dt[metric == (metric)]
        obs_val <- if (nrow(obs_row) > 0) obs_row$observed[1] else NA_real_

        plot_dt <- data.table(value = null_vals)

        p <- ggplot(plot_dt, aes(x = value)) +
          geom_histogram(bins = 30, fill = "gray70", color = "gray50", alpha = 0.7) +
          {if (!is.na(obs_val)) geom_vline(xintercept = obs_val, color = "red",
                                            linewidth = 1, linetype = "dashed")} +
          labs(
            title = paste0("Permutation Null: ", metric),
            subtitle = paste0(method, " / ", task,
                              if (!is.na(obs_val)) paste0("  |  observed = ", round(obs_val, 3)) else ""),
            x = metric, y = "Count"
          ) +
          theme_bench

        fname <- paste0("perm_null_", metric, "_", method, "_", task, ".png")
        ggsave(file.path(fig_outdir, fname), p, width = 5, height = 3.5, dpi = 200)
      }
    }
  }
  cat("  Permutation null histograms saved\n")
}

# =============================================================================
# 3. Biology Destruction Bar Chart
# =============================================================================
plot_biology_destruction <- function(results_dir, fig_outdir) {
  method_dirs <- list.dirs(results_dir, recursive = FALSE)
  all_summaries <- list()

  for (mdir in method_dirs) {
    method <- basename(mdir)
    task_dirs <- list.dirs(mdir, recursive = FALSE)

    for (tdir in task_dirs) {
      task <- basename(tdir)
      files <- list.files(file.path(tdir, "calibration"),
                          pattern = "^destruction_summary_", full.names = TRUE)
      for (f in files) {
        dt <- fread(f)
        dt$task <- task
        all_summaries[[length(all_summaries) + 1]] <- dt
      }
    }
  }

  if (length(all_summaries) == 0) {
    cat("  No biology destruction data found\n")
    return(invisible())
  }

  combined <- rbindlist(all_summaries, fill = TRUE)
  if (nrow(combined) == 0) return(invisible())

  # Retention rate
  if ("default_retention_rate" %in% names(combined)) {
    p1 <- ggplot(combined, aes(x = method, y = default_retention_rate, fill = task)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      labs(title = "Biology Destruction: Gene Retention Rate",
           subtitle = "Fraction of natively significant genes still significant after harmonization",
           x = NULL, y = "Retention Rate") +
      theme_bench +
      scale_fill_brewer(palette = "Set2")

    ggsave(file.path(fig_outdir, "biology_destruction_retention.png"),
           p1, width = 7, height = 4.5, dpi = 200)
  }

  # FC shrinkage
  if ("default_mean_fc_shrinkage" %in% names(combined)) {
    p2 <- ggplot(combined, aes(x = method, y = default_mean_fc_shrinkage, fill = task)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      labs(title = "Biology Destruction: FC Shrinkage",
           subtitle = "Mean fold-change shrinkage among natively significant genes",
           x = NULL, y = "Mean FC Shrinkage") +
      theme_bench +
      scale_fill_brewer(palette = "Set2")

    ggsave(file.path(fig_outdir, "biology_destruction_shrinkage.png"),
           p2, width = 7, height = 4.5, dpi = 200)
  }

  cat("  Biology destruction charts saved\n")
}

# =============================================================================
# 4. Residual Dependence Heatmap
# =============================================================================
plot_residual_heatmaps <- function(results_dir, fig_outdir) {
  method_dirs <- list.dirs(results_dir, recursive = FALSE)

  for (mdir in method_dirs) {
    method <- basename(mdir)
    task_dirs <- list.dirs(mdir, recursive = FALSE)

    for (tdir in task_dirs) {
      task <- basename(tdir)
      corr_files <- list.files(file.path(tdir, "calibration"),
                                pattern = "^residual_corr_matrix_", full.names = TRUE)

      for (cf in corr_files) {
        domain <- gsub("residual_corr_matrix_|\\.csv", "", basename(cf))
        dt <- fread(cf)

        if (ncol(dt) < 3) next

        rn <- dt[[1]]
        mat <- as.matrix(dt[, -1, with = FALSE])
        rownames(mat) <- rn

        # Melt for ggplot
        melt_dt <- data.table(
          row = rep(rn, each = ncol(mat)),
          col = rep(colnames(mat), times = nrow(mat)),
          value = as.vector(mat)
        )
        melt_dt[, row := factor(row, levels = rev(rn))]
        melt_dt[, col := factor(col, levels = colnames(mat))]

        p <- ggplot(melt_dt, aes(x = col, y = row, fill = value)) +
          geom_tile() +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                midpoint = 0, limits = c(-1, 1),
                                name = "Residual\nCorrelation") +
          labs(
            title = paste0("Residual Dependence: ", method, " / ", task, " / ", domain),
            x = NULL, y = NULL
          ) +
          theme_bench +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
            axis.text.y = element_text(size = 5)
          )

        fname <- paste0("residual_heatmap_", method, "_", task, "_", domain, ".png")
        ggsave(file.path(fig_outdir, fname), p,
               width = max(6, ncol(mat) * 0.15 + 2),
               height = max(6, nrow(mat) * 0.15 + 2),
               dpi = 150)
      }
    }
  }
  cat("  Residual dependence heatmaps saved\n")
}

# =============================================================================
# Run all
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("  CALIBRATION ANALYSIS FIGURES\n")
cat(strrep("=", 60), "\n\n")

if (!dir.exists(results_dir)) {
  cat("No benchmark results directory found at:", results_dir, "\n")
  quit(status = 0)
}

cat("1. FC scatterplots with calibration context\n")
tryCatch(plot_fc_scatter_calibrated(results_dir, fig_outdir),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("2. Permutation null histograms\n")
tryCatch(plot_permutation_histograms(results_dir, fig_outdir),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("3. Biology destruction charts\n")
tryCatch(plot_biology_destruction(results_dir, fig_outdir),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("4. Residual dependence heatmaps\n")
tryCatch(plot_residual_heatmaps(results_dir, fig_outdir),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("\nAll calibration figures saved to:", fig_outdir, "\n")
