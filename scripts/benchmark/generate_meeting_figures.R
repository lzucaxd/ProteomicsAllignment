#!/usr/bin/env Rscript
# =============================================================================
# Meeting Slide Figures — presentation-ready PDF/PNG for key insights
# =============================================================================
# Generates polished figures for committee/advisor meetings:
#   1. Contrast validation table
#   2. Concordance ceiling context plot
#   3. Cross-method comparison table
#   4. FC scatterplot panels (tasks x methods)
#   5. Fixed-basis PCA panels
#   6. Marker profile reference (points to existing base-R plots)
#   7. Subtype diagnostic narrative
#   8. Roadmap figure
#
# Usage: Rscript generate_meeting_figures.R --repo-root /path/to/repo
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

# ── Parse arguments ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

results_dir <- file.path(repo_root, "reports/benchmark_master/benchmark_results")
diag_dir    <- file.path(repo_root, "reports/benchmark_master/diagnostics")
meeting_dir <- file.path(repo_root, "reports/benchmark_master/meeting/figures")
dir.create(meeting_dir, recursive = TRUE, showWarnings = FALSE)

theme_slide <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

manifest <- data.table(figure = character(), file = character(), description = character())

add_to_manifest <- function(name, file, desc) {
  manifest <<- rbind(manifest, data.table(figure = name, file = file, description = desc))
}

# =============================================================================
# 1. Contrast Validation Table
# =============================================================================
generate_contrast_validation_table <- function() {
  diag_file <- file.path(diag_dir, "subtype_sign_diagnostic.csv")
  if (!file.exists(diag_file)) {
    cat("  [1] No subtype sign diagnostic found\n")
    return(invisible())
  }

  dt <- fread(diag_file)
  if (nrow(dt) == 0) return(invisible())

  dt[, status := ifelse(correct, "OK", "WRONG")]
  dt[, color_val := ifelse(correct, 1, -1)]

  p <- ggplot(dt, aes(x = domain, y = gene, fill = factor(correct))) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%+.2f", logFC)), size = 3.5) +
    scale_fill_manual(
      values = c("TRUE" = "#b8e6b8", "FALSE" = "#f5b7b1"),
      labels = c("TRUE" = "Correct", "FALSE" = "Wrong"),
      name = "Direction"
    ) +
    labs(
      title = "Contrast Direction Validation",
      subtitle = "Marker logFC signs in raw data (Luminal - Basal)",
      x = NULL, y = NULL
    ) +
    theme_slide +
    theme(panel.grid = element_blank())

  fname <- "01_contrast_validation_table.png"
  ggsave(file.path(meeting_dir, fname), p, width = 5, height = 5.5, dpi = 250)
  ggsave(file.path(meeting_dir, gsub("png$", "pdf", fname)), p, width = 5, height = 5.5)
  add_to_manifest("Contrast Validation", fname, "Marker logFC sign check per domain on raw data")
  cat("  [1] Contrast validation table saved\n")
}

# =============================================================================
# 2. Concordance Ceiling Context Plot
# =============================================================================
generate_ceiling_context <- function() {
  method_dirs <- list.dirs(results_dir, recursive = FALSE)
  rows <- list()

  for (mdir in method_dirs) {
    method <- basename(mdir)
    task_dirs <- list.dirs(mdir, recursive = FALSE)

    for (tdir in task_dirs) {
      task <- basename(tdir)

      # Observed FC correlation
      fc_file <- file.path(tdir, "representation_da", "fc_agreement.csv")
      if (!file.exists(fc_file)) next
      fc <- fread(fc_file)
      if (!all(c("logFC_cptac", "logFC_ccle") %in% names(fc))) next
      obs_r <- cor(fc$logFC_cptac, fc$logFC_ccle, use = "pairwise.complete.obs")

      # Ceiling (CPTAC + optional CCLE split-half)
      ceil_cptac <- NA_real_
      cf <- file.path(tdir, "calibration", "ceiling_summary_cptac.csv")
      if (file.exists(cf)) {
        ce <- fread(cf)
        if ("ceiling_fc_correlation" %in% names(ce)) ceil_cptac <- ce$ceiling_fc_correlation[1]
      } else {
        leg <- file.path(tdir, "calibration", "ceiling_summary.csv")
        if (file.exists(leg)) {
          ce <- fread(leg)
          if ("ceiling_fc_correlation" %in% names(ce)) ceil_cptac <- ce$ceiling_fc_correlation[1]
        }
      }
      ceil_ccle <- NA_real_
      cf2 <- file.path(tdir, "calibration", "ceiling_summary_ccle.csv")
      if (file.exists(cf2)) {
        ce2 <- fread(cf2)
        if ("ceiling_fc_correlation" %in% names(ce2)) ceil_ccle <- ce2$ceiling_fc_correlation[1]
      }

      # Null
      null_file <- file.path(tdir, "calibration", "observed_vs_null_summary.csv")
      null_mean <- NA_real_
      null_sd <- NA_real_
      if (file.exists(null_file)) {
        null_dt <- fread(null_file)
        fc_row <- null_dt[metric == "fc_correlation"]
        if (nrow(fc_row) > 0) {
          null_mean <- fc_row$null_mean[1]
          null_sd <- fc_row$null_sd[1]
        }
      }

      rows[[length(rows) + 1]] <- data.table(
        method = method, task = task,
        observed = obs_r, ceiling = ceil_cptac, ceiling_ccle = ceil_ccle,
        null_mean = null_mean, null_sd = null_sd
      )
    }
  }

  if (length(rows) == 0) {
    cat("  [2] No data for ceiling context plot\n")
    return(invisible())
  }

  dt <- rbindlist(rows)

  p <- ggplot(dt, aes(x = method, y = observed, color = task)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = null_mean - 2 * null_sd, ymax = null_mean + 2 * null_sd),
                   width = 0.15, linetype = "dotted", alpha = 0.6) +
    geom_point(aes(y = null_mean), shape = 4, size = 3, alpha = 0.6) +
    geom_point(aes(y = ceiling), shape = 17, size = 3, alpha = 0.8) +
    geom_point(data = dt[is.finite(ceiling_ccle)],
               aes(x = method, y = ceiling_ccle, color = task),
               shape = 15, size = 3, alpha = 0.75, position = position_nudge(x = 0.12)) +
    geom_hline(yintercept = 0, color = "gray60", linetype = "dashed") +
    facet_wrap(~ task, scales = "free_x") +
    labs(
      title = "Cross-Domain FC Correlation: Observed vs Ceiling vs Null",
      subtitle = "Triangle = CPTAC ceiling | Square (nudged) = CCLE split-half ceiling | X = null | Circle = observed",
      x = NULL, y = "FC Correlation (Pearson r)"
    ) +
    theme_slide +
    scale_color_brewer(palette = "Set1")

  fname <- "02_concordance_ceiling_context.png"
  ggsave(file.path(meeting_dir, fname), p, width = 9, height = 5, dpi = 250)
  ggsave(file.path(meeting_dir, gsub("png$", "pdf", fname)), p, width = 9, height = 5)
  add_to_manifest("Concordance Ceiling Context", fname,
                  "Observed FC correlation vs within-domain ceiling and permutation null")
  cat("  [2] Concordance ceiling context plot saved\n")
}

# =============================================================================
# 3. Comparison Table (ggplot2 heatmap-style)
# =============================================================================
generate_comparison_table <- function() {
  csv_file <- file.path(results_dir, "comparison_summary.csv")
  if (!file.exists(csv_file)) {
    cat("  [3] No comparison summary found\n")
    return(invisible())
  }

  dt <- fread(csv_file)
  if (nrow(dt) == 0) return(invisible())

  # Select key metrics to display
  key_metrics <- c("fc_correlation", "fc_same_dir_frac",
                    "struct_domain_r2_pc1", "struct_condition_r2_pc1",
                    "marker_sanity_cptac", "marker_sanity_ccle",
                    "concordance_ceiling_fc_corr", "calibrated_fc_corr",
                    "biology_destruction_retention")

  # Also check column name variants from the run
  da_variants <- c("da_fc_correlation", "da_same_direction_frac")
  available <- intersect(c(key_metrics, da_variants), names(dt))

  if (length(available) == 0) {
    cat("  [3] No displayable metrics in comparison summary\n")
    return(invisible())
  }

  id_cols <- intersect(c("method", "task"), names(dt))
  melt_dt <- melt(dt[, c(id_cols, available), with = FALSE],
                   id.vars = id_cols, variable.name = "metric", value.name = "value")

  melt_dt[, value := as.numeric(value)]
  melt_dt <- melt_dt[is.finite(value)]

  p <- ggplot(melt_dt, aes(x = method, y = metric, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", value)), size = 3) +
    facet_wrap(~ task, ncol = 1) +
    scale_fill_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850",
                          midpoint = 0.5, limits = c(-0.5, 1),
                          oob = scales::squish, name = "Value") +
    labs(
      title = "Method Comparison Summary",
      x = NULL, y = NULL
    ) +
    theme_slide +
    theme(axis.text.y = element_text(size = 9),
          panel.grid = element_blank())

  n_metrics <- length(unique(melt_dt$metric))
  fig_h <- max(4, n_metrics * 0.6 + 2)

  fname <- "03_comparison_table.png"
  ggsave(file.path(meeting_dir, fname), p, width = 9, height = fig_h, dpi = 250)
  ggsave(file.path(meeting_dir, gsub("png$", "pdf", fname)), p, width = 9, height = fig_h)
  add_to_manifest("Comparison Table", fname, "Heatmap of key metrics across methods and tasks")
  cat("  [3] Comparison table saved\n")
}

# =============================================================================
# 4. FC Scatterplot Panels (tasks x methods)
# =============================================================================
generate_fc_scatter_panels <- function() {
  method_dirs <- list.dirs(results_dir, recursive = FALSE)
  all_fc <- list()

  for (mdir in method_dirs) {
    method <- basename(mdir)
    task_dirs <- list.dirs(mdir, recursive = FALSE)

    for (tdir in task_dirs) {
      task <- basename(tdir)
      fc_file <- file.path(tdir, "representation_da", "fc_agreement.csv")
      if (!file.exists(fc_file)) next
      fc <- fread(fc_file)
      if (!all(c("logFC_cptac", "logFC_ccle") %in% names(fc))) next
      fc$method <- method
      fc$task <- task
      all_fc[[length(all_fc) + 1]] <- fc
    }
  }

  if (length(all_fc) == 0) {
    cat("  [4] No FC data for scatter panels\n")
    return(invisible())
  }

  combined <- rbindlist(all_fc, fill = TRUE)

  p <- ggplot(combined, aes(x = logFC_cptac, y = logFC_ccle)) +
    geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "gray80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "gray80", linewidth = 0.3) +
    geom_point(alpha = 0.1, size = 0.5, color = "steelblue") +
    facet_grid(task ~ method) +
    coord_fixed() +
    labs(
      title = "Cross-Domain Fold Change Agreement",
      x = "logFC (CPTAC)", y = "logFC (CCLE)"
    ) +
    theme_slide

  n_methods <- length(unique(combined$method))
  n_tasks <- length(unique(combined$task))

  fname <- "04_fc_scatter_panels.png"
  ggsave(file.path(meeting_dir, fname), p,
         width = min(14, 3.5 * n_methods + 1), height = 3.5 * n_tasks + 1, dpi = 200)
  ggsave(file.path(meeting_dir, gsub("png$", "pdf", fname)), p,
         width = min(14, 3.5 * n_methods + 1), height = 3.5 * n_tasks + 1)
  add_to_manifest("FC Scatter Panels", fname,
                  "Tasks x Methods grid of cross-domain fold change scatter")
  cat("  [4] FC scatter panels saved\n")
}

# =============================================================================
# 5. Fixed-Basis PCA Panels
# =============================================================================
generate_pca_panels <- function() {
  # Collect existing PCA plot PNGs
  pca_files <- list.files(results_dir, pattern = "pca_domain_.*\\.png$",
                           recursive = TRUE, full.names = TRUE)
  if (length(pca_files) == 0) {
    cat("  [5] No PCA plots found\n")
    return(invisible())
  }

  # Just create a manifest entry pointing to the structure/ dirs
  # since PCA plots are already method-specific PNG files
  add_to_manifest("PCA Structure Plots", "See structure/ subdirectories",
                  "PCA plots colored by domain and condition for each method x task")
  cat("  [5] PCA panels noted in manifest (existing plots in structure/ dirs)\n")
}

# =============================================================================
# 6. Marker Profiles (reference to existing)
# =============================================================================
generate_marker_profile_reference <- function() {
  profile_dir <- file.path(repo_root, "reports/benchmark_master/marker_profiles")
  if (dir.exists(profile_dir)) {
    profile_files <- list.files(profile_dir, pattern = "\\.png$", recursive = TRUE)
    if (length(profile_files) > 0) {
      add_to_manifest("Marker Profiles", "See marker_profiles/ directory",
                      paste0(length(profile_files), " profile plot files in marker_profiles/"))
      cat("  [6] Marker profiles noted in manifest (", length(profile_files), " files)\n")
    }
  } else {
    cat("  [6] No marker profile plots found\n")
  }
}

# =============================================================================
# 7. Subtype Diagnostic / Limitation Narrative
# =============================================================================
generate_subtype_diagnostic <- function() {
  text_lines <- c(
    "Breast Subtype Benchmark: Diagnostic Summary",
    "",
    "FINDING: The initial benchmark showed negative FC correlation (-0.36) and",
    "below-chance same-direction fraction (39.6%) for breast subtype across ALL",
    "harmonization methods, including raw unharmonized data.",
    "",
    "ROOT CAUSE: The Python DA engine (Welch t-tests) used metadata row order",
    "to determine group_a vs group_b, causing opposite contrast directions in",
    "CPTAC vs CCLE when the metadata listed conditions in different orders.",
    "",
    "RESOLUTION: Migrated DA to R limma with explicit contrast matrices.",
    "Factor levels are sorted alphabetically (Basal < Luminal), ensuring",
    "consistent Luminal - Basal contrast direction across all domains.",
    "",
    "The contrast_validation.R diagnostic now runs BEFORE any benchmark step",
    "to verify marker logFC signs against expected biology.",
    "",
    "REMAINING LIMITATIONS:",
    "- CCLE breast subtype has very few samples (4 per group)",
    "- CPTAC subtype assignment depends on mRNA-derived PAM50 mapping",
    "- Concordance ceiling for CCLE will be estimated via jackknife"
  )

  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste(text_lines, collapse = "\n"),
             hjust = 0.5, vjust = 0.5, size = 3.5, family = "mono") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    labs(title = "Subtype Benchmark: Sign-Flip Diagnostic") +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

  fname <- "07_subtype_diagnostic.png"
  ggsave(file.path(meeting_dir, fname), p, width = 9, height = 6, dpi = 200)
  ggsave(file.path(meeting_dir, gsub("png$", "pdf", fname)), p, width = 9, height = 6)
  add_to_manifest("Subtype Diagnostic", fname,
                  "Narrative explaining the sign-flip issue and resolution")
  cat("  [7] Subtype diagnostic narrative saved\n")
}

# =============================================================================
# 8. Roadmap Figure
# =============================================================================
generate_roadmap <- function() {
  roadmap_dt <- data.table(
    phase = c("Phase 0", "Phase 1", "Phase 2", "Phase 3", "Phase 4",
              "Next: ComBat", "Next: MOBER/VAE", "Next: DIA"),
    item = c("Contrast Validation", "Migrate DA to limma",
             "Calibration Infrastructure", "Updated Output Format",
             "Meeting Figures",
             "Add ComBat/removeBatchEffect", "Add neural harmonization",
             "Extend to DIA-based data"),
    status = c("Done", "Done", "Done", "Done", "Done",
               "Planned", "Planned", "Planned"),
    order = 1:8
  )

  roadmap_dt[, status_color := ifelse(status == "Done", "Complete", "Planned")]

  p <- ggplot(roadmap_dt, aes(x = order, y = 1, fill = status_color)) +
    geom_tile(width = 0.9, height = 0.6, color = "white", linewidth = 1) +
    geom_text(aes(label = paste0(phase, "\n", item)), size = 3, lineheight = 1.1) +
    scale_fill_manual(values = c("Complete" = "#b8e6b8", "Planned" = "#fce4b8"),
                       name = "Status") +
    labs(title = "Benchmark Infrastructure Roadmap") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "bottom"
    )

  fname <- "08_roadmap.png"
  ggsave(file.path(meeting_dir, fname), p, width = 12, height = 3, dpi = 200)
  ggsave(file.path(meeting_dir, gsub("png$", "pdf", fname)), p, width = 12, height = 3)
  add_to_manifest("Roadmap", fname, "Infrastructure roadmap: done and planned phases")
  cat("  [8] Roadmap figure saved\n")
}

# =============================================================================
# 9–12: Union vs intersection, stratified FC, volcano side-by-side, disconnect
# =============================================================================
generate_union_vs_intersection_bars <- function() {
  method_dirs <- list.dirs(results_dir, recursive = FALSE)
  rows <- list()
  for (mdir in method_dirs) {
    method <- basename(mdir)
    for (tdir in list.dirs(mdir, recursive = FALSE)) {
      task <- basename(tdir)
      cf <- file.path(tdir, "representation_da", "cross_domain_metrics.csv")
      if (!file.exists(cf)) next
      dt <- fread(cf)
      dt[, method := method][, task := task]
      rows[[length(rows) + 1]] <- dt
    }
  }
  if (length(rows) == 0) {
    cat("  [9] No cross_domain_metrics for union/intersection plot\n")
    return(invisible())
  }
  long <- rbindlist(rows)
  long <- long[gene_set %in% c("union", "intersection")]
  long[, gene_set := factor(gene_set, levels = c("union", "intersection"))]
  for (tk in unique(long$task)) {
    sub <- long[task == tk]
    if (nrow(sub) == 0) next
    p <- ggplot(sub, aes(x = method, y = fc_correlation, fill = gene_set)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      labs(
        title = paste("Cross-domain FC correlation:", tk),
        subtitle = "Lighter = union gene set; darker = intersection (both-domain coverage)",
        x = NULL, y = "Pearson r (logFC CPTAC vs CCLE)", fill = NULL
      ) +
      theme_slide +
      scale_fill_manual(values = c(union = "#9ecae1", intersection = "#08519c"))
    fname <- paste0("09_union_vs_intersection_", tk, ".pdf")
    ggsave(file.path(meeting_dir, fname), p, width = 8, height = 5)
    add_to_manifest(paste("Union vs Intersection", tk), fname, "FC correlation on union vs intersection genes")
  }
  cat("  [9] Union vs intersection bar charts saved\n")
}

generate_stratified_fc_bars <- function() {
  for (task in c("breast_subtype", "breast_vs_lung")) {
    cf <- file.path(diag_dir, paste0("fc_stratified_", task, ".csv"))
    if (!file.exists(cf)) next
    dt <- fread(cf)
    p <- ggplot(dt, aes(x = method, y = fc_correlation, fill = stratum)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      labs(
        title = paste("Stratified FC correlation:", task),
        x = NULL, y = "Pearson r", fill = "Significance stratum"
      ) +
      theme_slide +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
    fname <- paste0("10_fc_stratified_", task, ".pdf")
    ggsave(file.path(meeting_dir, fname), p, width = 9, height = 5)
    add_to_manifest(paste("Stratified FC", task), fname, "FC r by DA significance stratum")
  }
  cat("  [10] Stratified FC figures saved (if diagnostics exist)\n")
}

generate_volcano_comparison_task <- function(task) {
  da_c <- file.path(results_dir, "raw", task, "representation_da", "cptac", "da_limma_result.csv")
  da_e <- file.path(results_dir, "raw", task, "representation_da", "ccle", "da_limma_result.csv")
  if (!file.exists(da_c) || !file.exists(da_e)) return(NULL)
  plot_one <- function(path, title) {
    da <- fread(path)
    n_sig <- sum(da$adj.P.Val < 0.05, na.rm = TRUE)
    ggplot(da, aes(x = logFC, y = -log10(P.Value))) +
      geom_point(alpha = 0.12, size = 0.4, color = "grey45") +
      geom_point(data = da[da$adj.P.Val < 0.05], color = "steelblue", alpha = 0.35, size = 0.6) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      labs(title = title, subtitle = paste0(n_sig, " / ", nrow(da), " adj.p < 0.05")) +
      theme_minimal(base_size = 11)
  }
  p <- grid.arrange(plot_one(da_c, "CPTAC"), plot_one(da_e, "CCLE"), ncol = 2)
  fname <- paste0("11_volcano_", task, ".pdf")
  ggsave(file.path(meeting_dir, fname), p, width = 10, height = 4.5)
  add_to_manifest(paste("Volcano", task), fname, "Raw CPTAC vs CCLE volcanos side by side")
}

generate_volcano_comparisons <- function() {
  for (task in c("breast_subtype", "breast_vs_lung")) {
    tryCatch(generate_volcano_comparison_task(task),
             error = function(e) cat("  [11] volcano", task, ":", conditionMessage(e), "\n"))
  }
  cat("  [11] Volcano comparison panels done\n")
}

generate_disconnect_combined <- function() {
  csv_file <- file.path(results_dir, "comparison_summary.csv")
  if (!file.exists(csv_file)) {
    cat("  [12] No comparison_summary for disconnect plot\n")
    return(invisible())
  }
  dt <- fread(csv_file)
  need <- c("method", "task", "struct_domain_r2_pc1", "fc_correlation_intersection",
            "concordance_ceiling_fc_corr", "biology_destruction_retention")
  if (!all(need %in% names(dt))) {
    cat("  [12] comparison_summary missing columns for disconnect plot\n")
    return(invisible())
  }
  rows <- list()
  for (tk in unique(dt$task)) {
    raw_r <- dt[method == "raw" & task == tk]
    if (nrow(raw_r) != 1) next
    r_dr2 <- raw_r$struct_domain_r2_pc1[1]
    r_fc <- raw_r$fc_correlation_intersection[1]
    ceil <- raw_r$concordance_ceiling_fc_corr[1]
    denom <- ceil - r_fc
    for (m in unique(dt$method)) {
      if (m == "raw") next
      mr <- dt[method == m & task == tk]
      if (nrow(mr) != 1) next
      m_dr2 <- mr$struct_domain_r2_pc1[1]
      m_fc <- mr$fc_correlation_intersection[1]
      ret <- mr$biology_destruction_retention[1]
      geom_i <- if (is.finite(r_dr2) && r_dr2 != 0) (r_dr2 - m_dr2) / r_dr2 else NA_real_
      da_i <- if (is.finite(denom) && abs(denom) > 1e-8) (m_fc - r_fc) / denom else NA_real_
      bio_cost <- if (is.finite(ret)) 1 - ret else NA_real_
      rows[[length(rows) + 1]] <- data.table(
        method = m, task = tk,
        geom_improvement = geom_i, da_improvement = da_i,
        disconnect = geom_i - da_i, biology_cost = bio_cost
      )
    }
  }
  if (length(rows) == 0) return(invisible())
  plot_dt <- rbindlist(rows)
  p <- ggplot(plot_dt, aes(x = geom_improvement, y = da_improvement)) +
    geom_abline(slope = 1, intercept = 0, color = "grey55", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_point(aes(size = pmax(biology_cost, 0.01, na.rm = TRUE), color = task), alpha = 0.85) +
    scale_size_continuous(range = c(2, 8), name = "Biology cost\n(1 - retention)") +
    labs(
      title = "Geometry vs DA improvement (intersection FC)",
      subtitle = "Diagonal: matched gains; above diagonal: DA lags geometry",
      x = "Geometry improvement (domain R² PC1)", y = "DA improvement (normalized)"
    ) +
    theme_slide +
    facet_wrap(~ task)
  fname <- "12_disconnect_combined.pdf"
  ggsave(file.path(meeting_dir, fname), p, width = 9, height = 5)
  add_to_manifest("Disconnect combined", fname, "Geometry vs intersection-DA improvement")
  cat("  [12] Disconnect scatter saved\n")
}

generate_ccle_subtype_power_comparison <- function() {
  leg <- file.path(diag_dir, "legacy_ccle_da_limma_breast_subtype.csv")
  newp <- file.path(results_dir, "raw", "breast_subtype", "representation_da", "ccle", "da_limma_result.csv")
  if (!file.exists(newp)) {
    cat("  [13] Skip CCLE power: missing new CCLE DA\n")
    return(invisible())
  }
  volc <- function(path, stitle) {
    da <- fread(path)
    n_sig <- sum(da$adj.P.Val < 0.05, na.rm = TRUE)
    ggplot(da, aes(x = logFC, y = -log10(P.Value))) +
      geom_point(alpha = 0.12, size = 0.4, color = "grey45") +
      geom_point(data = da[da$adj.P.Val < 0.05], color = "steelblue", alpha = 0.35, size = 0.6) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      labs(title = stitle, subtitle = paste0(n_sig, " / ", nrow(da), " adj.p < 0.05")) +
      theme_minimal(base_size = 11)
  }
  fname <- "13_ccle_power_comparison_subtype.pdf"
  if (file.exists(leg)) {
    p <- grid.arrange(volc(leg, "CCLE (archived small-n)"), volc(newp, "CCLE (v2 expanded)"), ncol = 2)
  } else {
    p <- volc(newp, "CCLE subtype (v2 expanded lines)")
  }
  ggsave(file.path(meeting_dir, fname), p, width = 11, height = 4.5)
  add_to_manifest("CCLE power comparison", fname, "Optional legacy vs v2 CCLE subtype volcano")
  cat("  [13] CCLE power comparison saved\n")
}

# =============================================================================
# Run all
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("  MEETING SLIDE FIGURES\n")
cat(strrep("=", 60), "\n\n")

cat("[1] Contrast Validation Table\n")
tryCatch(generate_contrast_validation_table(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[2] Concordance Ceiling Context\n")
tryCatch(generate_ceiling_context(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[3] Comparison Table\n")
tryCatch(generate_comparison_table(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[4] FC Scatter Panels\n")
tryCatch(generate_fc_scatter_panels(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[5] PCA Panels\n")
tryCatch(generate_pca_panels(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[6] Marker Profiles\n")
tryCatch(generate_marker_profile_reference(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[7] Subtype Diagnostic\n")
tryCatch(generate_subtype_diagnostic(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[8] Roadmap\n")
tryCatch(generate_roadmap(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[9] Union vs intersection bars\n")
tryCatch(generate_union_vs_intersection_bars(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[10] Stratified FC bars\n")
tryCatch(generate_stratified_fc_bars(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[11] Volcano comparisons\n")
tryCatch(generate_volcano_comparisons(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[12] Disconnect scatter\n")
tryCatch(generate_disconnect_combined(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

cat("[13] CCLE subtype power (v2 vs optional legacy)\n")
tryCatch(generate_ccle_subtype_power_comparison(),
         error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))

# ── Save manifest ────────────────────────────────────────────────────────────
fwrite(manifest, file.path(meeting_dir, "figure_manifest.csv"))
cat("\nManifest saved:", file.path(meeting_dir, "figure_manifest.csv"), "\n")

# ── Generate HTML index ──────────────────────────────────────────────────────
html_lines <- c(
  "<!DOCTYPE html>",
  "<html><head><title>Meeting Figures</title>",
  "<style>body{font-family:sans-serif;max-width:1200px;margin:auto;padding:20px}",
  "h1{color:#333}h2{color:#555;margin-top:30px}",
  "img{max-width:100%;border:1px solid #ddd;margin:10px 0}",
  ".desc{color:#666;font-style:italic}</style></head><body>",
  "<h1>Benchmark Meeting Figures</h1>",
  paste0("<p>Generated: ", Sys.time(), "</p>")
)

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  html_lines <- c(html_lines,
    paste0("<h2>", row$figure, "</h2>"),
    paste0("<p class='desc'>", row$description, "</p>"),
    if (grepl("\\.(png|pdf)$", row$file))
      paste0("<img src='", row$file, "' alt='", row$figure, "'>")
    else
      paste0("<p>", row$file, "</p>")
  )
}

html_lines <- c(html_lines, "</body></html>")
writeLines(html_lines, file.path(meeting_dir, "index.html"))
cat("HTML index saved:", file.path(meeting_dir, "index.html"), "\n")

cat("\nDone.\n")
