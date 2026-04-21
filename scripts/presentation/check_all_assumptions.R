#!/usr/bin/env Rscript
# For every method × domain × task where limma DA exists (benchmark_results),
# generate diagnostic plots and a quantitative summary table.
#
# Output:
#   presentation_materials/figures/assumptions/
#     qq_residuals_{domain}_{task}.pdf
#     resid_fitted_{domain}_{task}.pdf
#     sa_{domain}_{task}.pdf
#   presentation_materials/tables/assumption_summary.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) {
  dirname(normalizePath(sub("^--file=", "", ff[1L])))
} else {
  normalizePath(file.path(getwd(), "scripts", "presentation"), mustWork = FALSE)
}
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

OUT_FIG <- file.path(PRES_OUT, "figures", "assumptions")
OUT_TAB <- file.path(PRES_OUT, "tables")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
tasks <- c("breast_subtype", "breast_vs_lung")
domains <- c("cptac", "ccle")
markers_subtype <- c("FOXA1", "GATA3", "EGFR", "KRT5")
markers_bvl <- c("NKX2-1", "SFTPB", "NAPSA", "GATA3")

load_mat <- function(path) {
  dt <- fread(path, check.names = FALSE)
  gn <- names(dt)[1L]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[gn]]
  colnames(mat) <- names(dt)[-1L]
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  mat
}

empty_panel <- function() {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "not quantified", size = 3, color = "grey50") +
    theme_void()
}

# QQ reference line through quartiles (y = sorted residuals vs theoretical N(0,1) quantiles)
qq_line_ab <- function(y_sorted) {
  probs <- c(0.25, 0.75)
  yq <- as.numeric(stats::quantile(y_sorted, probs, names = FALSE, na.rm = TRUE))
  xq <- stats::qnorm(probs)
  slope <- diff(yq) / diff(xq)
  intercept <- yq[1L] - slope * xq[1L]
  list(intercept = intercept, slope = slope)
}

generate_limma_residual_qq <- function(matrix_path, meta_path,
                                       method, task, domain,
                                       marker_genes) {
  if (!file.exists(matrix_path) || !file.exists(meta_path)) return(NULL)
  mat <- load_mat(matrix_path)
  meta <- fread(meta_path)
  setnames(meta, tolower(names(meta)))
  if (!all(c("sample_id", "domain", "condition") %in% names(meta))) return(NULL)

  dom_u <- toupper(domain)
  domain_meta <- meta[toupper(as.character(meta$domain)) == dom_u]
  domain_cols <- intersect(colnames(mat), domain_meta$sample_id)
  if (length(domain_cols) < 4L) return(NULL)

  domain_meta <- domain_meta[sample_id %in% domain_cols]
  domain_mat <- mat[, domain_cols, drop = FALSE]

  plots <- list()
  for (gene in marker_genes) {
    if (!gene %in% rownames(domain_mat)) next
    vals <- as.numeric(domain_mat[gene, , drop = TRUE])
    names(vals) <- colnames(domain_mat)
    df <- data.frame(
      sample_id = names(vals),
      abundance = vals,
      condition = domain_meta$condition[match(names(vals), domain_meta$sample_id)],
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$abundance) & !is.na(df$condition), , drop = FALSE]
    if (nrow(df) < 4L || length(unique(df$condition)) < 2L) next

    fit <- stats::lm(abundance ~ condition, data = df)
    df$residual <- stats::residuals(fit)
    n <- nrow(df)
    df <- df[order(df$residual), , drop = FALSE]
    df$theoretical <- stats::qnorm(stats::ppoints(n))
    ql <- qq_line_ab(df$residual)

    p <- ggplot(df, aes(x = theoretical, y = residual, color = condition)) +
      geom_point(size = 1.5, alpha = 0.7) +
      geom_abline(intercept = ql$intercept, slope = ql$slope, color = "red", linewidth = 0.4) +
      labs(title = gene, x = "", y = "") +
      theme_minimal(base_size = 9) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 9, face = "bold")
      )
    plots[[gene]] <- p
  }
  plots
}

generate_resid_vs_fitted <- function(matrix_path, meta_path,
                                       method, task, domain,
                                       marker_genes) {
  if (!file.exists(matrix_path) || !file.exists(meta_path)) return(NULL)
  mat <- load_mat(matrix_path)
  meta <- fread(meta_path)
  setnames(meta, tolower(names(meta)))
  if (!all(c("sample_id", "domain", "condition") %in% names(meta))) return(NULL)

  dom_u <- toupper(domain)
  domain_meta <- meta[toupper(as.character(meta$domain)) == dom_u]
  domain_cols <- intersect(colnames(mat), domain_meta$sample_id)
  if (length(domain_cols) < 4L) return(NULL)

  domain_meta <- domain_meta[sample_id %in% domain_cols]
  domain_mat <- mat[, domain_cols, drop = FALSE]

  plots <- list()
  for (gene in marker_genes) {
    if (!gene %in% rownames(domain_mat)) next
    vals <- as.numeric(domain_mat[gene, , drop = TRUE])
    df <- data.frame(
      abundance = vals,
      condition = domain_meta$condition[match(colnames(domain_mat), domain_meta$sample_id)],
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$abundance) & !is.na(df$condition), , drop = FALSE]
    if (nrow(df) < 4L) next

    fit <- stats::lm(abundance ~ condition, data = df)
    df$fitted <- stats::fitted(fit)
    df$residual <- stats::residuals(fit)

    p <- ggplot(df, aes(x = fitted, y = residual)) +
      geom_point(size = 1.2, alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.4, formula = y ~ x) +
      labs(title = gene, x = "", y = "") +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold"))
    plots[[gene]] <- p
  }
  plots
}

da_limma_path <- function(method, task, domain) {
  p1 <- file.path(
    REPO, "reports", "benchmark_master", "benchmark_results",
    method, task, "representation_da", domain, "da_limma_result.csv"
  )
  if (file.exists(p1)) return(p1)
  p2 <- file.path(
    REPO, "reports", "benchmark_master", "benchmark_results",
    method, task, "representation_da", paste0("da_", domain, ".csv")
  )
  if (file.exists(p2)) return(p2)
  character(0)
}

# ═══════════════════════════════════════════════════════════
# 1. QQ plots of lm(abundance ~ condition) residuals per marker
# ═══════════════════════════════════════════════════════════

for (task in tasks) {
  markers <- if (task == "breast_subtype") markers_subtype else markers_bvl
  meta_path <- file.path(REPO, "data", "processed", "union", paste0("sample_meta_", task, ".csv"))
  if (!file.exists(meta_path)) {
    meta_path <- file.path(REPO, "data", "processed", paste0("sample_meta_", task, ".csv"))
  }

  for (domain in domains) {
    all_plots <- list()
    row_labels <- character(0)

    for (method in methods) {
      mat_path <- file.path(REPO, "data", "processed", "methods", method, paste0("transformed_", task, ".csv"))
      gene_plots <- generate_limma_residual_qq(
        mat_path, meta_path, method, task, domain, markers
      )
      if (is.null(gene_plots) || length(gene_plots) == 0L) next

      for (g in markers) {
        if (!g %in% names(gene_plots)) gene_plots[[g]] <- empty_panel()
      }
      all_plots <- c(all_plots, gene_plots[markers])
      row_labels <- c(row_labels, method)
    }

    if (length(all_plots) > 0L) {
      n_methods <- length(row_labels)
      n_genes <- length(markers)
      grid <- plot_grid(plotlist = all_plots, nrow = n_methods, ncol = n_genes, align = "hv")

      row_label_plots <- lapply(row_labels, function(lab) {
        ggdraw() + draw_label(lab, size = 10, fontface = "bold", angle = 90)
      })
      label_col <- plot_grid(plotlist = row_label_plots, ncol = 1, align = "v")

      col_label_plots <- lapply(markers, function(g) {
        ggdraw() + draw_label(g, size = 10, fontface = "bold")
      })
      label_row <- plot_grid(plotlist = col_label_plots, nrow = 1)

      final <- plot_grid(
        plot_grid(ggdraw(), label_row, ncol = 2, rel_widths = c(0.08, 0.92)),
        plot_grid(label_col, grid, ncol = 2, rel_widths = c(0.08, 0.92)),
        ncol = 1, rel_heights = c(0.06, 0.94)
      )

      outpath <- file.path(OUT_FIG, sprintf("qq_residuals_%s_%s.pdf", domain, task))
      ggsave(outpath, final, width = 14, height = 3 * n_methods + 1, limitsize = FALSE)
      message("Wrote ", outpath)
    }
  }
}

# ═══════════════════════════════════════════════════════════
# 2. Residual vs fitted
# ═══════════════════════════════════════════════════════════

for (task in tasks) {
  markers <- if (task == "breast_subtype") markers_subtype else markers_bvl
  meta_path <- file.path(REPO, "data", "processed", "union", paste0("sample_meta_", task, ".csv"))
  if (!file.exists(meta_path)) {
    meta_path <- file.path(REPO, "data", "processed", paste0("sample_meta_", task, ".csv"))
  }

  for (domain in domains) {
    all_plots <- list()
    row_labels <- character(0)

    for (method in methods) {
      mat_path <- file.path(REPO, "data", "processed", "methods", method, paste0("transformed_", task, ".csv"))
      gene_plots <- generate_resid_vs_fitted(
        mat_path, meta_path, method, task, domain, markers
      )
      if (is.null(gene_plots) || length(gene_plots) == 0L) next

      for (g in markers) {
        if (!g %in% names(gene_plots)) gene_plots[[g]] <- empty_panel()
      }
      all_plots <- c(all_plots, gene_plots[markers])
      row_labels <- c(row_labels, method)
    }

    if (length(all_plots) > 0L) {
      n_methods <- length(row_labels)
      grid <- plot_grid(plotlist = all_plots, nrow = n_methods, ncol = length(markers), align = "hv")

      row_label_plots <- lapply(row_labels, function(lab) {
        ggdraw() + draw_label(lab, size = 10, fontface = "bold", angle = 90)
      })
      label_col <- plot_grid(plotlist = row_label_plots, ncol = 1, align = "v")
      col_label_plots <- lapply(markers, function(g) {
        ggdraw() + draw_label(g, size = 10, fontface = "bold")
      })
      label_row <- plot_grid(plotlist = col_label_plots, nrow = 1)
      final <- plot_grid(
        plot_grid(ggdraw(), label_row, ncol = 2, rel_widths = c(0.08, 0.92)),
        plot_grid(label_col, grid, ncol = 2, rel_widths = c(0.08, 0.92)),
        ncol = 1, rel_heights = c(0.06, 0.94)
      )

      outpath <- file.path(OUT_FIG, sprintf("resid_fitted_%s_%s.pdf", domain, task))
      ggsave(outpath, final, width = 14, height = 3 * n_methods + 1, limitsize = FALSE)
      message("Wrote ", outpath)
    }
  }
}

# ═══════════════════════════════════════════════════════════
# 3. SA-style plot: |logFC| vs AveExpr from limma topTable
# ═══════════════════════════════════════════════════════════

for (task in tasks) {
  for (domain in domains) {
    plots <- list()

    for (method in methods) {
      da_path <- da_limma_path(method, task, domain)
      if (length(da_path) == 0L) next
      da <- fread(da_path)
      if (!"AveExpr" %in% names(da)) next
      # Loess on >~8k points is very slow; subsample for plotting only
      if (nrow(da) > 8000L) {
        da <- da[sample.int(nrow(da), 8000L)]
      }

      da[, abs_logFC := abs(logFC)]
      yhi <- stats::quantile(da$abs_logFC, 0.99, na.rm = TRUE)
      if (!is.finite(yhi) || yhi <= 0) yhi <- max(da$abs_logFC, na.rm = TRUE)

      p <- ggplot(da, aes(x = AveExpr, y = abs_logFC)) +
        geom_point(alpha = 0.05, size = 0.3, color = "grey30") +
        geom_smooth(method = "loess", color = "red", linewidth = 0.6, se = FALSE, formula = y ~ x) +
        labs(title = method, x = "Mean expression (AveExpr)", y = "|logFC|") +
        coord_cartesian(ylim = c(0, yhi)) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold"))

      plots[[method]] <- p
    }

    if (length(plots) > 0L) {
      grid <- plot_grid(plotlist = plots, nrow = 1, align = "h")
      title_plot <- ggdraw() +
        draw_label(sprintf("SA-style plot: %s — %s", toupper(domain), task),
                   size = 12, fontface = "bold")
      final <- plot_grid(title_plot, grid, ncol = 1, rel_heights = c(0.07, 0.93))

      outpath <- file.path(OUT_FIG, sprintf("sa_%s_%s.pdf", domain, task))
      ggsave(outpath, final, width = 16, height = 4, limitsize = FALSE)
      message("Wrote ", outpath)
    }
  }
}

# ═══════════════════════════════════════════════════════════
# 4. Quantitative summary (marginal logFC distribution + AveExpr correlation)
# ═══════════════════════════════════════════════════════════

compute_assumption_stats <- function(da_path, method, task, domain) {
  da <- fread(da_path)
  lfc <- da$logFC[is.finite(da$logFC)]
  if (length(lfc) < 10L) return(NULL)

  lfc_sample <- if (length(lfc) > 5000L) sample(lfc, 5000L) else lfc
  sw <- tryCatch(
    stats::shapiro.test(lfc_sample),
    error = function(e) list(p.value = NA_real_)
  )

  skew <- mean((lfc - mean(lfc))^3) / (stats::sd(lfc)^3)
  kurt <- mean((lfc - mean(lfc))^4) / (stats::sd(lfc)^4) - 3

  mv_cor <- NA_real_
  if ("AveExpr" %in% names(da)) {
    mv_cor <- suppressWarnings(stats::cor(abs(da$logFC), da$AveExpr, use = "complete.obs"))
  }

  extreme_frac <- NA_real_
  if ("t" %in% names(da)) {
    extreme_frac <- mean(abs(da$t) > 5, na.rm = TRUE)
  }

  data.frame(
    method = method,
    task = task,
    domain = toupper(domain),
    n_genes = length(lfc),
    shapiro_p = signif(sw$p.value, 3),
    skewness = round(skew, 3),
    excess_kurtosis = round(kurt, 3),
    mean_var_cor = round(mv_cor, 3),
    frac_extreme_t = round(extreme_frac, 4),
    stringsAsFactors = FALSE
  )
}

all_stats <- data.frame()
for (method in methods) {
  for (task in tasks) {
    for (domain in domains) {
      da_path <- da_limma_path(method, task, domain)
      if (length(da_path) == 0L) next
      s <- compute_assumption_stats(da_path, method, task, domain)
      if (!is.null(s)) all_stats <- rbind(all_stats, s)
    }
  }
}

message("\n=== ASSUMPTION SUMMARY ===")
print(all_stats)

write.csv(all_stats, file.path(OUT_TAB, "assumption_summary.csv"), row.names = FALSE)

message("\n=== FLAGS ===")
for (i in seq_len(nrow(all_stats))) {
  r <- all_stats[i, ]
  flags <- character(0)
  if (!is.na(r$skewness) && abs(r$skewness) > 1.0) {
    flags <- c(flags, sprintf("skew=%.2f", r$skewness))
  }
  if (!is.na(r$excess_kurtosis) && abs(r$excess_kurtosis) > 3.0) {
    flags <- c(flags, sprintf("kurtosis=%.2f", r$excess_kurtosis))
  }
  if (!is.na(r$mean_var_cor) && abs(r$mean_var_cor) > 0.2) {
    flags <- c(flags, sprintf("mean-var r=%.2f", r$mean_var_cor))
  }
  if (!is.na(r$frac_extreme_t) && r$frac_extreme_t > 0.02) {
    flags <- c(flags, sprintf("%.1f%% extreme t", r$frac_extreme_t * 100))
  }

  if (length(flags) > 0L) {
    message(sprintf(
      "  WARNING %s / %s / %s: %s",
      r$method, r$task, r$domain, paste(flags, collapse = "; ")
    ))
  }
}

message("\nDone. Outputs: ", OUT_FIG, " and ", file.path(OUT_TAB, "assumption_summary.csv"))
