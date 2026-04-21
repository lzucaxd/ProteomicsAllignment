#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1]))) else normalizePath(file.path(getwd(), "scripts", "presentation"))
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

load_gene_matrix <- function(path) {
  dt <- fread(path)
  gn <- names(dt)[1L]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[gn]]
  colnames(mat) <- names(dt)[-1L]
  mat
}

generate_qq_by_plex <- function(mat, meta, gene, task) {
  if (!gene %in% rownames(mat)) return(NULL)
  vals <- as.numeric(mat[gene, , drop = TRUE])
  names(vals) <- colnames(mat)
  meta <- copy(meta)
  meta <- meta[sample_id %in% names(vals)]
  meta[, abundance := vals[sample_id]]
  meta <- meta[!is.na(abundance)]

  plex_col <- intersect(c("mixture", "plex", "Mixture", "run", "study_id"), names(meta))
  if (length(plex_col) == 0L) {
    meta[, plex_label := domain]
  } else {
    meta[, plex_label := as.factor(get(plex_col[1L]))]
  }

  meta <- meta[order(abundance)]
  meta[, theoretical := stats::qnorm((seq_len(.N) - 0.375) / (.N + 0.25))]

  ggplot(meta, aes(x = theoretical, y = abundance, color = plex_label)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(
      intercept = mean(meta$abundance, na.rm = TRUE),
      slope = stats::sd(meta$abundance, na.rm = TRUE),
      color = "red",
      linewidth = 0.5
    ) +
    labs(
      title = sprintf("QQ: %s (%s)", gene, task),
      x = "Theoretical quantiles (normal)",
      y = "Sample quantiles (abundance)",
      color = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

generate_resid_vs_fitted <- function(mat, meta, gene, task) {
  if (!gene %in% rownames(mat)) return(NULL)
  vals <- as.numeric(mat[gene, , drop = TRUE])
  names(vals) <- colnames(mat)
  meta <- copy(meta)
  meta <- meta[sample_id %in% names(vals)]
  meta[, abundance := vals[sample_id]]
  meta <- meta[!is.na(abundance)]
  meta[, condition := factor(condition)]

  fit <- stats::lm(abundance ~ condition, data = meta)
  meta[, fitted := stats::fitted(fit)]
  meta[, residual := stats::residuals(fit)]

  ggplot(meta, aes(x = fitted, y = residual)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.5, formula = y ~ x) +
    labs(
      title = sprintf("Residual vs Fitted: %s (%s)", gene, task),
      x = sprintf("Fitted (%s)", gene),
      y = "Residual"
    ) +
    theme_minimal(base_size = 12)
}

generate_sa_plot <- function(da_path, domain, task) {
  da <- fread(da_path)
  if (!all(c("AveExpr", "logFC") %in% names(da))) {
    message("SKIP SA plot: no AveExpr in ", da_path)
    return(NULL)
  }
  da[, abs_logFC := abs(logFC)]
  ggplot(da, aes(x = AveExpr, y = abs_logFC)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_smooth(method = "loess", color = "red", linewidth = 0.7, formula = y ~ x) +
    labs(
      title = sprintf("Effect size vs mean expression: %s %s", toupper(domain), task),
      x = "Average expression", y = "|logFC|"
    ) +
    theme_minimal(base_size = 12)
}

mat_path <- file.path(REPO, "data/processed/methods/raw/transformed_breast_subtype.csv")
meta_path <- file.path(REPO, "data/processed/union/sample_meta_breast_subtype.csv")
if (file.exists(mat_path) && file.exists(meta_path)) {
  mat <- load_gene_matrix(mat_path)
  meta <- fread(meta_path)
  for (gene in c("FOXA1", "GATA3", "EGFR", "KRT5")) {
    p <- generate_qq_by_plex(mat, meta, gene, "breast_subtype")
    if (!is.null(p)) {
      ggsave(file.path(PRES_OUT, sprintf("figures/qq_%s_subtype.pdf", gene)), p, width = 7, height = 5)
    }
  }
  for (gene in c("ESR1", "GATA3", "EGFR", "KRT5")) {
    p <- generate_resid_vs_fitted(mat, meta, gene, "breast_subtype")
    if (!is.null(p)) {
      ggsave(file.path(PRES_OUT, sprintf("figures/resid_fitted_%s.pdf", gene)), p, width = 6, height = 4)
    }
  }
}

for (domain in c("cptac", "ccle")) {
  da_path <- file.path(REPO, sprintf(
    "reports/benchmark_master/benchmark_results/raw/breast_subtype/representation_da/%s/da_limma_result.csv",
    domain
  ))
  if (file.exists(da_path)) {
    p <- generate_sa_plot(da_path, domain, "breast_subtype")
    if (!is.null(p)) {
      ggsave(file.path(PRES_OUT, sprintf("figures/sa_plot_%s_subtype.pdf", domain)), p, width = 6, height = 4)
    }
  }
}

message("Assumption-check figures written under ", file.path(PRES_OUT, "figures"))
