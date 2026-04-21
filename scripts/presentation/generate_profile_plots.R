#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1]))) else normalizePath(file.path(getwd(), "scripts", "presentation"))
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Install cowplot: install.packages('cowplot') — added to install_r_packages.R for this repo.")
}

load_gene_matrix <- function(path) {
  dt <- fread(path)
  gn <- names(dt)[1L]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[gn]]
  colnames(mat) <- names(dt)[-1L]
  mat
}

generate_profile_plot <- function(gene, mat, meta, task, method) {
  if (!gene %in% rownames(mat)) {
    message(sprintf("WARNING: %s not in matrix (%s / %s)", gene, method, task))
    return(NULL)
  }
  vals <- as.numeric(mat[gene, , drop = TRUE])
  names(vals) <- colnames(mat)
  meta <- copy(meta)
  meta <- meta[sample_id %in% names(vals)]
  if (nrow(meta) == 0L) return(NULL)
  meta[, abundance := vals[sample_id]]

  meta[, block := paste0(domain, "\n", condition)]
  block_order <- c("CPTAC\nLuminal", "CCLE\nLuminal", "CPTAC\nBasal", "CCLE\nBasal")
  meta[, block := factor(block, levels = block_order)]
  meta <- meta[!is.na(block)]

  medians <- meta[, .(abundance = median(abundance, na.rm = TRUE)), by = block]

  ggplot(meta, aes(x = block, y = abundance)) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 1.5, color = "grey50") +
    geom_crossbar(
      data = medians,
      aes(x = block, y = abundance, ymin = abundance, ymax = abundance),
      width = 0.5,
      color = "black",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey70") +
    labs(
      title = sprintf("%s — %s (%s)", gene, task, method),
      x = NULL, y = "Abundance (log2 scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
}

tasks <- "breast_subtype"
methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
markers <- c("ESR1", "GATA3", "EGFR", "KRT5", "FOXA1", "FOXC1")

for (task in tasks) {
  meta_path <- file.path(REPO, sprintf("data/processed/union/sample_meta_%s.csv", task))
  if (!file.exists(meta_path)) {
    message("SKIP task ", task, ": missing ", meta_path)
    next
  }
  meta <- fread(meta_path)

  for (method in methods) {
    mat_path <- file.path(REPO, sprintf("data/processed/methods/%s/transformed_%s.csv", method, task))
    if (!file.exists(mat_path)) {
      message("SKIP ", method, "/", task)
      next
    }
    mat <- load_gene_matrix(mat_path)
    for (gene in markers) {
      p <- generate_profile_plot(gene, mat, meta, task, method)
      if (!is.null(p)) {
        outp <- file.path(PRES_OUT, sprintf("figures/profile_%s_%s_%s.pdf", gene, method, task))
        ggsave(outp, p, width = 6, height = 4)
      }
    }
  }
}

for (gene in c("ESR1", "EGFR")) {
  plots <- list()
  mlabs <- character(0)
  meta_path <- file.path(REPO, "data/processed/union/sample_meta_breast_subtype.csv")
  if (!file.exists(meta_path)) next
  meta <- fread(meta_path)
  for (method in methods) {
    mat_path <- file.path(REPO, sprintf("data/processed/methods/%s/transformed_breast_subtype.csv", method))
    if (!file.exists(mat_path)) next
    mat <- load_gene_matrix(mat_path)
    p <- generate_profile_plot(gene, mat, meta, "breast_subtype", method)
    if (!is.null(p)) {
      plots[[length(plots) + 1L]] <- p
      mlabs <- c(mlabs, method)
    }
  }
  if (length(plots) > 0L) {
    combined <- cowplot::plot_grid(plotlist = plots, nrow = 1L, labels = mlabs)
    ggsave(
      file.path(PRES_OUT, sprintf("figures/profile_combined_%s_breast_subtype.pdf", gene)),
      combined,
      width = 5 * length(plots),
      height = 4
    )
  }
}

message("Profile plots done under ", file.path(PRES_OUT, "figures"))
