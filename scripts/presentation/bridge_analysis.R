#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
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

explain_bridge <- function(matrix_path, meta_path, method, task, repo, out_root) {
  mat <- load_gene_matrix(matrix_path)
  meta <- fread(meta_path)
  meta <- meta[sample_id %in% colnames(mat)]
  cptac_ids <- meta[domain == "CPTAC", sample_id]
  ccle_ids <- meta[domain == "CCLE", sample_id]

  rep_genes <- c("ESR1", "GATA3", "EGFR", "KRT5", "ACTB", "GAPDH")
  rep_genes <- rep_genes[rep_genes %in% rownames(mat)]

  results <- data.table()
  for (g in rep_genes) {
    vals_c <- as.numeric(mat[g, colnames(mat) %in% cptac_ids, drop = TRUE])
    vals_e <- as.numeric(mat[g, colnames(mat) %in% ccle_ids, drop = TRUE])
    results <- rbind(results, data.table(
      gene = g, domain = "CPTAC",
      mean = round(mean(vals_c, na.rm = TRUE), 3L),
      sd = round(stats::sd(vals_c, na.rm = TRUE), 3L),
      n_measured = sum(!is.na(vals_c))
    ))
    results <- rbind(results, data.table(
      gene = g, domain = "CCLE",
      mean = round(mean(vals_e, na.rm = TRUE), 3L),
      sd = round(stats::sd(vals_e, na.rm = TRUE), 3L),
      n_measured = sum(!is.na(vals_e))
    ))
  }

  wide <- pivot_wider(
    results,
    names_from = domain,
    values_from = c(mean, sd, n_measured)
  )
  wide[, location_shift := mean_CPTAC - mean_CCLE]
  wide[, scale_ratio := fifelse(sd_CCLE > 0, sd_CPTAC / sd_CCLE, NA_real_)]

  cat(sprintf("\n=== BRIDGE ANALYSIS: %s / %s ===\n", method, task))
  print(as.data.frame(wide))

  fwrite(wide, file.path(out_root, sprintf("tables/bridge_analysis_%s_%s.csv", method, task)))
  invisible(wide)
}

explain_bridge(
  file.path(REPO, "data/processed/methods/raw/transformed_breast_subtype.csv"),
  file.path(REPO, "data/processed/union/sample_meta_breast_subtype.csv"),
  "raw", "breast_subtype", REPO, PRES_OUT
)
explain_bridge(
  file.path(REPO, "data/processed/methods/bridge_shift/transformed_breast_subtype.csv"),
  file.path(REPO, "data/processed/union/sample_meta_breast_subtype.csv"),
  "bridge_shift", "breast_subtype", REPO, PRES_OUT
)
