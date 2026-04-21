#!/usr/bin/env Rscript
# QQ plots with points colored by TMT Mixture (plex) — Olga feedback item 1.
# Gene-level log2 abundances vs normal quantiles; one panel per marker gene.
#
# Usage (repo root):
#   Rscript --vanilla data/scripts/diagnostics_qq_by_mixture_benchmark_v1.R [out_dir]
#
# Requires: data/results/PDC000120/gene_matrix.csv, DA_subtype_tumor_only_basal_luminal_subset.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

find_data_root <- function() {
  wd <- getwd()
  if (file.exists(file.path(wd, "data", "results", "PDC000120", "gene_matrix.csv")))
    return(normalizePath(file.path(wd, "data")))
  if (file.exists(file.path(wd, "results", "PDC000120", "gene_matrix.csv")))
    return(normalizePath(wd))
  stop("Run from repo root.")
}
root <- find_data_root()
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[1L] else file.path(root, "..", "reports", "benchmark_v1", "diagnostics")
out_dir <- normalizePath(out_dir, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

des_path <- file.path(root, "results", "PDC000120", "DA_subtype_tumor_only_basal_luminal_subset.csv")
gm_path <- file.path(root, "results", "PDC000120", "gene_matrix.csv")
if (!file.exists(des_path) || !file.exists(gm_path))
  stop("Missing subset annotation or gene_matrix: ", des_path)

des <- fread(des_path)
des[, pam50 := trimws(as.character(pam50))]
des[pam50 %in% c("LumA", "LumB"), pam50 := "Luminal"]
des <- des[pam50 %in% c("Luminal", "Basal")]
id_col <- if ("matrix_sample_id" %in% names(des)) "matrix_sample_id" else names(des)[1L]
mix_col <- if ("mixture" %in% names(des)) "mixture" else NA_character_
if (is.na(mix_col)) stop("No mixture column in subset annotation")

gm <- fread(gm_path)
gene_col <- names(gm)[1L]
cols <- setdiff(names(gm), intersect(c(gene_col, "UniProtID", "uniprotid"), names(gm)))
des_lower <- tolower(trimws(des[[id_col]]))
col_lower <- tolower(trimws(cols))
mat_cols <- cols[match(des_lower, col_lower)]
ok <- !is.na(mat_cols)
des <- des[ok]
mat_cols <- mat_cols[ok]
mix <- des[[mix_col]]

MARKERS <- c("FOXA1", "GATA3", "EGFR", "KRT5")

plots <- list()
for (g in MARKERS) {
  r <- which(trimws(as.character(gm[[gene_col]])) == g)
  if (!length(r)) {
    warning("Gene not in matrix: ", g)
    next
  }
  y <- as.numeric(gm[r[1L], ..mat_cols][1L, ])
  n <- qqnorm(y, plot.it = FALSE)
  dt <- data.table(
    theoretical = n$x,
    sample = n$y,
    Mixture = factor(mix),
    gene = g
  )
  dt <- dt[is.finite(theoretical) & is.finite(sample)]
  p <- ggplot(dt, aes(theoretical, sample, color = Mixture)) +
    geom_point(alpha = 0.85, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "gray40", linewidth = 0.3) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(
      title = paste0("CPTAC QQ — ", g, " (log2 gene abundance)"),
      subtitle = paste0(
        "Mixture-balanced subtype subset; n=", nrow(des), " tumors; points colored by plex"
      ),
      x = "Theoretical normal quantiles",
      y = "Sample quantiles (log2)"
    )
  plots[[g]] <- p
}

if (!length(plots)) stop("No plots generated.")

suppressPackageStartupMessages(library(gridExtra))
plist <- unname(plots)
combined <- grid.arrange(grobs = plist, ncol = 2L)

ggsave(file.path(out_dir, "qq_by_plex.pdf"), combined, width = 12, height = 10, dpi = 200)
ggsave(file.path(out_dir, "qq_by_plex.png"), combined, width = 12, height = 10, dpi = 200)
message("Wrote ", file.path(out_dir, "qq_by_plex.pdf"))
