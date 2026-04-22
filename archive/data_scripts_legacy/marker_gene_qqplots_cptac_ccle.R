#!/usr/bin/env Rscript
# Per-gene QQ plots for subtype marker genes (CPTAC subset + CCLE 8 lines).
# QQ compares the empirical distribution of log2 abundances across samples to a
# normal distribution with matching mean and SD (standard qqnorm).
#
# Usage (from repo root):
#   Rscript data/scripts/marker_gene_qqplots_cptac_ccle.R [out_dir]
# Default out_dir: data/results/PDC000120/diagnostics/marker_qq_by_gene
#
# Notes for slides:
# - CPTAC: one point per tumor sample in the mixture-balanced subset (same samples as
#   MSstatsTMT subtype subset run).
# - CCLE: one point per cell line (n=8 total) — QQ is very sparse; interpret as a
#   requested visual, not a formal normality test.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Resolve .../data (repo root vs data/ as cwd)
find_data_root <- function() {
  wd <- getwd()
  if (file.exists(file.path(wd, "data", "results", "PDC000120", "gene_matrix.csv")))
    return(normalizePath(file.path(wd, "data")))
  if (file.exists(file.path(wd, "results", "PDC000120", "gene_matrix.csv")))
    return(normalizePath(wd))
  stop("Cannot find data/results/PDC000120/gene_matrix.csv (run from repo root or data/).")
}
root <- find_data_root()

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[1L] else file.path(root, "results", "PDC000120", "diagnostics", "marker_qq_by_gene")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Two Luminal-associated + two Basal-associated (gene symbols; must exist in each matrix)
LUM_GENES <- c("FOXA1", "GATA3")
BAS_GENES <- c("KRT5", "EGFR")
MARKERS <- c(LUM_GENES, BAS_GENES)

# --- CPTAC: subset design ---
subset_csv <- file.path(root, "results", "PDC000120", "DA_subtype_tumor_only_basal_luminal_subset.csv")
if (!file.exists(subset_csv))
  stop("Missing subset annotation: ", subset_csv)
des <- fread(subset_csv)
des[, pam50 := trimws(as.character(pam50))]
des[pam50 %in% c("LumA", "LumB"), pam50 := "Luminal"]
des <- des[pam50 %in% c("Luminal", "Basal")]
id_col <- if ("matrix_sample_id" %in% names(des)) "matrix_sample_id" else names(des)[1L]

gm <- fread(file.path(root, "results", "PDC000120", "gene_matrix.csv"))
gene_col <- names(gm)[1L]
cols <- setdiff(names(gm), intersect(c(gene_col, "UniProtID", "uniprotid"), names(gm)))
des_lower <- tolower(trimws(des[[id_col]]))
col_lower <- tolower(trimws(cols))
mat_cols <- cols[match(des_lower, col_lower)]
ok <- !is.na(mat_cols)
des <- des[ok]
mat_cols <- mat_cols[ok]

extract_row <- function(g, gm_dt, gene_col) {
  r <- which(trimws(as.character(gm_dt[[gene_col]])) == g)
  if (!length(r)) return(NULL)
  if (length(r) > 1L) r <- r[1L] # first row if duplicates
  as.numeric(gm_dt[r, ..mat_cols][1, ])
}

plot_qq_cohort <- function(vecs, cohort_title, subtitle) {
  rows <- lapply(seq_along(vecs), function(i) {
    y <- vecs[[i]]
    y <- y[is.finite(y)]
    data.table(gene = names(vecs)[i], value = y)
  })
  dt <- rbindlist(rows)
  dt[, gene := factor(gene, levels = MARKERS)]
  ggplot(dt, aes(sample = value)) +
    stat_qq(alpha = 0.85, size = 1.2) +
    stat_qq_line(color = "firebrick", linewidth = 0.4) +
    facet_wrap(~gene, ncol = 2, scales = "free") +
    theme_bw(base_size = 11) +
    labs(
      title = cohort_title,
      subtitle = subtitle,
      x = "Theoretical quantiles (normal)",
      y = "Sample quantiles (log2 abundance)"
    )
}

vecs_cptac <- setNames(vector("list", length(MARKERS)), MARKERS)
for (g in MARKERS) {
  v <- extract_row(g, gm, gene_col)
  if (is.null(v)) {
    warning("CPTAC: gene not found in gene_matrix: ", g)
    vecs_cptac[[g]] <- numeric(0)
  } else {
    vecs_cptac[[g]] <- v
  }
}

n_lum <- sum(des$pam50 == "Luminal")
n_bas <- sum(des$pam50 == "Basal")
p_cptac <- plot_qq_cohort(
  vecs_cptac,
  "CPTAC PDC000120 — QQ of log2 gene abundance (mixture-balanced subset)",
  paste0(
    "n = ", length(mat_cols), " tumors (", n_lum, " Luminal [LumA+LumB] + ", n_bas, " Basal); ",
    "same sample set as MSstatsTMT subtype subset; matrix: gene_matrix.csv"
  )
)

# --- CCLE corrected ---
gm_c <- fread(file.path(root, "results", "CCLE_corrected", "gene_matrix.csv"))
luminal <- c("MCF7", "T-47D", "CAMA-1", "ZR-75-1")
basal <- c("HCC 1806", "HCC1143", "HCC70", "MDA-MB-468")
samples <- c(luminal, basal)
miss <- setdiff(samples, names(gm_c))
if (length(miss)) stop("CCLE: missing columns: ", paste(miss, collapse = ", "))

vecs_ccle <- setNames(vector("list", length(MARKERS)), MARKERS)
for (g in MARKERS) {
  r <- which(trimws(as.character(gm_c[[1]])) == g)
  if (!length(r)) {
    warning("CCLE: gene not found: ", g)
    vecs_ccle[[g]] <- numeric(0)
    next
  }
  if (length(r) > 1L) r <- r[1L]
  vecs_ccle[[g]] <- suppressWarnings(as.numeric(gm_c[r, ..samples][1, ]))
}

p_ccle <- plot_qq_cohort(
  vecs_ccle,
  "CCLE (corrected gene matrix) — QQ of log2 gene abundance",
  "n = 8 cell lines (4 Luminal + 4 Basal); QQ is very sparse — show for completeness, not as formal diagnostics"
)

out_pdf1 <- file.path(out_dir, "marker_genes_QQ_CPTAC_subset.pdf")
out_pdf2 <- file.path(out_dir, "marker_genes_QQ_CCLE_corrected.pdf")
out_pdf_comb <- file.path(out_dir, "marker_genes_QQ_CPTAC_and_CCLE_combined.pdf")
# Use base pdf() — avoids cairo/X11 dependency issues on some macOS setups
ggsave(out_pdf1, p_cptac, width = 9, height = 7, device = grDevices::pdf)
ggsave(out_pdf2, p_ccle, width = 9, height = 7, device = grDevices::pdf)

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- p_cptac / p_ccle + plot_layout(heights = c(1, 1))
  ggsave(out_pdf_comb, combined, width = 10, height = 14, device = grDevices::pdf)
  message("Wrote combined: ", normalizePath(out_pdf_comb, mustWork = FALSE))
} else {
  message("Note: install.packages('patchwork') to also write a single combined PDF.")
}

message("Wrote: ", normalizePath(out_pdf1, mustWork = FALSE))
message("Wrote: ", normalizePath(out_pdf2, mustWork = FALSE))
