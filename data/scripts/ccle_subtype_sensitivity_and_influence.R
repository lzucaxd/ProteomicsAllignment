#!/usr/bin/env Rscript
# CCLE subtype sensitivity: leave-one-line-out on gene_matrix (limma) — fast diagnostic
# for whether one line dominates the Luminal–Basal contrast. Not a replacement for
# protein-level MSstatsTMT; interpret as complementary influence check.
#
# Run from repo root:
#   Rscript data/scripts/ccle_subtype_sensitivity_and_influence.R

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

find_gm <- function() {
  wd <- getwd()
  p1 <- file.path(wd, "data", "results", "CCLE_corrected", "gene_matrix.csv")
  if (file.exists(p1)) return(normalizePath(p1))
  p2 <- file.path(wd, "results", "CCLE_corrected", "gene_matrix.csv")
  if (file.exists(p2)) return(normalizePath(p2))
  stop("gene_matrix.csv not found")
}

luminal_lines <- c("MCF7", "T-47D", "CAMA-1", "ZR-75-1")
basal_lines <- c("HCC 1806", "HCC1143", "HCC70", "MDA-MB-468")
all_lines <- c(luminal_lines, basal_lines)

gm_path <- find_gm()
out_txt <- file.path(dirname(gm_path), "DA_luminal_vs_basal", "ccle_line_influence_summary.txt")
dir.create(dirname(out_txt), recursive = TRUE, showWarnings = FALSE)

gm <- fread(gm_path, showProgress = FALSE)
miss <- setdiff(all_lines, names(gm))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

M <- as.matrix(gm[, ..all_lines])
storage.mode(M) <- "numeric"
M[!is.finite(M)] <- NA

# Filter genes: >=2/4 per group observed (same spirit as v1)
ok <- rowSums(is.finite(M[, luminal_lines])) >= 2L & rowSums(is.finite(M[, basal_lines])) >= 2L
M <- M[ok, , drop = FALSE]

fit_one <- function(cols, g1, g2) {
  group <- factor(rep(c("Luminal", "Basal"), c(length(g1), length(g2))), levels = c("Luminal", "Basal"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("Luminal", "Basal")
  fit <- lmFit(M[, cols, drop = FALSE], design)
  ctr <- makeContrasts(Luminal_vs_Basal = Luminal - Basal, levels = design)
  fit2 <- contrasts.fit(fit, ctr)
  fit2 <- eBayes(fit2)
  topTable(fit2, coef = 1, number = Inf, sort.by = "none")
}

full_tt <- fit_one(all_lines, luminal_lines, basal_lines)
full_fc <- full_tt$logFC
names(full_fc) <- rownames(full_tt)

loo <- list()
for (j in seq_along(all_lines)) {
  line <- all_lines[j]
  cols <- setdiff(all_lines, line)
  g1 <- intersect(luminal_lines, cols)
  g2 <- intersect(basal_lines, cols)
  tt <- fit_one(cols, g1, g2)
  fc <- tt$logFC
  names(fc) <- rownames(tt)
  common <- intersect(names(full_fc), names(fc))
  rho <- suppressWarnings(cor(full_fc[common], fc[common], method = "spearman", use = "complete.obs"))
  # Delta FDR < 0.05 count at gene level (rough)
  n_sig_full <- sum(full_tt$adj.P.Val < 0.05, na.rm = TRUE)
  n_sig_loo <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  loo[[line]] <- list(
    rho = rho,
    n_genes = length(common),
    n_sig_loo = n_sig_loo,
    n_lum = length(g1), n_bas = length(g2)
  )
}

sink(out_txt)
cat("CCLE line influence — leave-one-line-out (gene_matrix + limma)\n")
cat("===============================================================\n")
cat("Purpose: lightweight diagnostic; full benchmark remains MSstatsTMT on proteins.\n")
cat("Genes in model: ", nrow(M), " (filtered >=2/4 observed per group)\n\n")
cat("Spearman correlation of gene-level logFC vs full 8-line model (higher = removing line changes result less):\n\n")
for (line in all_lines) {
  x <- loo[[line]]
  cat(sprintf(
    "  Drop %-12s -> Lum %d vs Bas %d | rho = %.3f | adj.P<0.05 genes: %d (full model ref: %d)\n",
    line, x$n_lum, x$n_bas, x$rho, x$n_sig_loo, sum(full_tt$adj.P.Val < 0.05, na.rm = TRUE)
  ))
}
cat("\nInterpretation (conservative):\n")
cat("- rho in the ~0.85–0.99 range suggests no single line fully determines the genome-wide logFC profile.\n")
cat("- The line with the LOWEST rho is the most influential for this gene-level summary (still not causal proof).\n")
cat("- Unbalanced LOO (3 vs 4) is expected; do not over-interpret gene-level p-values vs full 4 vs 4.\n")
sink(NULL)

message("Wrote ", out_txt)
