#!/usr/bin/env Rscript
# Legacy / optional: limma on CCLE gene_matrix only (gene-level), no MSstatsTMT groupComparisonTMT.
# Primary CCLE subtype DA is ccle_DA_luminal_basal_v1.R (MSstatsTMT on protein_summary).
#
# Run from repo root: Rscript data/scripts/ccle_DA_luminal_basal_limma_gene_matrix.R

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Install ggrepel: install.packages(\"ggrepel\")")
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
root <- getwd()
find_gm <- function() {
  candidates <- c(
    file.path(root, "data", "results", "CCLE_corrected", "gene_matrix.csv"),
    file.path(root, "data", "results", "CCLE", "gene_matrix.csv")
  )
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  NULL
}
if (length(args) >= 1L) {
  gm_path <- args[1L]
} else {
  gm_path <- find_gm()
  if (is.null(gm_path)) stop("Cannot find gene_matrix.csv")
}
out_dir <- if (length(args) >= 2L) args[2L] else file.path(dirname(gm_path), "DA_luminal_vs_basal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

luminal <- c("MCF7", "T-47D", "CAMA-1", "ZR-75-1")
basal   <- c("HCC 1806", "HCC1143", "HCC70", "MDA-MB-468")
samples <- c(luminal, basal)

gm <- fread(gm_path, showProgress = FALSE)
miss <- setdiff(samples, names(gm))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

uid <- gm[[2]]
M <- as.matrix(gm[, ..samples])
storage.mode(M) <- "numeric"
rownames(M) <- uid
M[!is.finite(M)] <- NA

ok <- rowSums(is.finite(M[, luminal])) >= 2L & rowSums(is.finite(M[, basal])) >= 2L
M <- M[ok, , drop = FALSE]

group <- factor(rep(c("Luminal", "Basal"), c(4L, 4L)), levels = c("Luminal", "Basal"))
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Luminal", "Basal")
fit <- lmFit(M, design)
contr <- makeContrasts(Luminal_vs_Basal = Luminal - Basal, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
tt <- as.data.table(topTable(fit2, coef = 1, number = Inf, sort.by = "P"), keep.rownames = "UniProtID")
fwrite(tt, file.path(out_dir, "DA_luminal_vs_basal_limma.csv"))
message("Wrote ", file.path(out_dir, "DA_luminal_vs_basal_limma.csv"))
