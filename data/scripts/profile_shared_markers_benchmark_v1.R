#!/usr/bin/env Rscript
# Profile plots (log2 abundance by subtype) for genes expected same-direction in both cohorts.
# Output: reports/benchmark_v1/diagnostics/shared_protein_profiles/*.pdf
#
# Uses gene_matrix (gene-level), not raw protein — labeled explicitly in plot subtitles.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

find_root <- function() {
  wd <- getwd()
  if (file.exists(file.path(wd, "data", "results", "PDC000120", "gene_matrix.csv")))
    return(normalizePath(file.path(wd, "data")))
  stop("Run from repo root.")
}
root <- find_root()
out_base <- file.path(root, "..", "reports", "benchmark_v1", "diagnostics_feedback", "shared_protein_profiles")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

des <- fread(file.path(root, "results", "PDC000120", "DA_subtype_tumor_only_basal_luminal_subset.csv"))
des[, pam50 := trimws(as.character(pam50))]
des[pam50 %in% c("LumA", "LumB"), pam50 := "Luminal"]
des <- des[pam50 %in% c("Luminal", "Basal")]
id_col <- "matrix_sample_id"
gm <- fread(file.path(root, "results", "PDC000120", "gene_matrix.csv"))
gc <- names(gm)[1L]
cols <- setdiff(names(gm), gc)
des_l <- tolower(trimws(des[[id_col]]))
col_l <- tolower(trimws(cols))
mc <- cols[match(des_l, col_l)]
des <- des[!is.na(mc)]

gm_c <- fread(file.path(root, "results", "CCLE_corrected", "gene_matrix.csv"))
lum <- c("MCF7", "T-47D", "CAMA-1", "ZR-75-1")
bas <- c("HCC 1806", "HCC1143", "HCC70", "MDA-MB-468")

genes <- c("FOXA1", "EGFR", "KRT5")
for (gene in genes) {
  r <- which(trimws(as.character(gm[[gc]])) == gene)
  if (!length(r)) next
  y <- as.numeric(gm[r[1L], ..mc][1L, ])
  dt <- data.table(abundance = y, Subtype = des$pam50, Cohort = "CPTAC")
  rc <- which(trimws(as.character(gm_c[[gc]])) == gene)
  if (length(rc)) {
    sc <- intersect(names(gm_c), c(lum, bas))
    y2 <- as.numeric(gm_c[rc[1L], ..sc][1L, ])
    dt2 <- data.table(
      abundance = y2,
      Subtype = rep(c("Luminal", "Basal"), c(4L, 4L)),
      Cohort = "CCLE"
    )
    dt <- rbind(dt, dt2)
  }
  p <- ggplot(dt, aes(Subtype, abundance, fill = Subtype)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(width = 0.12, size = 1, alpha = 0.6) +
    facet_wrap(~Cohort, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste0(gene, " — log2 abundance by subtype (gene matrix)"),
      subtitle = "CPTAC: mixture-balanced tumors; CCLE: 4+4 lines; not raw reporter-ion scale",
      y = "log2 abundance",
      x = NULL
    )
  ggsave(file.path(out_base, paste0("profile_", gene, ".pdf")), p, width = 8, height = 4, dpi = 200)
}

message("Wrote profiles to ", out_base)
