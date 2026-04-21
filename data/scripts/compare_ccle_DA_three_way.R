#!/usr/bin/env Rscript
# Compare luminal vs basal limma results: Table S2 vs old gene_matrix vs CCLE_corrected
suppressPackageStartupMessages(library(data.table))

strip_acc <- function(x) {
  x <- as.character(x)
  vapply(strsplit(sub("^sp\\|", "", x), "|", fixed = TRUE), function(z) z[[1]], character(1))
}

sum_dt <- fread(file.path("data", "results", "CCLE", "ccle_sum", "DA_luminal_vs_basal_table_s2", "DA_luminal_vs_basal_limma.csv"))
old_dt <- fread(file.path("data", "results", "CCLE", "DA_luminal_vs_basal", "DA_luminal_vs_basal_limma.csv"))
new_dt <- fread(file.path("data", "results", "CCLE_corrected", "DA_luminal_vs_basal", "DA_luminal_vs_basal_limma.csv"))

sum_dt[, id := strip_acc(Protein_Id)]
old_dt[, id := strip_acc(UniProtID)]
new_dt[, id := strip_acc(UniProtID)]
# One row per accession (isoforms map to same id — keep best adj.P)
sum_dt <- sum_dt[order(adj.P.Val)][, .SD[1], by = id]
old_dt <- old_dt[order(adj.P.Val)][, .SD[1], by = id]
new_dt <- new_dt[order(adj.P.Val)][, .SD[1], by = id]

summ <- function(nm, dt) {
  cat("\n=== ", nm, " ===\n", sep = "")
  cat("n rows:", nrow(dt), "\n")
  cat("adj.P < 0.05:", sum(dt$adj.P.Val < 0.05, na.rm = TRUE), "\n")
  cat("adj.P < 0.1:", sum(dt$adj.P.Val < 0.1, na.rm = TRUE), "\n")
  cat("P < 0.05:", sum(dt$P.Value < 0.05, na.rm = TRUE), "\n")
}

summ("Table S2 (ccle_sum)", sum_dt)
summ("Old MSstatsTMT gene_matrix (data/results/CCLE)", old_dt)
summ("CCLE_corrected gene_matrix", new_dt)

sig <- function(dt, thr = 0.05) dt[adj.P.Val < thr, id]
s05_sum <- sig(sum_dt)
s05_old <- sig(old_dt)
s05_new <- sig(new_dt)
cat("\n--- Overlap FDR 5% (UniProt accession, isoform suffix stripped) ---\n")
cat("Table S2 vs corrected: intersection", length(intersect(s05_sum, s05_new)), "\n")
cat("Old vs corrected: intersection", length(intersect(s05_old, s05_new)), "\n")
cat("Table S2 vs old: intersection", length(intersect(s05_sum, s05_old)), "\n")

m2 <- merge(
  merge(sum_dt[, .(id, logFC_s2 = logFC)], old_dt[, .(id, logFC_old = logFC)], by = "id"),
  new_dt[, .(id, logFC_new = logFC)],
  by = "id"
)
cat("\nCommon proteins (all three):", nrow(m2), "\n")
cat("cor(logFC) Table S2 vs corrected:", cor(m2$logFC_s2, m2$logFC_new, use = "pair"), "\n")
cat("cor(logFC) old vs corrected:", cor(m2$logFC_old, m2$logFC_new, use = "pair"), "\n")
cat("cor(logFC) Table S2 vs old:", cor(m2$logFC_s2, m2$logFC_old, use = "pair"), "\n")

out <- file.path("data", "results", "CCLE_corrected", "DA_luminal_vs_basal", "compare_to_table_s2_and_old_ccle.txt")
lines <- c(
  "Luminal vs Basal (4+4 lines) — three-way comparison",
  "==================================================",
  "",
  "Table S2: ccle_sum/.../Table_S2 derived protein matrix (paper Table S2, CCLE-code columns).",
  "Old: MSstatsTMT gene_matrix under data/results/CCLE/ (pre bridge/mixture fix).",
  "New: data/results/CCLE_corrected/gene_matrix.csv (corrected converter + MSstatsTMT).",
  "",
  "All use same limma contrast: Luminal - Basal; gene/protein filtering differs slightly per build.",
  ""
)
cat(paste(lines, collapse = "\n"), "\n")
sink(out)
summ("Table S2 (ccle_sum)", sum_dt)
summ("Old MSstatsTMT gene_matrix (data/results/CCLE)", old_dt)
summ("CCLE_corrected gene_matrix", new_dt)
cat("\n--- Overlap FDR 5% ---\n")
cat("Table S2 vs corrected:", length(intersect(s05_sum, s05_new)), "\n")
cat("Old vs corrected:", length(intersect(s05_old, s05_new)), "\n")
cat("Table S2 vs old:", length(intersect(s05_sum, s05_old)), "\n")
cat("\nCommon proteins (all three):", nrow(m2), "\n")
cat("cor(logFC) Table S2 vs corrected:", cor(m2$logFC_s2, m2$logFC_new, use = "pair"), "\n")
cat("cor(logFC) old vs corrected:", cor(m2$logFC_old, m2$logFC_new, use = "pair"), "\n")
cat("cor(logFC) Table S2 vs old:", cor(m2$logFC_s2, m2$logFC_old, use = "pair"), "\n")
sink()
message("Wrote ", out)
