#!/usr/bin/env Rscript
# Build reports/subtype_marker_panel.csv from benchmark DA tables (reproducible).
suppressPackageStartupMessages(library(data.table))

repo <- getwd()
if (!file.exists(file.path(repo, "data", "results", "PDC000120", "DA_subtype_subset_runs"))) {
  repo <- normalizePath(file.path(getwd(), ".."))
}
if (!file.exists(file.path(repo, "data", "results", "PDC000120", "DA_subtype_subset_runs"))) {
  repo <- normalizePath(file.path(getwd(), "..", ".."))
}

path_cptac <- file.path(repo, "data/results/PDC000120/DA_subtype_subset_runs/DA_MSstatsTMT_Luminal_vs_Basal.csv")
path_ccle_s2 <- file.path(repo, "data/results/CCLE/ccle_sum/DA_luminal_vs_basal_table_s2/DA_luminal_vs_basal_limma.csv")
path_ccle_cor <- file.path(repo, "data/results/CCLE_corrected/DA_luminal_vs_basal/DA_luminal_vs_basal_limma.csv")
stopifnot(file.exists(path_cptac))
stopifnot(file.exists(path_ccle_cor))
has_s2 <- file.exists(path_ccle_s2)

extract_gene_uniprot <- function(uid) {
  uid <- as.character(uid)
  g <- sub("^[^|]*\\|[^|]*\\|([A-Z0-9]+)(?:-[0-9]+)?_HUMAN.*", "\\1", uid, perl = TRUE)
  ifelse(grepl("^[A-Z0-9]+$", g), g, NA_character_)
}

cpt <- fread(path_cptac, showProgress = FALSE)
setDT(cpt)
adj <- if ("adj.pvalue" %in% names(cpt)) "adj.pvalue" else "adj.P.Val"
cpt[, gene := trimws(as.character(Gene_symbol))]
cpt <- cpt[nzchar(gene)][order(get(adj))][, .SD[1], by = gene]

load_ccle_dedup <- function(path) {
  cc <- fread(path, showProgress = FALSE)
  setDT(cc)
  id_col <- if ("Protein_Id" %in% names(cc)) "Protein_Id" else "UniProtID"
  cc[, gene := extract_gene_uniprot(get(id_col))]
  cc[!is.na(gene)][order(adj.P.Val)][, .SD[1], by = gene]
}
cc_s2 <- if (has_s2) load_ccle_dedup(path_ccle_s2) else NULL
cc_cor <- load_ccle_dedup(path_ccle_cor)

markers <- data.table(
  Gene = c("ESR1", "GATA3", "FOXA1", "PGR", "KRT18", "KRT5", "KRT14", "KRT17", "EGFR", "FOXC1"),
  Expected_subtype = c(rep("Luminal", 5L), rep("Basal", 5L))
)

interpret_one <- function(exp, logfc, sig) {
  if (is.na(logfc)) return("absent from table after dedup")
  dir_ok <- (exp == "Luminal" && logfc > 0) || (exp == "Basal" && logfc < 0)
  if (!sig) {
    if (dir_ok) return("direction OK; not FDR-significant")
    return("not FDR-significant; direction wrong or near null")
  }
  if (dir_ok) return("direction OK; FDR < 0.05")
  return("FDR < 0.05; direction opposite to expectation")
}

markers[, `:=`(
  CPTAC_logFC = NA_real_,
  CPTAC_FDR = NA_real_,
  CPTAC_sig = NA_character_,
  CCLE_logFC = NA_real_,
  CCLE_FDR = NA_real_,
  CCLE_sig = NA_character_,
  CCLE_data_source = NA_character_
)]

for (i in seq_len(nrow(markers))) {
  g <- markers$Gene[i]
  r1 <- cpt[gene == g]
  r2s2 <- if (!is.null(cc_s2)) cc_s2[gene == g] else NULL
  r2cor <- cc_cor[gene == g]
  if (nrow(r1)) {
    markers[i, `:=`(
      CPTAC_logFC = r1$log2FC[1],
      CPTAC_FDR = r1[[adj]][1],
      CPTAC_sig = fifelse(r1[[adj]][1] < 0.05, "yes", "no")
    )]
  }
  if (!is.null(r2s2) && nrow(r2s2)) {
    markers[i, `:=`(
      CCLE_logFC = r2s2$logFC[1],
      CCLE_FDR = r2s2$adj.P.Val[1],
      CCLE_sig = fifelse(r2s2$adj.P.Val[1] < 0.05, "yes", "no"),
      CCLE_data_source = "CCLE_Table_S2"
    )]
  } else if (nrow(r2cor)) {
    markers[i, `:=`(
      CCLE_logFC = r2cor$logFC[1],
      CCLE_FDR = r2cor$adj.P.Val[1],
      CCLE_sig = fifelse(r2cor$adj.P.Val[1] < 0.05, "yes", "no"),
      CCLE_data_source = "CCLE_corrected (fallback; gene absent in Table S2 limma)"
    )]
  } else {
    markers[i, CCLE_data_source := "absent (both Table S2 and CCLE_corrected)"]
  }
}

markers[, CPTAC_direction := fifelse(
  is.na(CPTAC_logFC), "absent",
  fifelse(CPTAC_logFC > 0, "higher in Luminal (Luminal - Basal > 0)", "higher in Basal (Luminal - Basal < 0)")
)]
markers[, CCLE_direction := fifelse(
  is.na(CCLE_logFC), "absent (not in CCLE limma after gene filter / mapping)",
  fifelse(CCLE_logFC > 0, "higher in Luminal (Luminal - Basal > 0)", "higher in Basal (Luminal - Basal < 0)")
)]
markers[is.na(CCLE_sig), CCLE_sig := "n/a"]

markers[, Short_interpretation := paste0(
  "CPTAC: ", mapply(interpret_one, Expected_subtype, CPTAC_logFC, CPTAC_sig == "yes"),
  " | CCLE: ", mapply(interpret_one, Expected_subtype, CCLE_logFC, CCLE_sig == "yes")
)]

out <- markers[, .(
  Gene,
  Expected_subtype,
  CPTAC_direction,
  CPTAC_sig,
  CCLE_direction,
  CCLE_sig,
  CCLE_data_source,
  Short_interpretation
)]
fwrite(out, file.path(repo, "reports/subtype_marker_panel.csv"))

meta <- c(
  "subtype_marker_panel.csv",
  "==============",
  paste("CPTAC source:", path_cptac),
  paste("CCLE Table S2 (preferred when gene present):", path_ccle_s2),
  paste("CCLE_corrected fallback:", path_ccle_cor),
  "Rule: use Table S2 limma row if gene maps; else CCLE_corrected limma (see CCLE_data_source).",
  "Contrast (both): Luminal − Basal; positive logFC/log2FC ⇒ higher in Luminal.",
  "FDR: adj.pvalue (CPTAC) / adj.P.Val (CCLE) < 0.05.",
  "Gene dedup: one row per gene (min FDR per gene).",
  ""
)
writeLines(meta, file.path(repo, "reports/subtype_marker_panel_meta.txt"))

txt <- c(
  capture.output(print(out, row.names = FALSE)),
  "",
  paste(meta, collapse = "\n")
)
writeLines(txt, file.path(repo, "reports/subtype_marker_panel.txt"))

message("Wrote reports/subtype_marker_panel.csv and .txt")
