#!/usr/bin/env Rscript
# Step 0 (v2): Process ccle_breast_subtype_annotations_v2.csv — HER2 exclusion summary,
# CAL120 dedup, validation against CCLE gene matrix columns, disputed-line notes.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

raw_path <- file.path(repo_root, "data/ccle/ccle_breast_subtype_annotations_v2.csv")
out_proc <- file.path(repo_root, "data/processed/ccle_breast_subtype_annotation_processed.csv")
out_notes <- file.path(repo_root, "data/processed/ccle_breast_subtype_annotation_notes.txt")
out_val <- file.path(repo_root, "reports/benchmark_master/diagnostics/ccle_subtype_v2_validation.txt")
ccle_mat <- file.path(repo_root, "data/results/CCLE_corrected/gene_matrix.csv")

dir.create(dirname(out_proc), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_val), recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(raw_path))
ann_raw <- fread(raw_path)

cat("=== Raw annotation summary (BvL_group x subtype_detail) ===\n")
print(table(ann_raw$BvL_group, ann_raw$subtype_detail, useNA = "ifany"))

ann_lb <- ann_raw[BvL_group %in% c("Luminal", "Basal")]
cat("\nAfter HER2 exclusion (Luminal vs Basal task):\n")
print(table(ann_lb$BvL_group))

cal120 <- ann_lb[cell_line == "CAL120"]
cat("\nCAL120 entries (plexes to average in matrix if >1 column exists):\n")
print(cal120[, .(cell_line, plex, BvL_group, column_id)])

ann_dedup <- ann_lb[, .(
  BvL_group = BvL_group[1],
  subtype_detail = subtype_detail[1],
  plexes = paste(unique(plex), collapse = "+"),
  column_ids = paste(column_id, collapse = ";"),
  n_plexes = .N
), by = cell_line]

cat("\n=== Final cell line list (Basal + Luminal, HER2 excluded) ===\n")
cat("Luminal lines:", sum(ann_dedup$BvL_group == "Luminal"), "\n")
cat("Basal lines:", sum(ann_dedup$BvL_group == "Basal"), "\n")
cat("Total lines:", nrow(ann_dedup), "\n\n")
print(ann_dedup[order(BvL_group, cell_line)])

fwrite(ann_dedup, out_proc)
cat("\nWrote:", out_proc, "\n")

note_lines <- c(
  "Annotation notes (v2)",
  "--------------------",
  "HCC1500: Luminal (LA) per ER+ / ATCC / Kao2009; Neve2006 basalB is disputed — kept Luminal.",
  "MDAMB453: HER2 in BvL_group — excluded from Luminal vs Basal; included in BvL Breast group when using full v2 for CCLE breast.",
  "HER2 lines excluded from subtype contrast: AU565, HCC1954, HCC2218, JIMT1, MDAMB453 (per annotation file).",
  ""
)
writeLines(note_lines, out_notes)
cat("Wrote:", out_notes, "\n")

# --- 0b: validate column_id patterns against CCLE matrix header ---
ccle_cols <- strsplit(readLines(ccle_mat, n = 1L), ",", fixed = TRUE)[[1]]
ccle_cols <- trimws(ccle_cols)
ccle_cols <- ccle_cols[!tolower(ccle_cols) %in% c("genesymbol", "uniprotd", "gene")]

val_con <- file(out_val, open = "wt")
writeLines(c("=== CCLE subtype v2 — column_id validation ===", paste("Matrix:", ccle_mat), ""), val_con)

for (i in seq_len(nrow(ann_dedup))) {
  line <- ann_dedup$cell_line[i]
  ids <- strsplit(ann_dedup$column_ids[i], ";", fixed = TRUE)[[1]]
  ids <- trimws(ids)
  found_each <- vapply(ids, function(cid) any(grepl(cid, ccle_cols, fixed = TRUE)), logical(1L))
  if (all(found_each)) {
    writeLines(sprintf("FOUND: %s (column_id patterns OK)", line), val_con)
  } else {
    miss <- ids[!found_each]
    writeLines(sprintf("MISSING pattern: %s — %s", line, paste(miss, collapse = ", ")), val_con)
  }
}
close(val_con)
cat("Wrote:", out_val, "\n")
cat("\nprocess_ccle_annotations_v2.R complete.\n")
