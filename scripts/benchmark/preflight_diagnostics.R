#!/usr/bin/env Rscript
# Step 0: Pre-flight diagnostics (union matrices). ~1–2 min.
# Usage: Rscript preflight_diagnostics.R [--repo-root DIR] [--processed-dir data/processed_union]

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
proc_dir_rel <- "data/processed_union"
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--processed-dir" && i < length(args)) proc_dir_rel <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)
proc_dir <- normalizePath(file.path(repo_root, proc_dir_rel), mustWork = FALSE)
diag_dir <- file.path(repo_root, "reports/benchmark_master/diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

audit_gene_coverage <- function(matrix_path, annotation_path, task_name, cov_threshold = 0.30) {
  mat <- as.matrix(fread(matrix_path), rownames = 1)
  ann <- fread(annotation_path)
  if (nrow(mat) < ncol(mat)) mat <- t(mat)

  cptac_samples <- ann[toupper(domain) == "CPTAC", sample_id]
  ccle_samples  <- ann[toupper(domain) == "CCLE", sample_id]
  cptac_cols <- intersect(colnames(mat), cptac_samples)
  ccle_cols  <- intersect(colnames(mat), ccle_samples)

  cptac_coverage <- if (length(cptac_cols)) rowMeans(!is.na(mat[, cptac_cols, drop = FALSE])) else rep(0, nrow(mat))
  ccle_coverage  <- if (length(ccle_cols)) rowMeans(!is.na(mat[, ccle_cols, drop = FALSE])) else rep(0, nrow(mat))

  gene_audit <- data.table(
    gene = rownames(mat),
    cptac_coverage = round(cptac_coverage, 3),
    ccle_coverage  = round(ccle_coverage, 3)
  )
  gene_audit[, category := fcase(
    cptac_coverage > cov_threshold & ccle_coverage > cov_threshold, "both_domains",
    cptac_coverage > cov_threshold & ccle_coverage <= cov_threshold, "cptac_only",
    cptac_coverage <= cov_threshold & ccle_coverage > cov_threshold, "ccle_only",
    default = "low_coverage_both"
  )]

  cat("\n=== Gene Coverage Audit:", task_name, "(threshold", cov_threshold, "per domain) ===\n")
  print(table(gene_audit$category))
  cat("Total genes:", nrow(gene_audit), "\n")
  cat("Genes in both_domains:", sum(gene_audit$category == "both_domains"), "\n\n")

  out_csv <- file.path(diag_dir, paste0("gene_coverage_audit_", task_name, ".csv"))
  fwrite(gene_audit, out_csv)
  cat("Saved:", out_csv, "\n")
  invisible(gene_audit)
}

trace_ccle_subtype <- function(meta_path, out_txt, repo_root) {
  ann <- fread(meta_path)
  ccle_sub <- ann[toupper(domain) == "CCLE" & tolower(condition) %in% c("basal", "luminal")]
  lines <- c(
    "CCLE subtype samples (Basal/Luminal):",
    capture.output(print(ccle_sub[, .(sample_id, condition)])),
    "",
    paste0("n CCLE subtype samples: ", nrow(ccle_sub)),
    "",
    "Expected line name check (partial match):"
  )
  v2 <- file.path(repo_root, "data/ccle/ccle_breast_subtype_annotations_v2.csv")
  if (file.exists(v2)) {
    v2t <- fread(v2)
    if (!"BvL_group" %in% names(v2t)) {
      expected <- c("HCC70", "HCC1806", "HCC1143", "MDA-MB-468",
                    "CAMA-1", "MCF7", "T-47D", "ZR-75-1")
    } else {
      v2t <- v2t[tolower(BvL_group) %in% c("basal", "luminal")]
      expected <- unique(v2t$cell_line)
    }
    lines <- c(lines, paste("(from", v2, ")"))
  } else {
    expected <- c("HCC70", "HCC1806", "HCC1143", "MDA-MB-468",
                  "CAMA-1", "MCF7", "T-47D", "ZR-75-1")
  }
  for (line in expected) {
    found <- grep(line, ccle_sub$sample_id, ignore.case = TRUE, value = TRUE)
    if (length(found) == 0) lines <- c(lines, sprintf("MISSING: %s", line))
    else lines <- c(lines, sprintf("FOUND: %s -> %s", line, paste(found, collapse = ", ")))
  }
  cat(paste(lines, collapse = "\n"), "\n")
  writeLines(lines, out_txt)
  cat("Saved:", out_txt, "\n")
}

check_markers_file <- function(matrix_path, task_name, markers_named) {
  mat <- as.matrix(fread(matrix_path), rownames = 1)
  if (nrow(mat) < ncol(mat)) mat <- t(mat)
  genes <- rownames(mat)
  out_lines <- character()
  push <- function(...) out_lines <<- c(out_lines, sprintf(...))
  push("\n=== Marker Check: %s ===", task_name)
  for (i in seq_along(markers_named)) {
    m <- names(markers_named)[i]
    note <- markers_named[i]
    if (m %in% genes) {
      push("FOUND: %s (%s)", m, note)
    } else {
      cand <- grep(gsub("-", ".", fixed = TRUE, m), genes, ignore.case = TRUE, value = TRUE)
      if (length(cand) == 0) cand <- grep(m, genes, ignore.case = TRUE, value = TRUE)
      if (length(cand) > 0) push("ALIAS?: %s -> %s", m, paste(cand, collapse = ", "))
      else push("MISSING: %s (%s)", m, note)
    }
  }
  cat(paste(out_lines, collapse = "\n"), "\n")
  writeLines(out_lines, file.path(diag_dir, paste0("marker_presence_", task_name, ".txt")))
}

bvl_composition <- function(meta_path, out_csv) {
  ann <- fread(meta_path)
  tb <- table(ann$domain, ann$condition)
  con_txt <- file.path(diag_dir, "bvl_composition_console.txt")
  zz <- file(con_txt, open = "wt")
  writeLines(c("\n=== BvL Sample Composition ===", capture.output(print(tb))), zz)
  ccle <- ann[toupper(domain) == "CCLE"]
  writeLines(c(
    "",
    paste("CCLE Breast:", sum(ccle$condition == "Breast", na.rm = TRUE)),
    paste("CCLE Lung:", sum(ccle$condition == "Lung", na.rm = TRUE))
  ), zz)
  close(zz)
  cat("\n=== BvL Sample Composition ===\n")
  print(tb)
  cat("\nCCLE Breast:", sum(ccle$condition == "Breast", na.rm = TRUE), "\n")
  cat("CCLE Lung:", sum(ccle$condition == "Lung", na.rm = TRUE), "\n")
  tbdf <- as.data.frame.matrix(tb)
  tbdf$domain <- rownames(tbdf)
  fwrite(as.data.table(tbdf), out_csv)
  cat("Saved:", out_csv, con_txt, "\n")
}

cat("=== PREFLIGHT DIAGNOSTICS ===\n")
cat("Processed dir:", proc_dir, "\n")

if (!dir.exists(proc_dir)) stop("Processed dir not found: ", proc_dir)

subtype_mat <- file.path(proc_dir, "shared_gene_matrix_breast_subtype.csv")
subtype_meta <- file.path(proc_dir, "sample_meta_breast_subtype.csv")
bvl_mat <- file.path(proc_dir, "shared_gene_matrix_breast_vs_lung.csv")
bvl_meta <- file.path(proc_dir, "sample_meta_breast_vs_lung.csv")

if (file.exists(subtype_mat) && file.exists(subtype_meta)) {
  audit_gene_coverage(subtype_mat, subtype_meta, "breast_subtype", 0.30)
  trace_ccle_subtype(subtype_meta, file.path(diag_dir, "ccle_sample_trace.txt"), repo_root)
  check_markers_file(subtype_mat, "breast_subtype", c(
    ESR1 = "up_luminal", PGR = "up_luminal", GATA3 = "up_luminal", FOXA1 = "up_luminal",
    EGFR = "up_basal", KRT5 = "up_basal", KRT17 = "up_basal", FOXC1 = "up_basal"
  ))
}

if (file.exists(bvl_mat) && file.exists(bvl_meta)) {
  audit_gene_coverage(bvl_mat, bvl_meta, "breast_vs_lung", 0.30)
  bvl_composition(bvl_meta, file.path(diag_dir, "bvl_composition.csv"))
  check_markers_file(bvl_mat, "breast_vs_lung", c(
    GATA3 = "up_breast", FOXA1 = "up_breast", ESR1 = "up_breast",
    "NKX2-1" = "up_lung", SFTPB = "up_lung", NAPSA = "up_lung"
  ))
}

cat("\nPreflight complete. Outputs in:", diag_dir, "\n")
