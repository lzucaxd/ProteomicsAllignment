#!/usr/bin/env Rscript
# =============================================================================
# Master orchestrator — run full benchmark for all representations and tasks
# =============================================================================
# Usage:
#   Rscript scripts/benchmark/run_all_benchmarks.R \
#     --cptac_breast  data/results/PDC000120/gene_matrix.csv \
#     --cptac_lung    data/results/PDC000153/gene_matrix.csv \
#     --ccle           data/results/CCLE_corrected/gene_matrix.csv \
#     --subtype_map   data/results/PDC000120/gene_matrix_subtype_mapping.csv \
#     --ccle_sample   data/ccle_peptide/sample_info_ccle.csv \
#     --outdir        reports/benchmark_master
#     --subtype_union_meta  (optional) defaults to data/processed/union/sample_meta_breast_subtype.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ── Parse arguments ─────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) args[idx + 1] else default
}

CPTAC_BREAST  <- parse_arg("--cptac_breast",  "data/results/PDC000120/gene_matrix.csv")
CPTAC_LUNG    <- parse_arg("--cptac_lung",    "data/results/PDC000153/gene_matrix.csv")
CCLE_MATRIX   <- parse_arg("--ccle",          "data/results/CCLE_corrected/gene_matrix.csv")
SUBTYPE_MAP   <- parse_arg("--subtype_map",   "data/results/PDC000120/gene_matrix_subtype_mapping.csv")
CCLE_SAMPLE   <- parse_arg("--ccle_sample",   "data/ccle_peptide/sample_info_ccle.csv")
OUTDIR        <- parse_arg("--outdir",        "reports/benchmark_master")
SUBTYPE_UNION_META <- parse_arg("--subtype_union_meta", NA_character_)

# ── Source helpers ──────────────────────────────────────────────────────────
script_dir <- dirname(sys.frame(1)$ofile %||% "scripts/benchmark/run_all_benchmarks.R")
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
union_default <- file.path(repo_root, "data/processed/union/sample_meta_breast_subtype.csv")
union_alt <- file.path(repo_root, "data/processed/sample_meta_breast_subtype.csv")
if (is.na(SUBTYPE_UNION_META) || !nzchar(SUBTYPE_UNION_META)) {
  SUBTYPE_UNION_META <- if (file.exists(union_default)) union_default else union_alt
}
for (f in c("subset_strategies.R", "evaluation_helpers.R", "benchmark_runner.R",
            "task_breast_subtype.R", "task_breast_vs_lung.R",
            "diagnostics.R", "native_domain_da.R")) {
  fp <- file.path(script_dir, f)
  if (file.exists(fp)) source(fp) else warning("Missing: ", fp)
}

# =============================================================================
# STEP 1 — Build subsets
# =============================================================================
message("\n", strrep("=", 60), "\n  STEP 1: Building subsets\n", strrep("=", 60))

# Task A: Breast subtype
subtype_cptac <- build_subtype_subset_cptac(
  SUBTYPE_MAP, outdir = file.path(OUTDIR, "tasks", "breast_subtype")
)
ccle_sub_dir <- file.path(OUTDIR, "tasks", "breast_subtype")
if (file.exists(SUBTYPE_UNION_META)) {
  subtype_ccle <- build_subtype_subset_ccle(
    CCLE_MATRIX,
    union_meta_path = SUBTYPE_UNION_META,
    outdir = ccle_sub_dir
  )
} else {
  warning("Union subtype metadata not found (tried: ", union_default, " / ", union_alt,
          ") — using legacy 8-line CCLE subtype panel.")
  subtype_ccle <- build_subtype_subset_ccle(
    CCLE_MATRIX,
    union_meta_path = NULL,
    ccle_sample_info_path = CCLE_SAMPLE,
    outdir = ccle_sub_dir
  )
}

# Task B: Breast vs lung
bvl_cptac <- build_breast_vs_lung_subset_cptac(
  CPTAC_BREAST, CPTAC_LUNG,
  outdir = file.path(OUTDIR, "tasks", "breast_vs_lung")
)
bvl_ccle <- build_breast_vs_lung_subset_ccle(
  CCLE_SAMPLE, outdir = file.path(OUTDIR, "tasks", "breast_vs_lung")
)

# =============================================================================
# STEP 2 — Native-domain DA
# =============================================================================
message("\n", strrep("=", 60), "\n  STEP 2: Native-domain DA\n", strrep("=", 60))

native_dir <- file.path(OUTDIR, "native_domain_da")

# Task A
tryCatch({
  run_native_subtype_da_cptac(
    gene_matrix_path = CPTAC_BREAST,
    subtype_subset = subtype_cptac$subset,
    outdir = file.path(native_dir, "breast_subtype", "cptac")
  )
}, error = function(e) warning("Native CPTAC subtype DA failed: ", e$message))

tryCatch({
  run_native_subtype_da_ccle(
    gene_matrix_path = CCLE_MATRIX,
    subtype_subset = subtype_ccle,
    outdir = file.path(native_dir, "breast_subtype", "ccle")
  )
}, error = function(e) warning("Native CCLE subtype DA failed: ", e$message))

# Task B — CCLE breast vs lung
tryCatch({
  run_native_breast_vs_lung_da(
    domain = "CCLE",
    gene_matrix_path = CCLE_MATRIX,
    sample_subset = bvl_ccle$subset,
    outdir = file.path(native_dir, "breast_vs_lung", "ccle")
  )
}, error = function(e) warning("Native CCLE breast vs lung DA failed: ", e$message))

# Task B — CPTAC breast vs lung requires merged matrix
tryCatch({
  gm_breast <- fread(CPTAC_BREAST, header = TRUE)
  gm_lung <- fread(CPTAC_LUNG, header = TRUE)
  shared_genes <- intersect(gm_breast$GeneSymbol, gm_lung$GeneSymbol)
  message("  CPTAC breast+lung: ", length(shared_genes), " shared genes")

  breast_cols <- setdiff(names(gm_breast), c("GeneSymbol", "UniProtID"))
  lung_cols <- setdiff(names(gm_lung), c("GeneSymbol", "UniProtID"))

  merged <- merge(
    gm_breast[GeneSymbol %in% shared_genes, c("GeneSymbol", breast_cols), with = FALSE],
    gm_lung[GeneSymbol %in% shared_genes, c("GeneSymbol", lung_cols), with = FALSE],
    by = "GeneSymbol"
  )
  mat_merged <- as.matrix(merged[, -1])
  rownames(mat_merged) <- merged$GeneSymbol

  run_native_breast_vs_lung_da(
    domain = "CPTAC",
    gene_matrix_path = NULL,
    sample_subset = bvl_cptac$subset,
    outdir = file.path(native_dir, "breast_vs_lung", "cptac")
  )
}, error = function(e) warning("Native CPTAC breast vs lung DA failed: ", e$message))

# =============================================================================
# STEP 3 — Representation-level DA for each method
# =============================================================================
message("\n", strrep("=", 60), "\n  STEP 3: Representation-level DA\n", strrep("=", 60))

# ── Helper: load and prepare a combined matrix ──────────────────────────────
load_raw_combined_matrix <- function(cptac_path, ccle_path) {
  gm_cptac <- fread(cptac_path, header = TRUE)
  gm_ccle <- fread(ccle_path, header = TRUE)
  shared <- intersect(gm_cptac$GeneSymbol, gm_ccle$GeneSymbol)
  message("  Shared genes: ", length(shared))

  cptac_cols <- setdiff(names(gm_cptac), c("GeneSymbol", "UniProtID"))
  ccle_cols <- setdiff(names(gm_ccle), c("GeneSymbol", "UniProtID"))

  merged <- merge(
    gm_cptac[GeneSymbol %in% shared, c("GeneSymbol", cptac_cols), with = FALSE],
    gm_ccle[GeneSymbol %in% shared, c("GeneSymbol", ccle_cols), with = FALSE],
    by = "GeneSymbol"
  )
  mat <- as.matrix(merged[, -1])
  rownames(mat) <- merged$GeneSymbol
  list(matrix = mat, cptac_cols = cptac_cols, ccle_cols = ccle_cols)
}

# ── Build sample metadata for each task ─────────────────────────────────────
build_subtype_meta <- function(cptac_cols, ccle_cols, subtype_cptac, subtype_ccle) {
  sc <- as.data.table(subtype_cptac)
  sid_col <- if ("matrix_sample_id" %in% names(sc)) "matrix_sample_id" else "sample_id"
  cptac_meta <- sc[get(sid_col) %in% cptac_cols, .(
    sample_id = get(sid_col),
    domain = "CPTAC",
    condition = subtype
  )]
  if (!"mixture" %in% names(sc)) cptac_meta[, mixture := NA_character_]
  else cptac_meta[, mixture := sc[match(cptac_meta$sample_id, sc[[sid_col]]), mixture]]

  # sample_id is already the gene-matrix column name (union build)
  ccle_meta <- data.table(
    sample_id = subtype_ccle$sample_id,
    domain = "CCLE",
    condition = subtype_ccle$subtype,
    mixture = NA_character_
  )
  ccle_meta <- ccle_meta[sample_id %in% ccle_cols]

  rbind(cptac_meta, ccle_meta, fill = TRUE)
}

build_bvl_meta <- function(cptac_cols, ccle_cols, bvl_cptac, bvl_ccle) {
  cptac_meta <- as.data.table(bvl_cptac)[sample_id %in% cptac_cols, .(
    sample_id, domain = "CPTAC", condition = cancer_type
  )]
  ccle_meta <- as.data.table(bvl_ccle)[sample_id %in% ccle_cols, .(
    sample_id, domain = "CCLE", condition = cancer_type
  )]
  rbind(cptac_meta, ccle_meta, fill = TRUE)
}

# ── Run benchmark for raw representation ────────────────────────────────────
run_representation <- function(rep_name, matrix, sample_meta_subtype,
                                sample_meta_bvl, feature_meta = NULL) {
  message("\n--- Representation: ", rep_name, " ---")

  # Task A: Breast subtype
  tryCatch({
    run_benchmark(
      matrix = matrix,
      sample_meta = sample_meta_subtype,
      feature_meta = feature_meta,
      representation_name = rep_name,
      task_name = "breast_subtype"
    )
    generate_diagnostics(
      matrix = matrix,
      sample_meta = sample_meta_subtype,
      representation_name = rep_name,
      task_name = "breast_subtype",
      outdir = file.path(OUTDIR, "representation_level_da", rep_name, "breast_subtype"),
      marker_genes = c("FOXA1", "GATA3", "KRT5", "KRT14", "KRT17",
                       "EGFR", "ESR1", "PGR", "ERBB2", "CDH1")
    )
  }, error = function(e) warning(rep_name, " breast_subtype failed: ", e$message))

  # Task B: Breast vs lung
  tryCatch({
    run_benchmark(
      matrix = matrix,
      sample_meta = sample_meta_bvl,
      feature_meta = feature_meta,
      representation_name = rep_name,
      task_name = "breast_vs_lung"
    )
    generate_diagnostics(
      matrix = matrix,
      sample_meta = sample_meta_bvl,
      representation_name = rep_name,
      task_name = "breast_vs_lung",
      outdir = file.path(OUTDIR, "representation_level_da", rep_name, "breast_vs_lung"),
      marker_genes = c("NKX2-1", "SFTPB", "SFTPC", "NAPSA",
                       "GATA3", "FOXA1", "ESR1", "KRT19",
                       "EGFR", "ERBB2", "CDH1", "VIM")
    )
  }, error = function(e) warning(rep_name, " breast_vs_lung failed: ", e$message))
}

# ---- Raw representation ----
raw <- load_raw_combined_matrix(CPTAC_BREAST, CCLE_MATRIX)
meta_subtype <- build_subtype_meta(raw$cptac_cols, raw$ccle_cols,
                                    subtype_cptac$subset, subtype_ccle)
meta_bvl <- build_bvl_meta(raw$cptac_cols, raw$ccle_cols,
                            bvl_cptac$subset, bvl_ccle$subset)
run_representation("raw", raw$matrix, meta_subtype, meta_bvl)

# ---- Bridge-aware representation ----
bridge_matrix_path <- file.path("data", "results", "bridge_aware", "combined_matrix.csv")
if (file.exists(bridge_matrix_path)) {
  message("\nLoading bridge-aware matrix...")
  ba <- fread(bridge_matrix_path, header = TRUE)
  ba_mat <- as.matrix(ba[, -1])
  rownames(ba_mat) <- ba[[1]]
  run_representation("bridge_aware", ba_mat, meta_subtype, meta_bvl)
} else {
  warning("Bridge-aware matrix not found at ", bridge_matrix_path,
          ". Skipping. Run scripts/methods/run_bridge_aware_representation.R first.")
}

# ---- Celligner representation ----
celligner_matrix_path <- file.path("data", "results", "celligner", "combined_output.csv")
if (file.exists(celligner_matrix_path)) {
  message("\nLoading Celligner matrix...")
  cl <- fread(celligner_matrix_path, header = TRUE)
  cl_mat <- as.matrix(cl[, -1])
  rownames(cl_mat) <- cl[[1]]
  run_representation("celligner", cl_mat, meta_subtype, meta_bvl)
} else {
  warning("Celligner matrix not found at ", celligner_matrix_path,
          ". Skipping. Run scripts/methods/run_celligner_representation.py first.")
}

# =============================================================================
# STEP 4 — Aggregate comparison table
# =============================================================================
message("\n", strrep("=", 60), "\n  STEP 4: Aggregating comparison table\n", strrep("=", 60))

agreement_files <- Sys.glob(file.path(OUTDIR, "representation_level_da",
                                       "*", "*", "cross_domain_agreement.csv"))
if (length(agreement_files) > 0) {
  all_agreement <- rbindlist(lapply(agreement_files, fread))
  fwrite(all_agreement, file.path(OUTDIR, "benchmark_method_comparison.tsv"), sep = "\t")
  message("Comparison table: ", file.path(OUTDIR, "benchmark_method_comparison.tsv"))
  print(all_agreement)
} else {
  message("No cross-domain agreement files found yet.")
}

# =============================================================================
# STEP 5 — Summary
# =============================================================================
message("\n", strrep("=", 60))
message("  BENCHMARK COMPLETE")
message(strrep("=", 60))
message("Outputs: ", OUTDIR)
message("  tasks/                         ← subset definitions")
message("  native_domain_da/              ← MSstatsTMT / native limma results")
message("  representation_level_da/       ← per-method limma results + diagnostics")
message("  benchmark_method_comparison.tsv ← side-by-side metrics")
