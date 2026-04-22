#!/usr/bin/env Rscript
# =============================================================================
# Run Sample Profile Plots for All Current Benchmark Representations
# =============================================================================
# This script wires the method-agnostic plotting system to the current set of
# benchmark representations and tasks. To add a new method, add it to the
# `representations` list below — no changes to the plotting code needed.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

.local_args <- commandArgs(trailingOnly = FALSE)
.local_file <- .local_args[startsWith(.local_args, "--file=")]
if (length(.local_file)) {
  .local_bench <- dirname(normalizePath(sub("^--file=", "", .local_file[1L])))
} else {
  .local_bench <- normalizePath(file.path(getwd(), "scripts", "benchmark"), mustWork = FALSE)
}
source(file.path(.local_bench, "harmonize_paths.R"))
REPO <- harmonize_repo_root()
source(file.path(REPO, "scripts/benchmark/sample_profile_plots.R"))

OUTBASE <- file.path(REPO, "reports/benchmark_master/marker_profiles")

# ─── Helper: load gene matrix (genes×samples) ───────────────────────────
load_gm <- function(path) {
  dt <- fread(path, header = TRUE)
  id_cols <- intersect(c("GeneSymbol", "UniProtID", "Gene"), names(dt))
  scols <- setdiff(names(dt), id_cols)
  mat <- as.matrix(dt[, ..scols])
  rownames(mat) <- as.character(dt[[names(dt)[1]]])
  mat
}

# ═══════════════════════════════════════════════════════════════════════════
# STEP 1: Load representations
# ═══════════════════════════════════════════════════════════════════════════
cat("Loading representations...\n")

gm_breast <- load_gm(file.path(REPO, "data/results/PDC000120/gene_matrix.csv"))
gm_lung   <- load_gm(file.path(REPO, "data/results/PDC000153/gene_matrix.csv"))
gm_ccle   <- load_gm(file.path(REPO, "data/results/CCLE_corrected/gene_matrix.csv"))

# Raw: intersection of CPTAC breast + CCLE (subtype task)
shared_bc <- intersect(rownames(gm_breast), rownames(gm_ccle))
raw_subtype_mat <- cbind(gm_breast[shared_bc, ], gm_ccle[shared_bc, ])

# Raw: intersection of CPTAC breast + lung + CCLE (BvL task)
shared_blc <- Reduce(intersect, list(rownames(gm_breast), rownames(gm_lung), rownames(gm_ccle)))
raw_bvl_mat <- cbind(gm_breast[shared_blc, ], gm_lung[shared_blc, ], gm_ccle[shared_blc, ])

# Bridge-aware
gm_bshift <- load_gm(file.path(REPO, "reports/benchmark_master/methods/bridge_aware/bridge_aware_shift_only_matrix.csv"))
gm_bscale <- load_gm(file.path(REPO, "reports/benchmark_master/methods/bridge_aware/bridge_aware_shift_scale_matrix.csv"))

# Celligner
cell_dt <- fread(file.path(REPO, "reports/benchmark_master/celligner_all/celligner_aligned_matrix.csv"))
gm_cell <- t(as.matrix(cell_dt[, -1, with = FALSE]))
colnames(gm_cell) <- cell_dt[[1]]
rownames(gm_cell) <- names(cell_dt)[-1]

cat("  Raw (subtype):", ncol(raw_subtype_mat), "samples ×", nrow(raw_subtype_mat), "genes\n")
cat("  Raw (bvl):    ", ncol(raw_bvl_mat), "samples ×", nrow(raw_bvl_mat), "genes\n")
cat("  Bridge shift: ", ncol(gm_bshift), "samples ×", nrow(gm_bshift), "genes\n")
cat("  Bridge scale: ", ncol(gm_bscale), "samples ×", nrow(gm_bscale), "genes\n")
cat("  Celligner:    ", ncol(gm_cell), "samples ×", nrow(gm_cell), "genes\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 2: Build metadata for both tasks
# ═══════════════════════════════════════════════════════════════════════════
cat("\nBuilding task metadata...\n")

# ── Subtype metadata ─────────────────────────────────────────────────────
sm <- fread(file.path(REPO, "data/annotations/cptac/PDC000120/gene_matrix_subtype_mapping.csv"))
st_col <- if ("sample_type" %in% names(sm)) "sample_type" else "sample_type_if_available"

gm_cols <- colnames(gm_breast)
gm_cols_lower <- tolower(gm_cols)
spam <- sm[tolower(get(st_col)) == "sample" &
            tolower(pam50) %in% c("basal", "luma", "lumb") &
            exists_in_gene_matrix == TRUE]
spam[, subtype := ifelse(tolower(pam50) == "basal", "Basal", "Luminal")]
spam[, matched_col := {
  idx <- match(tolower(matrix_sample_id), gm_cols_lower)
  fifelse(is.na(idx), NA_character_, gm_cols[idx])
}, by = seq_len(nrow(spam))]
tumors <- unique(spam[!is.na(matched_col)], by = "matched_col")
cptac_sub <- tumors[, .(sample_id = matched_col, condition = subtype, domain = "CPTAC")]

ccle_basal <- c("HCC70", "HCC1806", "HCC1143", "MDA-MB-468")
ccle_luminal <- c("CAMA-1", "MCF7", "T-47D", "ZR-75-1")
match_ccle <- function(lines, cols) {
  out <- character(0)
  for (ln in lines) {
    pat <- gsub("-", ".", ln, fixed = TRUE)
    m <- grep(pat, cols, ignore.case = TRUE, value = TRUE)
    if (length(m) >= 1) out <- c(out, m[1])
  }
  out
}
ccle_b_ids <- match_ccle(ccle_basal, colnames(gm_ccle))
ccle_l_ids <- match_ccle(ccle_luminal, colnames(gm_ccle))
ccle_sub <- data.table(
  sample_id = c(ccle_b_ids, ccle_l_ids),
  condition = c(rep("Basal", length(ccle_b_ids)), rep("Luminal", length(ccle_l_ids))),
  domain = "CCLE"
)
meta_subtype <- rbind(cptac_sub[, .(sample_id, condition, domain)], ccle_sub)
cat("  Subtype: ", nrow(meta_subtype), "samples\n")

# ── Breast vs lung metadata ─────────────────────────────────────────────
ccle_info <- as.data.table(read.csv(
  file.path(REPO, "data/ccle_peptide/sample_info_ccle.csv"),
  stringsAsFactors = FALSE, fill = TRUE))
ccle_info <- ccle_info[nchar(Cell.Line) > 0]
setnames(ccle_info, "Cell.Line", "Cell Line", skip_absent = TRUE)
setnames(ccle_info, "Tissue.of.Origin", "Tissue of Origin", skip_absent = TRUE)

tissue_map <- setNames(ccle_info[["Tissue of Origin"]], ccle_info[["Cell Line"]])
ccle_cols <- colnames(gm_ccle)
normalize_name <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
ccle_col_tissue <- rep(NA_character_, length(ccle_cols))
for (i in seq_along(ccle_cols)) {
  col <- ccle_cols[i]
  if (col %in% names(tissue_map)) { ccle_col_tissue[i] <- tissue_map[[col]]; next }
  col_norm <- normalize_name(col)
  for (nm in names(tissue_map)) {
    if (normalize_name(nm) == col_norm) { ccle_col_tissue[i] <- tissue_map[[nm]]; break }
  }
}
names(ccle_col_tissue) <- ccle_cols
ccle_breast_ids <- ccle_cols[!is.na(ccle_col_tissue) & tolower(ccle_col_tissue) == "breast"]
ccle_lung_ids <- ccle_cols[!is.na(ccle_col_tissue) & tolower(ccle_col_tissue) == "lung"]

meta_bvl <- rbind(
  data.table(sample_id = colnames(gm_breast), condition = "Breast", domain = "CPTAC"),
  data.table(sample_id = colnames(gm_lung), condition = "Lung", domain = "CPTAC"),
  data.table(sample_id = ccle_breast_ids, condition = "Breast", domain = "CCLE"),
  data.table(sample_id = ccle_lung_ids, condition = "Lung", domain = "CCLE")
)
cat("  BvL:     ", nrow(meta_bvl), "samples\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 3: Define marker panels
# ═══════════════════════════════════════════════════════════════════════════

subtype_markers <- c("ESR1", "PGR", "GATA3", "FOXA1", "EGFR", "KRT5", "KRT17",
                      "FOXC1", "ERBB2", "CDH1", "KRT14")

bvl_markers <- c("NKX2-1", "SFTPB", "NAPSA", "TTF1",
                  "GATA3", "FOXA1", "ESR1", "KRT19",
                  "EGFR", "ERBB2", "CDH1", "VIM",
                  "KRT5", "KRT7", "MUC1")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4: Assemble method-agnostic representation collections
# ═══════════════════════════════════════════════════════════════════════════
cat("\nAssembling representations...\n")

# For subtype task: Raw uses intersection of CPTAC breast + CCLE
reps_subtype <- list(
  raw           = list(matrix = raw_subtype_mat, name = "Raw"),
  bridge_shift  = list(matrix = gm_bshift,       name = "Bridge Shift-Only"),
  bridge_scale  = list(matrix = gm_bscale,       name = "Bridge Shift+Scale"),
  celligner     = list(matrix = gm_cell,          name = "Celligner")
)

# For BvL task: Raw uses intersection of CPTAC breast + lung + CCLE
reps_bvl <- list(
  raw           = list(matrix = raw_bvl_mat,  name = "Raw"),
  bridge_shift  = list(matrix = gm_bshift,    name = "Bridge Shift-Only"),
  bridge_scale  = list(matrix = gm_bscale,    name = "Bridge Shift+Scale"),
  celligner     = list(matrix = gm_cell,       name = "Celligner")
)

# ═══════════════════════════════════════════════════════════════════════════
# STEP 5: Audit marker availability
# ═══════════════════════════════════════════════════════════════════════════
cat("\nAuditing marker availability...\n")

audit_markers <- function(markers, reps, label) {
  rows <- list()
  for (g in markers) {
    for (nm in names(reps)) {
      rows[[length(rows) + 1]] <- data.table(
        task = label, marker = g, method = nm,
        present = g %in% rownames(reps[[nm]]$matrix)
      )
    }
  }
  rbindlist(rows)
}

audit_sub <- audit_markers(subtype_markers, reps_subtype, "breast_subtype")
audit_bvl <- audit_markers(bvl_markers, reps_bvl, "breast_vs_lung")
audit_all <- rbind(audit_sub, audit_bvl)

# Which markers are present in at least one method?
sub_present <- audit_sub[, .(any_present = any(present)), by = marker]
bvl_present <- audit_bvl[, .(any_present = any(present)), by = marker]

cat("  Subtype markers available in >= 1 method:",
    sum(sub_present$any_present), "/", nrow(sub_present), "\n")
cat("  BvL markers available in >= 1 method:",
    sum(bvl_present$any_present), "/", nrow(bvl_present), "\n")

# Filter to plottable markers
subtype_plot <- sub_present[any_present == TRUE, marker]
bvl_plot <- bvl_present[any_present == TRUE, marker]

# Save marker selection summary
selection_dir <- file.path(OUTBASE)
dir.create(selection_dir, recursive = TRUE, showWarnings = FALSE)

sel_lines <- c(
  "# Marker Selection Summary",
  paste0("\nGenerated: ", Sys.time()),
  "",
  "## Breast Subtype Markers",
  "",
  sprintf("| Marker | %s |",
          paste(names(reps_subtype), collapse = " | "))
)
for (g in subtype_markers) {
  row <- sprintf("| %s |", g)
  for (nm in names(reps_subtype)) {
    p <- audit_sub[marker == g & method == nm, present]
    row <- paste0(row, ifelse(p, " present |", " **ABSENT** |"))
  }
  sel_lines <- c(sel_lines, row)
}
sel_lines <- c(sel_lines, "",
  paste0("Plotted: ", paste(subtype_plot, collapse = ", ")),
  paste0("Skipped: ", paste(setdiff(subtype_markers, subtype_plot), collapse = ", ")),
  "",
  "## Breast vs Lung Markers",
  "",
  sprintf("| Marker | %s |",
          paste(names(reps_bvl), collapse = " | "))
)
for (g in bvl_markers) {
  row <- sprintf("| %s |", g)
  for (nm in names(reps_bvl)) {
    p <- audit_bvl[marker == g & method == nm, present]
    row <- paste0(row, ifelse(p, " present |", " **ABSENT** |"))
  }
  sel_lines <- c(sel_lines, row)
}
sel_lines <- c(sel_lines, "",
  paste0("Plotted: ", paste(bvl_plot, collapse = ", ")),
  paste0("Skipped: ", paste(setdiff(bvl_markers, bvl_plot), collapse = ", "))
)
writeLines(sel_lines, file.path(selection_dir, "marker_selection_summary.md"))

# ═══════════════════════════════════════════════════════════════════════════
# STEP 6: Generate plots
# ═══════════════════════════════════════════════════════════════════════════

# ── Task A: Breast subtype ───────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("  TASK A: Breast Subtype Profile Plots\n")
cat(strrep("=", 60), "\n")

sub_outdir <- file.path(OUTBASE, "breast_subtype")
res_sub <- make_sample_profile_plots(
  representations = reps_subtype,
  task_meta = meta_subtype,
  markers = subtype_plot,
  task_name = "breast_subtype",
  outdir = sub_outdir,
  block_order = c("CPTAC Luminal", "CPTAC Basal", "CCLE Luminal", "CCLE Basal"),
  fmt = c("png", "pdf")
)

# ── Task B: Breast vs lung ──────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("  TASK B: Breast vs Lung Profile Plots\n")
cat(strrep("=", 60), "\n")

bvl_outdir <- file.path(OUTBASE, "breast_vs_lung")
res_bvl <- make_sample_profile_plots(
  representations = reps_bvl,
  task_meta = meta_bvl,
  markers = bvl_plot,
  task_name = "breast_vs_lung",
  outdir = bvl_outdir,
  block_order = c("CPTAC Breast", "CPTAC Lung", "CCLE Breast", "CCLE Lung"),
  fmt = c("png", "pdf")
)

# ═══════════════════════════════════════════════════════════════════════════
# STEP 7: Save audit
# ═══════════════════════════════════════════════════════════════════════════
fwrite(audit_all, file.path(OUTBASE, "marker_availability_audit.csv"))

# ═══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("  SAMPLE PROFILE PLOTS COMPLETE\n")
cat(strrep("=", 60), "\n\n")

cat("Files created:\n")
cat("  scripts/benchmark/sample_profile_plots.R (generic plotting library)\n")
cat("  scripts/benchmark/run_sample_profile_plots.R (this runner)\n")
cat("  reports/benchmark_master/marker_profiles/marker_selection_summary.md\n")
cat("  reports/benchmark_master/marker_profiles/marker_availability_audit.csv\n")

n_sub <- length(res_sub$grid$plottable)
n_bvl <- length(res_bvl$grid$plottable)
cat("\nPlots generated (combined figures, not per-marker files):\n")
cat("  Breast subtype:", n_sub, "markers ×", length(reps_subtype), "methods → 2 files (grid + boxplot)\n")
cat("  Breast vs lung:", n_bvl, "markers ×", length(reps_bvl), "methods → 2 files (grid + boxplot)\n")

cat("\nMethods handled:", length(reps_subtype), "\n")
cat("  ", paste(sapply(reps_subtype, `[[`, "name"), collapse = ", "), "\n")

cat("\nOutput locations:\n")
cat("  ", sub_outdir, "\n")
cat("  ", bvl_outdir, "\n")

cat("\nNote: Bridge matrices currently only contain CPTAC PDC000120 (breast).\n")
cat("  BvL bridge plots will show CCLE samples only for the CPTAC Lung block.\n")
cat("  Once PDC000153 lung data is processed, re-run to include CPTAC lung.\n")
