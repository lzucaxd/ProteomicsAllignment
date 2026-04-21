#!/usr/bin/env Rscript
# =============================================================================
# Run Polished Marker Profile Plots (Meeting / Paper Quality)
# =============================================================================
# Sources the same data as run_sample_profile_plots.R but uses the redesigned
# polished plotting functions: faint sample traces + bold median, larger panels,
# max 4 markers per page, clear block shading and separators.
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
source(file.path(REPO, "scripts/benchmark/polished_profile_plots.R"))

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
gm_ccle   <- load_gm(file.path(REPO, "data/results/CCLE_corrected/gene_matrix.csv"))

# Lung: try to load if available
lung_path <- file.path(REPO, "data/results/PDC000153/gene_matrix.csv")
gm_lung <- tryCatch(load_gm(lung_path), error = function(e) NULL)
has_lung <- !is.null(gm_lung)
if (has_lung) cat("  Loaded CPTAC lung (PDC000153)\n") else cat("  CPTAC lung not yet available\n")

# Raw combined matrices
shared_bc <- intersect(rownames(gm_breast), rownames(gm_ccle))
raw_subtype_mat <- cbind(gm_breast[shared_bc, ], gm_ccle[shared_bc, ])

if (has_lung) {
  shared_blc <- Reduce(intersect, list(rownames(gm_breast), rownames(gm_lung), rownames(gm_ccle)))
  raw_bvl_mat <- cbind(gm_breast[shared_blc, ], gm_lung[shared_blc, ], gm_ccle[shared_blc, ])
} else {
  shared_blc <- intersect(rownames(gm_breast), rownames(gm_ccle))
  raw_bvl_mat <- cbind(gm_breast[shared_blc, ], gm_ccle[shared_blc, ])
}

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
sm <- fread(file.path(REPO, "data/results/PDC000120/gene_matrix_subtype_mapping.csv"))
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

meta_bvl_rows <- list(
  data.table(sample_id = colnames(gm_breast), condition = "Breast", domain = "CPTAC"),
  data.table(sample_id = ccle_breast_ids, condition = "Breast", domain = "CCLE"),
  data.table(sample_id = ccle_lung_ids, condition = "Lung", domain = "CCLE")
)
if (has_lung) {
  meta_bvl_rows <- c(list(meta_bvl_rows[[1]]),
    list(data.table(sample_id = colnames(gm_lung), condition = "Lung", domain = "CPTAC")),
    meta_bvl_rows[2:3])
}
meta_bvl <- rbindlist(meta_bvl_rows)
cat("  BvL:     ", nrow(meta_bvl), "samples\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 3: Define markers
# ═══════════════════════════════════════════════════════════════════════════

subtype_markers <- c("ESR1", "PGR", "GATA3", "FOXA1", "EGFR", "KRT5", "KRT17", "FOXC1")

bvl_markers <- c("NKX2-1", "SFTPB", "NAPSA", "GATA3", "FOXA1", "ESR1",
                  "EGFR", "ERBB2", "CDH1", "KRT7", "MUC1")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 4: Assemble method representations
# ═══════════════════════════════════════════════════════════════════════════
cat("\nAssembling representations...\n")

reps_subtype <- list(
  raw           = list(matrix = raw_subtype_mat, name = "Raw"),
  bridge_shift  = list(matrix = gm_bshift,       name = "Bridge Shift-Only"),
  bridge_scale  = list(matrix = gm_bscale,       name = "Bridge Shift+Scale"),
  celligner     = list(matrix = gm_cell,          name = "Celligner")
)

reps_bvl <- list(
  raw           = list(matrix = raw_bvl_mat,  name = "Raw"),
  bridge_shift  = list(matrix = gm_bshift,    name = "Bridge Shift-Only"),
  bridge_scale  = list(matrix = gm_bscale,    name = "Bridge Shift+Scale"),
  celligner     = list(matrix = gm_cell,       name = "Celligner")
)

# ═══════════════════════════════════════════════════════════════════════════
# STEP 5: Generate polished figures
# ═══════════════════════════════════════════════════════════════════════════

# ── Task A: Breast subtype ───────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("  TASK A: Breast Subtype — Polished Profiles\n")
cat(strrep("=", 60), "\n")

sub_pol <- file.path(OUTBASE, "breast_subtype/polished")
res_sub <- make_polished_sample_profile_plots(
  representations = reps_subtype,
  task_meta = meta_subtype,
  markers = subtype_markers,
  task_name = "breast subtype",
  outdir = sub_pol,
  block_order = c("CPTAC Luminal", "CCLE Luminal", "CPTAC Basal", "CCLE Basal"),
  markers_per_page = 4,
  fmt = c("png", "pdf")
)

# ── Task B: Breast vs lung ──────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("  TASK B: Breast vs Lung — Polished Profiles\n")
cat(strrep("=", 60), "\n")

bvl_pol <- file.path(OUTBASE, "breast_vs_lung/polished")
res_bvl <- make_polished_sample_profile_plots(
  representations = reps_bvl,
  task_meta = meta_bvl,
  markers = bvl_markers,
  task_name = "breast vs lung",
  outdir = bvl_pol,
  block_order = c("CPTAC Breast", "CCLE Breast", "CPTAC Lung", "CCLE Lung"),
  markers_per_page = 4,
  fmt = c("png", "pdf")
)

# ═══════════════════════════════════════════════════════════════════════════
# STEP 6: Move old grids to QC subdirectories
# ═══════════════════════════════════════════════════════════════════════════
cat("\nOrganising old grid plots as QC...\n")

move_to_qc <- function(task_dir) {
  qc_dir <- file.path(task_dir, "qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  old <- list.files(task_dir, pattern = "^(profile_grid|boxplot_summary)",
                    full.names = TRUE)
  for (f in old) {
    dest <- file.path(qc_dir, basename(f))
    file.rename(f, dest)
  }
  length(old)
}

n_moved_sub <- move_to_qc(file.path(OUTBASE, "breast_subtype"))
n_moved_bvl <- move_to_qc(file.path(OUTBASE, "breast_vs_lung"))
cat("  Moved", n_moved_sub, "old grid files to breast_subtype/qc/\n")
cat("  Moved", n_moved_bvl, "old grid files to breast_vs_lung/qc/\n")

# ═══════════════════════════════════════════════════════════════════════════
# STEP 7: Write polished plot README
# ═══════════════════════════════════════════════════════════════════════════

readme_lines <- c(
  "# Polished Marker Profile Plots",
  "",
  "## What Changed",
  "",
  "The original grid plots packed all markers into a single giant figure with tiny cells.",
  "The polished version splits markers into pages of 3-4 markers each, uses much larger",
  "panels, and applies a **faint trace + bold median** visual style:",
  "",
  "| Feature | Old grid | Polished |",
  "|---------|----------|----------|",
  "| Panels per page | 10+ markers × 4 methods | 3-4 markers × N methods |",
  "| Individual samples | small coloured dots | faint translucent dots |",
  "| Summary statistic | none | bold block-median segment |",
  "| Block separation | light background only | shaded bands + dashed separators |",
  "| Readability | too small for meetings | presentation/paper ready |",
  "",
  "## How to Interpret",
  "",
  "Each row is a marker gene. Each column is a method/representation.",
  "",
  "Within each panel:",
  "",
  "- **Faint dots** = individual sample abundance values, ordered by biological block.",
  "- **Bold horizontal segments** = block-wise median abundance. These are the",
  "  summary traces that matter most for comparing group behaviour across methods.",
  "- **Thin dashed connector** between block medians shows the overall cross-block",
  "  trend for the marker in that representation.",
  "- **Background shading** identifies the four biological blocks (e.g. CPTAC Luminal,",
  "  CPTAC Basal, CCLE Luminal, CCLE Basal).",
  "- **Vertical dashed lines** separate the blocks.",
  "",
  "### What \"Good\" Looks Like",
  "",
  "- For a known luminal marker (e.g. ESR1): the median segment for Luminal blocks should",
  "  be clearly higher than for Basal blocks, and this pattern should be consistent across",
  "  methods and across domains (CPTAC and CCLE).",
  "- A well-harmonised method should reduce the vertical gap between CPTAC and CCLE",
  "  medians while preserving the biological contrast.",
  "",
  "### What to Watch Out For",
  "",
  "- **Flat medians**: if all four block medians sit at the same level, the biological",
  "  contrast has been lost (bad for alignment).",
  "- **Large domain gap**: if CPTAC medians are systematically 5-10 units above CCLE",
  "  (or vice versa) regardless of biology, harmonisation has not helped.",
  "- **Panels marked N/A**: marker is absent from that method's matrix (often Celligner,",
  "  which uses a reduced gene set).",
  "",
  "## Output Structure",
  "",
  "```",
  "marker_profiles/",
  "  breast_subtype/",
  "    polished/           <- main polished figures (meeting/paper)",
  "    qc/                 <- old dense grids + boxplots (internal QC)",
  "  breast_vs_lung/",
  "    polished/           <- main polished figures",
  "    qc/                 <- old dense grids + boxplots",
  "```",
  "",
  "## Reproducibility",
  "",
  "Regenerate with:",
  "```bash",
  "Rscript --vanilla scripts/benchmark/run_polished_profile_plots.R",
  "```"
)

writeLines(readme_lines, file.path(OUTBASE, "polished_profile_plot_readme.md"))

# ═══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("  POLISHED PROFILE PLOTS COMPLETE\n")
cat(strrep("=", 60), "\n\n")

cat("Polished figures created:\n")
for (p in c(res_sub$paths, res_bvl$paths)) cat("  ", p, "\n")

cat("\nBreast subtype markers plotted:", paste(res_sub$plottable, collapse = ", "), "\n")
if (length(res_sub$skipped) > 0)
  cat("  Skipped:", paste(res_sub$skipped, collapse = ", "), "\n")

cat("Breast vs lung markers plotted:", paste(res_bvl$plottable, collapse = ", "), "\n")
if (length(res_bvl$skipped) > 0)
  cat("  Skipped:", paste(res_bvl$skipped, collapse = ", "), "\n")

cat("\nMethods:", paste(sapply(reps_subtype, `[[`, "name"), collapse = ", "), "\n")
cat("CPTAC lung data:", ifelse(has_lung, "included", "not yet available"), "\n")
