#!/usr/bin/env Rscript
# =============================================================================
# Subtype-aware differential abundance — CPTAC breast (PDC000120) using MSstatsTMT
# Tumor-only, PAM50-labeled samples; contrast e.g. Basal vs LumA.
# =============================================================================
# Inputs:
#   - results/PDC000120/protein_summary.tsv (from MSstatsTMT proteinSummarization)
#   - results/PDC000120/DA_subtype_tumor_only.csv
#   - results/PDC000120/annotation_filled_corrected.csv (optional, for mapping)
#
# Usage (from project root or data/):
#   Rscript data/scripts/DA_subtype_MSstatsTMT_PDC000120.R
#   # or: cd data && Rscript --vanilla scripts/DA_subtype_MSstatsTMT_PDC000120.R
#
# Optional environment variables:
#   PDC_SUBTYPE_ANNOT — path to annotation CSV (default: results/PDC000120/DA_subtype_tumor_only.csv)
#   PDC_MSSTATSTMT_OUT_DIR — where to write DA outputs (default: same as results/PDC000120)
# (protein_summary.tsv can be large; loading may take a minute.)
#
# Outputs (names use contrast_name, e.g. Luminal_vs_Basal):
#   - results/PDC000120/DA_MSstatsTMT_<contrast>.csv
#   - results/PDC000120/volcano_MSstatsTMT_<contrast>.pdf
#   - results/PDC000120/DA_MSstatsTMT_<contrast>_summary.txt
#   - results/PDC000120/<contrast>_mixture_balance.csv
# =============================================================================

# -----------------------------------------------------------------------------
# CONFIG — change for other contrasts
# -----------------------------------------------------------------------------
# Primary contrast: Luminal - Basal (positive log2FC => higher in Luminal).
subtype_a   <- "Luminal"   # numerator in estimand; LumA + LumB pooled to "Luminal"
subtype_b   <- "Basal"
contrast_name <- "Luminal_vs_Basal"

# Coverage filter: require non-missing in at least this fraction of samples (before DA)
pct_overall_min <- 0.35   # 35% of all samples (30–40%)
pct_group_min   <- 0.25   # 25% of each subtype group (20–30%)

# For other contrasts, set e.g.:
# subtype_a <- "Basal";  subtype_b <- "LumA"; contrast_name <- "Basal_vs_LumA"
# subtype_a <- "Her2";   subtype_b <- "LumA"; contrast_name <- "Her2_vs_LumA"
# subtype_a <- "LumA";   subtype_b <- "LumB"; contrast_name <- "LumA_vs_LumB"

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("MSstatsTMT", quietly = TRUE))
    BiocManager::install("MSstatsTMT", update = FALSE, ask = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel", repos = "https://cloud.r-project.org")
})

library(data.table)
library(MSstatsTMT)
library(ggplot2)
library(ggrepel)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# Paths
DATA_DIR <- getwd()
if (!file.exists(file.path(DATA_DIR, "results", "PDC000120", "protein_summary.tsv")))
  DATA_DIR <- file.path(getwd(), "data")
if (!file.exists(file.path(DATA_DIR, "results", "PDC000120", "protein_summary.tsv")))
  stop("Cannot find results/PDC000120/protein_summary.tsv. Run from project root or data/.")

RESULTS_DIR <- file.path(DATA_DIR, "results", "PDC000120")
PROT_SUM_PATH    <- file.path(RESULTS_DIR, "protein_summary.tsv")
SUBTYPE_ANNOT_PATH <- Sys.getenv("PDC_SUBTYPE_ANNOT", unset = "")
if (!nzchar(SUBTYPE_ANNOT_PATH)) {
  SUBTYPE_ANNOT_PATH <- file.path(RESULTS_DIR, "DA_subtype_tumor_only.csv")
} else if (!file.exists(SUBTYPE_ANNOT_PATH)) {
  stop("PDC_SUBTYPE_ANNOT file not found: ", SUBTYPE_ANNOT_PATH)
}
ANNOT_PATH       <- file.path(RESULTS_DIR, "annotation_filled_corrected.csv")
OUT_DIR <- Sys.getenv("PDC_MSSTATSTMT_OUT_DIR", unset = "")
if (!nzchar(OUT_DIR)) {
  OUT_DIR <- RESULTS_DIR
} else {
  OUT_DIR <- normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE)
}
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("Subtype annotation: ", SUBTYPE_ANNOT_PATH)
message("Output directory: ", OUT_DIR)

trim_header <- function(dt) {
  setnames(dt, trimws(gsub("^\uFEFF", "", names(dt))))
  invisible(dt)
}

# =============================================================================
# STEP 1 — Build subtype annotation and filter protein summary
# =============================================================================
message("========== STEP 1: Build subtype annotation for protein summary ==========")

# 1.1 Load subtype design (tumor-only, PAM50)
design <- fread(SUBTYPE_ANNOT_PATH)
trim_header(design)
id_col    <- names(design)[grepl("matrix_sample_id|bioreplicate", names(design), ignore.case = TRUE)][1]
pam50_col <- names(design)[grepl("pam50", names(design), ignore.case = TRUE)][1]
mix_col   <- names(design)[grepl("mixture", names(design), ignore.case = TRUE)][1]
if (is.na(id_col))    id_col    <- "matrix_sample_id"
if (is.na(pam50_col)) pam50_col <- "pam50"
if (is.na(mix_col))   mix_col   <- "mixture"

design[, (pam50_col) := trimws(as.character(get(pam50_col)))]
design <- design[get(pam50_col) != "" & !is.na(get(pam50_col))]

# Pool LumA + LumB to Luminal when Luminal is in the contrast
if ("Luminal" %in% c(subtype_a, subtype_b) || contrast_name == "Luminal_vs_Basal") {
  design[get(pam50_col) %in% c("LumA", "LumB"), (pam50_col) := "Luminal"]
}

# Keep only the two subtypes for this contrast
design <- design[get(pam50_col) %in% c(subtype_a, subtype_b)]
if (id_col != "matrix_sample_id" && id_col %in% names(design)) setnames(design, id_col, "matrix_sample_id")
if (pam50_col != "pam50" && pam50_col %in% names(design)) setnames(design, pam50_col, "pam50")
if (length(mix_col) && !is.na(mix_col) && mix_col != "mixture" && mix_col %in% names(design)) setnames(design, mix_col, "mixture")
if (!"mixture" %in% names(design)) design[, mixture := NA_character_]

# Map BioReplicate (case-insensitive) -> PAM50 and mixture
design[, br_lower := tolower(trimws(matrix_sample_id))]
biorep_to_pam50  <- setNames(design$pam50,  design$br_lower)
biorep_to_mixture <- setNames(design$mixture, design$br_lower)
keep_bioreps <- design$br_lower
message("  Subtype design: ", nrow(design), " samples (", subtype_a, " ", sum(design$pam50 == subtype_a),
        ", ", subtype_b, " ", sum(design$pam50 == subtype_b), ")")

# 1.2 Load protein summary
message("  Loading protein_summary.tsv (this may take a moment for large files)...")
prot <- fread(PROT_SUM_PATH)
trim_header(prot)

req_cols <- c("Protein", "Abundance", "BioReplicate", "Condition", "Mixture", "TechRepMixture", "Run", "Channel")
missing_cols <- setdiff(req_cols, names(prot))
if (length(missing_cols) > 0) stop("protein_summary.tsv missing columns: ", paste(missing_cols, collapse = ", "))

prot[, BioReplicate := trimws(as.character(BioReplicate))]
prot[, Condition_orig := trimws(as.character(Condition))]

# Exclude Norm / reference channels from comparison (do not use them in testing)
prot[, br_lower := tolower(BioReplicate)]
prot <- prot[br_lower %in% keep_bioreps]
prot <- prot[tolower(Condition_orig) != "norm"]
message("  After excluding Norm/reference and keeping only subtype samples: ", nrow(prot), " rows")

# Remap Condition to PAM50 subtype
prot[, Condition := biorep_to_pam50[br_lower]]
prot <- prot[!is.na(Condition)]
prot[, br_lower := NULL][, Condition_orig := NULL]

# Drop NA/blank Abundance
prot <- prot[!is.na(Abundance) & is.finite(as.numeric(Abundance))]
# Coerce Abundance to numeric if needed
prot[, Abundance := as.numeric(Abundance)]
prot <- prot[is.finite(Abundance)]

n_rows_retained <- nrow(prot)
n_uniq_biorep  <- prot[, uniqueN(BioReplicate)]
n_a <- prot[Condition == subtype_a, uniqueN(BioReplicate)]
n_b <- prot[Condition == subtype_b, uniqueN(BioReplicate)]
message("  Protein-level rows retained: ", n_rows_retained)
message("  Unique BioReplicate (samples): ", n_uniq_biorep)
message("  ", subtype_a, " samples: ", n_a)
message("  ", subtype_b, " samples: ", n_b)

# Mixtures per subtype
if ("Mixture" %in% names(prot) && prot[, any(!is.na(Mixture))]) {
  mix_per_subtype <- prot[, .(n_mixtures = uniqueN(Mixture)), by = Condition]
  for (i in seq_len(nrow(mix_per_subtype))) {
    message("  Mixtures in ", mix_per_subtype$Condition[i], ": ", mix_per_subtype$n_mixtures[i])
  }
}

if (n_rows_retained == 0)
  stop("No protein-level rows left after filtering. Check that BioReplicate in protein_summary matches matrix_sample_id in DA_subtype_tumor_only (case-insensitive).")
if (n_a == 0 || n_b == 0)
  stop("One or both subtypes have zero samples. Cannot run contrast.")

# -----------------------------------------------------------------------------
# STEP 1c — Group-wise coverage filter (before DA)
# -----------------------------------------------------------------------------
message("\n========== STEP 1c: Coverage filter (non-missing per group) ==========")
n_tot_samples <- prot[, uniqueN(BioReplicate)]
n_a_samples   <- prot[Condition == subtype_a, uniqueN(BioReplicate)]
n_b_samples   <- prot[Condition == subtype_b, uniqueN(BioReplicate)]
n_a_samples   <- max(1, n_a_samples)
n_b_samples   <- max(1, n_b_samples)

# Per protein: count unique BioReplicate with non-NA, finite Abundance (overall and per Condition)
cov <- prot[, .(
  n_obs_overall = uniqueN(BioReplicate[is.finite(Abundance) & !is.na(Abundance)]),
  n_obs_a       = uniqueN(BioReplicate[Condition == subtype_a & is.finite(Abundance) & !is.na(Abundance)]),
  n_obs_b       = uniqueN(BioReplicate[Condition == subtype_b & is.finite(Abundance) & !is.na(Abundance)])
), by = Protein]
cov[, pct_overall := n_obs_overall / n_tot_samples]
cov[, pct_a       := n_obs_a / n_a_samples]
cov[, pct_b       := n_obs_b / n_b_samples]
keep_proteins <- cov[(pct_overall >= pct_overall_min) & (pct_a >= pct_group_min) & (pct_b >= pct_group_min), Protein]
prot <- prot[Protein %in% keep_proteins]
n_rows_retained <- nrow(prot)
n_proteins_retained <- length(keep_proteins)
message("  Require non-missing in >= ", round(100 * pct_overall_min), "% overall and >= ", round(100 * pct_group_min), "% per subtype.")
message("  Proteins retained after coverage filter: ", n_proteins_retained, " (dropped ", cov[, uniqueN(Protein)] - n_proteins_retained, ")")

# =============================================================================
# STEP 2 — Mixture / subtype balance
# =============================================================================
message("\n========== STEP 2: Mixture / subtype balance ==========")

balance_dt <- NULL
n_mixtures_both <- NA_integer_
if ("Mixture" %in% names(prot) && prot[, uniqueN(Mixture) > 1]) {
  balance_dt <- prot[, .N, by = .(Mixture, Condition)]
  balance_wide <- dcast(balance_dt, Mixture ~ Condition, value.var = "N", fill = 0)
  ca <- balance_wide[[subtype_a]]
  cb <- balance_wide[[subtype_b]]
  if (is.null(ca)) ca <- rep(0, nrow(balance_wide))
  if (is.null(cb)) cb <- rep(0, nrow(balance_wide))
  n_mixtures_both <- sum((as.numeric(ca) > 0) & (as.numeric(cb) > 0))
  message("  Mixtures containing both ", subtype_a, " and ", subtype_b, ": ", n_mixtures_both, " / ", nrow(balance_wide))
  out_balance <- file.path(OUT_DIR, paste0(contrast_name, "_mixture_balance.csv"))
  fwrite(balance_wide, out_balance)
  message("  Wrote ", out_balance)

  if (n_mixtures_both < nrow(balance_wide) / 2) {
    message("  WARNING: Design may be confounded by mixture — many plexes have only one subtype.")
  }
} else {
  message("  Mixture column missing or single level; skipping balance table.")
}

# =============================================================================
# STEP 3 — Run MSstatsTMT groupComparisonTMT
# =============================================================================
message("\n========== STEP 3: Run MSstatsTMT groupComparisonTMT ==========")

# Contrast row: beta_Luminal - beta_Basal. Column order (Basal, Luminal) matches common level ordering.
contrast_matrix <- matrix(c(-1, 1), nrow = 1,
                          dimnames = list(paste0(subtype_a, "-", subtype_b), c("Basal", "Luminal")))

res_dt <- NULL
gc_result <- NULL
tryCatch({
  gc_result <- groupComparisonTMT(
    data = list(ProteinLevelData = as.data.frame(prot)),
    contrast.matrix = contrast_matrix,
    moderated = TRUE,
    adj.method = "BH",
    remove_norm_channel = TRUE,
    remove_empty_channel = TRUE,
    use_log_file = TRUE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = file.path(OUT_DIR, paste0("MSstatsTMT_groupComparison_", contrast_name, ".log"))
  )
  res_dt <- as.data.table(gc_result$ComparisonResult)
  message("  groupComparisonTMT succeeded.")
}, error = function(e) {
  message("  groupComparisonTMT failed: ", conditionMessage(e))
  stop("MSstatsTMT groupComparisonTMT failed. Check mixture balance and that both subtypes appear in multiple runs. ", conditionMessage(e))
})

# Standardize column names
if (!"log2FC" %in% names(res_dt) && "logFC" %in% names(res_dt)) res_dt[, log2FC := logFC]
if (!"adj.pvalue" %in% names(res_dt) && "Adjusted.Pvalue" %in% names(res_dt)) res_dt[, adj.pvalue := Adjusted.Pvalue]
if (!"adj.pvalue" %in% names(res_dt) && "adj.P.Val" %in% names(res_dt)) res_dt[, adj.pvalue := adj.P.Val]
if (!"Protein" %in% names(res_dt) && "ProteinName" %in% names(res_dt)) res_dt[, Protein := ProteinName]
if (!"pvalue" %in% names(res_dt) && "P.Value" %in% names(res_dt)) res_dt[, pvalue := P.Value]
res_dt[, contrast := contrast_name]
res_dt[, method_used := "MSstatsTMT"]

# -----------------------------------------------------------------------------
# Map Protein (RefSeq NP_ / accession) to gene symbol — presentation table
# -----------------------------------------------------------------------------
message("\n========== Gene symbol mapping ==========")
res_dt[, Gene_symbol := NA_character_]
# Named directions match biology: Luminal_up / Basal_up (not subtype_a / subtype_b strings).
res_dt[, direction := fcase(
  log2FC > 0 & adj.pvalue < 0.05, "Luminal_up",
  log2FC < 0 & adj.pvalue < 0.05, "Basal_up",
  default = "NS"
)]
suppressPackageStartupMessages({
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    tryCatch({
      library(org.Hs.eg.db)
      # Strip RefSeq version suffix (e.g. NP_123.2 -> NP_123) for org.Hs.eg.db lookup
      res_dt[, Protein_strip := sub("\\.[0-9]+$", "", as.character(Protein))]
      prots <- unique(res_dt$Protein_strip)
      prots <- prots[grepl("^NP_|^XP_|^NM_", prots)]
      if (length(prots) == 0) prots <- unique(res_dt$Protein_strip)
      sym_map <- suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db,
        keys = prots,
        columns = "SYMBOL",
        keytype = "REFSEQ"
      ))
      if (is.data.frame(sym_map) && nrow(sym_map) > 0 && "REFSEQ" %in% names(sym_map)) {
        sym_map <- as.data.table(sym_map)
        sym_map <- sym_map[!is.na(SYMBOL) & nzchar(SYMBOL)]
        if (nrow(sym_map) > 0) {
          sym_map <- sym_map[, .(Gene_symbol = SYMBOL[1L]), by = REFSEQ]
          setnames(sym_map, "REFSEQ", "Protein_strip")
          res_dt[sym_map, Gene_symbol := i.Gene_symbol, on = "Protein_strip"]
          n_mapped <- res_dt[, sum(!is.na(Gene_symbol) & nzchar(Gene_symbol))]
          message("  Mapped ", n_mapped, " / ", nrow(res_dt), " proteins to gene symbols (org.Hs.eg.db REFSEQ).")
        }
      }
      if ("Protein_strip" %in% names(res_dt)) res_dt[, Protein_strip := NULL]
    }, error = function(e) message("  Gene symbol mapping skipped: ", conditionMessage(e)))
  } else {
    message("  org.Hs.eg.db not installed; run BiocManager::install(\"org.Hs.eg.db\") for symbol mapping.")
  }
})
# Fallback: if Protein is already a gene symbol (e.g. no NP_ prefix), use as-is
res_dt[is.na(Gene_symbol) | !nzchar(Gene_symbol), Gene_symbol := fifelse(
  grepl("^[A-Z0-9]+$", substr(Protein, 1, 2)) & !grepl("^NP_|^NM_|^XP_", as.character(Protein)),
  as.character(Protein), NA_character_
)]

# Presentation table: Protein_accession, Gene_symbol, log2FC, FDR, direction (sorted by FDR)
tab_present <- res_dt[order(adj.pvalue), .(
  Protein_accession = Protein,
  Gene_symbol = fifelse(is.na(Gene_symbol) | !nzchar(Gene_symbol), "", Gene_symbol),
  log2FC = round(log2FC, 4),
  FDR = round(adj.pvalue, 6),
  direction
)]
out_genes <- file.path(OUT_DIR, paste0("DA_MSstatsTMT_", contrast_name, "_gene_symbols.csv"))
fwrite(tab_present, out_genes)
message("  Wrote ", out_genes, " (Protein_accession, Gene_symbol, log2FC, FDR, direction)")

# -----------------------------------------------------------------------------
# Subtype marker sanity check
# -----------------------------------------------------------------------------
markers_basal_up   <- c("KRT5", "KRT14", "KRT17", "EGFR", "FOXC1")
markers_luminal_up <- c("ESR1", "GATA3", "FOXA1", "KRT18", "PGR")
marker_dt <- res_dt[Gene_symbol %in% c(markers_basal_up, markers_luminal_up)]
marker_dt[, expected := fifelse(Gene_symbol %in% markers_basal_up, "Basal_up", "Luminal_up")]
marker_dt[, observed := fcase(log2FC > 0, "Luminal_up", log2FC < 0, "Basal_up", default = "NS")]
marker_dt[, direction_ok := (expected == observed)]
out_markers <- file.path(OUT_DIR, paste0("DA_MSstatsTMT_", contrast_name, "_marker_sanity.csv"))
fwrite(marker_dt[, .(Protein, Gene_symbol, log2FC, adj.pvalue, expected, observed, direction_ok)], out_markers)
message("  Wrote ", out_markers, " (marker sanity check)")

# =============================================================================
# STEP 4 — Outputs
# =============================================================================
message("\n========== STEP 4: Outputs ==========")

out_csv <- file.path(OUT_DIR, paste0("DA_MSstatsTMT_", contrast_name, ".csv"))
res_dt_out <- copy(res_dt)
if ("Gene_symbol" %in% names(res_dt_out)) setcolorder(res_dt_out, c(setdiff(names(res_dt_out), "Gene_symbol"), "Gene_symbol"))
fwrite(res_dt_out, out_csv)
message("  Wrote ", out_csv)

# Summary stats
n_sig <- res_dt[adj.pvalue < 0.05 & abs(log2FC) > 1, .N]
group_a_up <- res_dt[log2FC > 0][order(-log2FC)]
group_b_up <- res_dt[log2FC < 0][order(log2FC)]

out_summary <- file.path(OUT_DIR, paste0("DA_MSstatsTMT_", contrast_name, "_summary.txt"))
sink(out_summary)
cat("Subtype DA (MSstatsTMT) summary — ", contrast_name, "\n", sep = "")
cat("============================================\n\n")
cat("Subtype group sizes (unique BioReplicate):\n")
cat("  ", subtype_a, ": ", n_a, "\n", sep = "")
cat("  ", subtype_b, ": ", n_b, "\n\n", sep = "")
cat("Coverage filter: >= ", round(100 * pct_overall_min), "% non-missing overall, >= ", round(100 * pct_group_min), "% per subtype.\n", sep = "")
cat("Protein-level rows used: ", n_rows_retained, "\n", sep = "")
cat("Proteins tested (after filter): ", nrow(res_dt), "\n", sep = "")
if (!is.null(balance_dt)) {
  cat("Number of mixtures: ", balance_dt[, uniqueN(Mixture)], "\n", sep = "")
  cat("Mixtures containing both subtypes: ", n_mixtures_both, "\n\n", sep = "")
}
cat("Significant (FDR < 0.05, |log2FC| > 1): ", n_sig, "\n\n", sep = "")
cat("Top 10 ", subtype_a, "-up proteins (with gene symbol):\n", sep = "")
top_a <- group_a_up[1:min(10, nrow(group_a_up)), .(Protein, Gene_symbol, log2FC, adj.pvalue)]
if (!"Gene_symbol" %in% names(top_a)) top_a[, Gene_symbol := NA_character_]
print(top_a)
cat("\nTop 10 ", subtype_b, "-up proteins (with gene symbol):\n", sep = "")
top_b <- group_b_up[1:min(10, nrow(group_b_up)), .(Protein, Gene_symbol, log2FC, adj.pvalue)]
if (!"Gene_symbol" %in% names(top_b)) top_b[, Gene_symbol := NA_character_]
print(top_b)
cat("\n--- Subtype marker sanity check ---\n")
cat("Expected ", subtype_a, "-up: ", paste(markers_basal_up, collapse = ", "), "\n", sep = "")
cat("Expected ", subtype_b, "-up: ", paste(markers_luminal_up, collapse = ", "), "\n\n", sep = "")
for (i in seq_len(nrow(marker_dt))) {
  m <- marker_dt[i]
  ok <- if (m$direction_ok) "OK" else "WRONG"
  cat("  ", m$Gene_symbol, " (", m$Protein, "): log2FC = ", round(m$log2FC, 3), ", FDR = ", format(m$adj.pvalue, digits = 3),
      ", expected ", m$expected, ", observed ", m$observed, " [", ok, "]\n", sep = "")
}
n_ok <- marker_dt[, sum(direction_ok)]
cat("\n  Direction match: ", n_ok, " / ", nrow(marker_dt), " markers.\n", sep = "")
if (!is.null(balance_dt) && !is.na(n_mixtures_both) && balance_dt[, uniqueN(Mixture)] > 1 && n_mixtures_both < balance_dt[, uniqueN(Mixture)] / 2)
  cat("\nWARNING: Contrast may be confounded by mixture (many plexes have only one subtype).\n")
cat("\nInterpretation: positive log2FC = higher abundance in Luminal than Basal (Luminal - Basal).\n")
sink(NULL)
message("  Wrote ", out_summary)

# =============================================================================
# STEP 5 — Volcano plot
# =============================================================================
message("\n========== STEP 5: Volcano plot ==========")

# Use Protein as label (often RefSeq/UniProt; optional gene symbol mapping can be added)
plot_dt <- copy(res_dt)
plot_dt[, neglog10FDR := -log10(pmax(adj.pvalue, 1e-300))]
plot_dt[, sig := abs(log2FC) > 1 & adj.pvalue < 0.05]
top_n <- 20L
plot_dt[, rank := rank(adj.pvalue, ties.method = "first")]
plot_dt[, label := fifelse(sig & rank <= top_n, as.character(Protein), NA_character_)]

# Use Gene_symbol for labels when available (presentation-friendly)
plot_dt[!is.na(Gene_symbol) & nzchar(Gene_symbol), label := fifelse(sig & rank <= top_n, Gene_symbol, label)]
# Highlight subtype markers
all_markers <- c(markers_basal_up, markers_luminal_up)
plot_dt[Gene_symbol %in% all_markers, label := Gene_symbol]

out_pdf <- file.path(OUT_DIR, paste0("volcano_MSstatsTMT_", contrast_name, ".pdf"))
p <- ggplot(plot_dt, aes(x = log2FC, y = neglog10FDR)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.1) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55") +
  ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3, box.padding = 0.4, min.segment.length = 0) +
  labs(
    title = paste0("MSstatsTMT ", contrast_name, " (PDC000120)"),
    x = "log2FC (Luminal - Basal)",
    y = "-log10(FDR)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(out_pdf, p, width = 8.5, height = 6.5)
message("  Wrote ", out_pdf)

# =============================================================================
# Done
# =============================================================================
message("\n========== Done ==========")
message("Contrast: ", contrast_name, " (positive log2FC = ", subtype_a, " up)")
message("Significant (FDR<0.05, |log2FC|>1): ", n_sig)
message("Outputs: ", out_csv, ", ", out_pdf, ", ", out_summary)
