#!/usr/bin/env Rscript
# CCLE Luminal vs Basal using MSstats (not MSstatsTMT) on the summarized gene matrix.
# Mirrors data/scripts/DA_subtype_MSstats_PDC000120.R: long format -> dataProcess -> groupComparison.
#
# Contrast: Luminal_vs_Basal with positive log2FC = higher in Luminal lines.
#
# Run from repo root:
#   Rscript data/scripts/ccle_DA_luminal_basal_MSstats.R

suppressPackageStartupMessages({
  if (!requireNamespace("MSstats", quietly = TRUE))
    BiocManager::install("MSstats", update = FALSE, ask = FALSE)
  if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", update = FALSE, ask = FALSE)
  library(MSstats)
  library(limma)
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Paths
# Optional args: [gene_matrix.csv] [out_dir]
# Default gene matrix: data/results/CCLE_corrected/gene_matrix.csv (falls back to CCLE/).
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
root <- getwd()
find_gm <- function() {
  candidates <- c(
    file.path(root, "data", "results", "CCLE_corrected", "gene_matrix.csv"),
    file.path(root, "data", "results", "CCLE", "gene_matrix.csv"),
    file.path(normalizePath(file.path(getwd(), "..")), "data", "results", "CCLE_corrected", "gene_matrix.csv"),
    file.path(normalizePath(file.path(getwd(), "..")), "data", "results", "CCLE", "gene_matrix.csv"),
    file.path(normalizePath(file.path(getwd(), "..", "..")), "data", "results", "CCLE_corrected", "gene_matrix.csv"),
    file.path(normalizePath(file.path(getwd(), "..", "..")), "data", "results", "CCLE", "gene_matrix.csv")
  )
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  NULL
}
if (length(args) >= 1L) {
  gm_path <- args[1L]
  if (!file.exists(gm_path)) stop("Gene matrix not found: ", gm_path)
  gm_path <- normalizePath(gm_path)
} else {
  gm_path <- find_gm()
  if (is.null(gm_path)) stop("Cannot find data/results/CCLE_corrected/gene_matrix.csv (or CCLE fallback); pass path as arg1.")
}
out_dir <- if (length(args) >= 2L) args[2L] else file.path(dirname(gm_path), "DA_luminal_vs_basal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Gene matrix: ", gm_path)
message("Output dir: ", normalizePath(out_dir, mustWork = FALSE))

luminal_lines <- c("MCF7", "T-47D", "CAMA-1", "ZR-75-1")
basal_lines   <- c("HCC 1806", "HCC1143", "HCC70", "MDA-MB-468")
samples       <- c(luminal_lines, basal_lines)

group1 <- "Luminal"  # numerator in contrast: positive log2FC = Luminal higher
group2 <- "Basal"
contrast_name <- "Luminal_vs_Basal"

pct_overall_min <- 0.35
pct_group_min   <- 0.25

trim_header <- function(dt) {
  setnames(dt, trimws(gsub("^\uFEFF", "", names(dt))))
  invisible(dt)
}

message("========== CCLE MSstats (matrix-derived) ==========")
message("Using matrix: ", gm_path)

gm <- fread(gm_path, showProgress = FALSE)
trim_header(gm)
miss <- setdiff(samples, names(gm))
if (length(miss)) stop("Missing columns in gene matrix: ", paste(miss, collapse = ", "))

# Use UniProt column as protein ID (matches pipeline convention for CCLE)
id_col <- if ("UniProtID" %in% names(gm)) "UniProtID" else names(gm)[2]
gm_mat <- as.matrix(gm[, ..samples])
storage.mode(gm_mat) <- "double"
rownames(gm_mat) <- gm[[id_col]]

n_tot <- ncol(gm_mat)
n_g1 <- sum(samples %in% luminal_lines)
n_g2 <- sum(samples %in% basal_lines)
is_grp1 <- samples %in% luminal_lines
is_grp2 <- samples %in% basal_lines

keep_genes <- vapply(seq_len(nrow(gm_mat)), function(i) {
  x <- as.numeric(gm_mat[i, ])
  valid <- !is.na(x) & is.finite(x)
  pct_overall <- sum(valid) / n_tot
  pct_grp1   <- sum(valid[is_grp1]) / max(1L, n_g1)
  pct_grp2   <- sum(valid[is_grp2]) / max(1L, n_g2)
  (pct_overall >= pct_overall_min) && (pct_grp1 >= pct_group_min) && (pct_grp2 >= pct_group_min)
}, logical(1))

gm_mat <- gm_mat[keep_genes, , drop = FALSE]
# MSstats requires unique ProteinName keys
rn <- rownames(gm_mat)
if (any(duplicated(rn))) {
  message("Deduplicating ", sum(duplicated(rn)), " duplicate protein row names (keep first).")
  gm_mat <- gm_mat[!duplicated(rn), , drop = FALSE]
}
n_genes <- nrow(gm_mat)
message("Genes after coverage filter: ", n_genes)

design_dt <- data.table(
  matrix_sample_id = samples,
  pam50 = ifelse(samples %in% luminal_lines, "Luminal", "Basal"),
  matrix_col = samples
)
design_dt[, Run := matrix_col]
design_dt[, BioReplicate := matrix_col]
design_dt[, Condition := pam50]
cond_per_col <- setNames(design_dt$Condition, design_dt$matrix_col)

proteins <- rownames(gm_mat)
long_list <- lapply(seq_len(n_genes), function(i) {
  log2val <- as.numeric(gm_mat[i, ])
  # Unique synthetic peptide per protein — shared "PROTEIN_SUM" can break MSstats joins
  pep <- paste0("SYN_", i)
  data.table(
    ProteinName      = proteins[i],
    PeptideSequence  = pep,
    PrecursorCharge  = 0L,
    FragmentIon      = "NA",
    ProductCharge    = 0L,
    IsotopeLabelType = "L",
    Condition        = cond_per_col[colnames(gm_mat)],
    BioReplicate     = colnames(gm_mat),
    Run              = colnames(gm_mat),
    Intensity        = 2^log2val
  )
})
msstats_long <- rbindlist(long_list)
msstats_long <- msstats_long[!is.na(Intensity) & is.finite(Intensity)]
# Collapse accidental duplicate (ProteinName, Run) keys
msstats_long <- msstats_long[, .(Intensity = mean(Intensity)), by = .(
  ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge,
  IsotopeLabelType, Condition, BioReplicate, Run
)]
message("Long-format rows: ", nrow(msstats_long))

method_used <- "limma"
res_dt <- NULL
msstats_ok <- FALSE

tryCatch({
  processed <- dataProcess(
    as.data.frame(msstats_long),
    normalization = FALSE,
    summaryMethod = "TMP",
    MBimpute = FALSE,
    use_log_file = FALSE
  )
  # Contrast columns must match the *level order* MSstats uses (usually alphabetical: Basal, Luminal)
  lev <- sort(unique(as.character(msstats_long$Condition)))
  contrast_mat <- matrix(0, nrow = 1, ncol = length(lev))
  colnames(contrast_mat) <- lev
  rownames(contrast_mat) <- contrast_name
  # Luminal - Basal: for lev == c("Basal","Luminal") use (-1, 1)
  contrast_mat[1, ] <- vapply(lev, function(x) {
    if (x == "Luminal") 1 else if (x == "Basal") -1 else NA_real_
  }, numeric(1))
  comparison <- groupComparison(contrast.matrix = contrast_mat, data = processed)
  res_ms <- as.data.frame(comparison$ComparisonResult)
  res_ms$ProteinName <- as.character(res_ms$ProteinName)
  if (!"log2FC" %in% names(res_ms) && "logFC" %in% names(res_ms)) res_ms$log2FC <- res_ms$logFC
  if (!"adj.pvalue" %in% names(res_ms) && "Adjusted.Pvalue" %in% names(res_ms)) res_ms$adj.pvalue <- res_ms$Adjusted.Pvalue
  if (!"pvalue" %in% names(res_ms) && "P.Value" %in% names(res_ms)) res_ms$pvalue <- res_ms$P.Value
  res_dt <- as.data.table(res_ms)
  method_used <- "MSstats"
  msstats_ok <- TRUE
  message("MSstats groupComparison succeeded.")
}, error = function(e) {
  message("MSstats failed: ", conditionMessage(e))
  message("Falling back to limma on the same filtered matrix.")
})

if (!msstats_ok) {
  group <- factor(rep(c("Luminal", "Basal"), c(4L, 4L)), levels = c("Luminal", "Basal"))
  design_limma <- model.matrix(~ 0 + group)
  colnames(design_limma) <- c("Luminal", "Basal")
  fit <- lmFit(gm_mat, design_limma)
  ctr <- limma::makeContrasts(Luminal_vs_Basal = Luminal - Basal, levels = design_limma)
  fit2 <- contrasts.fit(fit, ctr)
  fit2 <- eBayes(fit2)
  res_limma <- topTable(fit2, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
  res_dt <- data.table(
    ProteinName = rownames(res_limma),
    log2FC = res_limma$logFC,
    pvalue = res_limma$P.Value,
    adj.pvalue = res_limma$adj.P.Val
  )
  method_used <- "limma_fallback"
}

res_dt[, contrast := contrast_name]
res_dt[, method_used := method_used]
if ("ProteinName" %in% names(res_dt)) setnames(res_dt, "ProteinName", "ProteinID")

extract_symbol <- function(p) {
  p <- as.character(p)
  m <- regexec("^[^|]*\\|[^|]*\\|([A-Z0-9]+)_HUMAN", p)
  r <- regmatches(p, m)
  if (length(r[[1]]) >= 2) return(r[[1]][2])
  NA_character_
}
res_dt[, Gene_symbol := vapply(ProteinID, extract_symbol, character(1))]

out_csv <- file.path(out_dir, paste0("DA_MSstats_", contrast_name, ".csv"))
fwrite(res_dt, out_csv)
message("Wrote ", out_csv)

# Directionality (literature markers)
markers_luminal <- c("ESR1", "GATA3", "FOXA1", "KRT18", "PGR")
markers_basal   <- c("KRT5", "KRT14", "KRT17", "EGFR", "FOXC1")
mk <- res_dt[Gene_symbol %in% c(markers_luminal, markers_basal)]
mk[, expected := fifelse(Gene_symbol %in% markers_luminal, "higher_in_Luminal", "higher_in_Basal")]
mk[, observed_sign := ifelse(log2FC > 0, "Luminal_higher", ifelse(log2FC < 0, "Basal_higher", "zero"))]
mk[, direction_ok := fifelse(Gene_symbol %in% markers_luminal, log2FC > 0, log2FC < 0)]
fwrite(
  mk[, .(Gene_symbol, ProteinID, log2FC, pvalue, adj.pvalue, expected, observed_sign, direction_ok)],
  file.path(out_dir, "directionality_markers_MSstats.csv")
)

writeLines(c(
  "Directionality — MSstats (matrix-derived pipeline)",
  "==================================================",
  "Contrast name: Luminal_vs_Basal",
  "  log2FC = estimate for Luminal minus Basal (same convention as PDC MSstats script).",
  "  Positive log2FC => higher abundance in Luminal lines on average.",
  "",
  "Implementation: summarized gene matrix converted to long format (Intensity = 2^log2),",
  "dataProcess(normalization=FALSE), groupComparison. If MSstats errors, limma on the same matrix is used.",
  "",
  "See directionality_markers_MSstats.csv for literature markers vs sign of log2FC."
), file.path(out_dir, "directionality_key_MSstats.txt"))

# Volcano (ggplot)
pv <- if ("pvalue" %in% names(res_dt)) res_dt$pvalue else res_dt[[grep("^P", names(res_dt), value = TRUE)[1]]]
res_dt[, p_plot := pmax(pv, 1e-300, na.rm = TRUE)]
p_volc <- ggplot(res_dt, aes(log2FC, -log10(p_plot))) +
  geom_point(aes(colour = adj.pvalue < 0.05), alpha = 0.35, size = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey55")) +
  theme_bw() +
  labs(
    title = paste0("CCLE Luminal vs Basal (", method_used, ")"),
    subtitle = "Positive log2FC = Luminal > Basal; matrix-derived MSstats pipeline",
    x = "log2FC (Luminal - Basal)",
    y = "-log10 p",
    colour = "adj.P < 0.05"
  )
ggsave(file.path(out_dir, paste0("volcano_MSstats_", contrast_name, ".pdf")), p_volc, width = 7.5, height = 5)

sink(file.path(out_dir, paste0("DA_MSstats_", contrast_name, "_summary.txt")))
cat("CCLE — MSstats-style DA\n")
cat("Method: ", method_used, "\n", sep = "")
cat("Contrast: positive log2FC = ", group1, " > ", group2, "\n", sep = "")
cat("Genes tested: ", nrow(res_dt), "\n", sep = "")
if (!msstats_ok) {
  cat("\nNote: MSstats::groupComparison did not return usable results (see console); limma fallback on the same filtered matrix.\n")
  cat("Contrast unchanged: Luminal - Basal.\n")
}
sink()

message("Done. Outputs in ", out_dir)
