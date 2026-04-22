#!/usr/bin/env Rscript
# =============================================================================
# CCLE breast: Luminal vs Basal — MSstatsTMT::groupComparisonTMT on protein_summary.tsv
# (mirrors data/scripts/DA_subtype_MSstatsTMT_PDC000120.R).
#
# Design: exactly 8 cell lines — 4 Luminal (MCF7, T-47D, CAMA-1, ZR-75-1) vs
#         4 Basal (HCC 1806, HCC1143, HCC70, MDA-MB-468). One BioReplicate per line.
#
# Optional gene-matrix PCA (visualization only): loads gene_matrix.csv in same folder.
#
# Usage (repo root):
#   Rscript data/scripts/ccle_DA_luminal_basal_v1.R
# Optional: [protein_summary.tsv] [out_dir]
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("MSstatsTMT", quietly = TRUE))
    BiocManager::install("MSstatsTMT", update = FALSE, ask = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel", repos = "https://cloud.r-project.org")
  library(data.table)
  library(MSstatsTMT)
  library(ggplot2)
  library(ggrepel)
})

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# --- Paths ---
find_data_root <- function() {
  wd <- getwd()
  if (file.exists(file.path(wd, "data", "results", "CCLE_corrected", "protein_summary.tsv")))
    return(normalizePath(file.path(wd, "data")))
  if (file.exists(file.path(wd, "results", "CCLE_corrected", "protein_summary.tsv")))
    return(normalizePath(wd))
  stop("Cannot find data/results/CCLE_corrected/protein_summary.tsv (run from repo root or data/).")
}
DATA_DIR <- find_data_root()
CCLE_RES <- file.path(DATA_DIR, "results", "CCLE_corrected")
PROT_SUM_PATH <- file.path(CCLE_RES, "protein_summary.tsv")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && file.exists(args[1L])) PROT_SUM_PATH <- normalizePath(args[1L])
out_dir <- if (length(args) >= 2L) args[2L] else file.path(CCLE_RES, "DA_luminal_vs_basal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
OUT_DIR <- normalizePath(out_dir, mustWork = FALSE)

# --- Contrast (same estimand as CPTAC) ---
subtype_a <- "Luminal"
subtype_b <- "Basal"
contrast_name <- "Luminal_vs_Basal"
pct_overall_min <- 0.35
pct_group_min <- 0.25

luminal_lines <- c("MCF7", "T-47D", "CAMA-1", "ZR-75-1")
basal_lines <- c("HCC 1806", "HCC1143", "HCC70", "MDA-MB-468")
keep_lines <- c(luminal_lines, basal_lines)

message("protein_summary: ", PROT_SUM_PATH)
message("Output: ", OUT_DIR)
message("Subtype lines: 4 Luminal + 4 Basal (8 BioReplicates total).")

trim_header <- function(dt) {
  setnames(dt, trimws(gsub("^\uFEFF", "", names(dt))))
  invisible(dt)
}

# =============================================================================
# Load protein summary and restrict to 8 lines
# =============================================================================
message("Loading protein_summary.tsv (large file; may take a few minutes)...")
prot <- fread(PROT_SUM_PATH, showProgress = TRUE)
trim_header(prot)

req_cols <- c("Protein", "Abundance", "BioReplicate", "Condition", "Mixture", "TechRepMixture", "Run", "Channel")
missing_cols <- setdiff(req_cols, names(prot))
if (length(missing_cols) > 0) stop("protein_summary.tsv missing columns: ", paste(missing_cols, collapse = ", "))

prot[, BioReplicate := trimws(as.character(BioReplicate))]
prot[, br_lower := tolower(BioReplicate)]
# Map BioReplicate (case-insensitive) to Luminal / Basal
line_to_subtype <- c(
  setNames(rep("Luminal", 4L), tolower(luminal_lines)),
  setNames(rep("Basal", 4L), tolower(basal_lines))
)
prot <- prot[br_lower %in% names(line_to_subtype)]
prot[, Condition := line_to_subtype[br_lower]]
# Drop bridge pool rows if present
prot <- prot[tolower(trimws(as.character(BioReplicate))) != "pool"]

prot <- prot[!is.na(Abundance) & is.finite(as.numeric(Abundance))]
prot[, Abundance := as.numeric(Abundance)]
prot[, br_lower := NULL]

n_a <- prot[Condition == subtype_a, uniqueN(BioReplicate)]
n_b <- prot[Condition == subtype_b, uniqueN(BioReplicate)]
message("  Luminal lines (unique BioReplicate): ", n_a, " ; Basal: ", n_b)
if (n_a != 4L || n_b != 4L) {
  message("  WARNING: Expected exactly 4 Luminal + 4 Basal lines. Found: ", subtype_a, "=", n_a, ", ", subtype_b, "=", n_b)
}

# =============================================================================
# Coverage filter (same rule as PDC000120 MSstatsTMT script)
# =============================================================================
message("Coverage filter (>= ", 100 * pct_overall_min, "% overall, >= ", 100 * pct_group_min, "% per subtype)...")
n_tot_samples <- prot[, uniqueN(BioReplicate)]
n_a_samples <- max(1L, prot[Condition == subtype_a, uniqueN(BioReplicate)])
n_b_samples <- max(1L, prot[Condition == subtype_b, uniqueN(BioReplicate)])

cov <- prot[, .(
  n_obs_overall = uniqueN(BioReplicate[is.finite(Abundance) & !is.na(Abundance)]),
  n_obs_a = uniqueN(BioReplicate[Condition == subtype_a & is.finite(Abundance) & !is.na(Abundance)]),
  n_obs_b = uniqueN(BioReplicate[Condition == subtype_b & is.finite(Abundance) & !is.na(Abundance)])
), by = Protein]
cov[, pct_overall := n_obs_overall / n_tot_samples]
cov[, pct_a := n_obs_a / n_a_samples]
cov[, pct_b := n_obs_b / n_b_samples]
keep_proteins <- cov[(pct_overall >= pct_overall_min) & (pct_a >= pct_group_min) & (pct_b >= pct_group_min), Protein]
prot <- prot[Protein %in% keep_proteins]
message("  Proteins retained: ", length(keep_proteins))

if (nrow(prot) == 0L) stop("No protein rows left after filtering.")

# Mixture balance table (informational)
if ("Mixture" %in% names(prot) && prot[, uniqueN(Mixture) > 1]) {
  balance_dt <- prot[, .N, by = .(Mixture, Condition)]
  balance_wide <- dcast(balance_dt, Mixture ~ Condition, value.var = "N", fill = 0)
  fwrite(balance_wide, file.path(OUT_DIR, paste0(contrast_name, "_mixture_balance.csv")))
}

# =============================================================================
# MSstatsTMT groupComparisonTMT (same arguments as CPTAC)
# =============================================================================
message("Running groupComparisonTMT ...")
contrast_matrix <- matrix(c(-1, 1), nrow = 1,
                          dimnames = list(paste0(subtype_a, "-", subtype_b), c("Basal", "Luminal")))

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
if (!"log2FC" %in% names(res_dt) && "logFC" %in% names(res_dt)) res_dt[, log2FC := logFC]
if (!"adj.pvalue" %in% names(res_dt) && "Adjusted.Pvalue" %in% names(res_dt)) res_dt[, adj.pvalue := Adjusted.Pvalue]
if (!"adj.pvalue" %in% names(res_dt) && "adj.P.Val" %in% names(res_dt)) res_dt[, adj.pvalue := adj.P.Val]
if (!"Protein" %in% names(res_dt) && "ProteinName" %in% names(res_dt)) res_dt[, Protein := ProteinName]
if (!"pvalue" %in% names(res_dt) && "P.Value" %in% names(res_dt)) res_dt[, pvalue := P.Value]
res_dt[, contrast := contrast_name]
res_dt[, method_used := "MSstatsTMT"]

# --- Gene symbols: UniProt from sp|ACC| or tr|ACC| ---
res_dt[, Gene_symbol := NA_character_]
extract_acc <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexec("(?:^#*)?(?:sp|tr)\\|([^|]+)", x))
  vapply(m, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
}
res_dt[, uniprot_acc := extract_acc(Protein)]
res_dt[, uniprot_acc := sub("-[0-9]+$", "", uniprot_acc)]

suppressPackageStartupMessages({
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    tryCatch({
      library(org.Hs.eg.db)
      accs <- unique(res_dt[!is.na(uniprot_acc), uniprot_acc])
      sym_map <- suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db,
        keys = accs,
        columns = "SYMBOL",
        keytype = "UNIPROT"
      ))
      if (is.data.frame(sym_map) && nrow(sym_map) > 0) {
        sym_map <- as.data.table(sym_map)
        sym_map <- sym_map[!is.na(SYMBOL) & nzchar(SYMBOL)]
        sym_map <- sym_map[, .(Gene_symbol = SYMBOL[1L]), by = UNIPROT]
        setnames(sym_map, "UNIPROT", "uniprot_acc")
        res_dt[sym_map, Gene_symbol := i.Gene_symbol, on = "uniprot_acc"]
      }
      message("  Mapped ", res_dt[!is.na(Gene_symbol) & nzchar(Gene_symbol), .N], " / ", nrow(res_dt), " proteins to gene symbols.")
    }, error = function(e) message("  Gene symbol mapping skipped: ", conditionMessage(e)))
  } else {
    message("  Install org.Hs.eg.db for gene symbol mapping (optional).")
  }
})
res_dt[, uniprot_acc := NULL]

markers_basal_up <- c("KRT5", "KRT14", "KRT17", "EGFR", "FOXC1")
markers_luminal_up <- c("ESR1", "GATA3", "FOXA1", "KRT18", "PGR")
marker_dt <- res_dt[Gene_symbol %in% c(markers_basal_up, markers_luminal_up)]
marker_dt[, expected := fifelse(Gene_symbol %in% markers_basal_up, "Basal_up", "Luminal_up")]
marker_dt[, observed := fcase(log2FC > 0, "Luminal_up", log2FC < 0, "Basal_up", default = "NS")]
marker_dt[, direction_ok := (expected == observed)]

# --- Write outputs (names aligned with CPTAC) ---
out_csv <- file.path(OUT_DIR, paste0("DA_MSstatsTMT_", contrast_name, ".csv"))
fwrite(res_dt, out_csv)
message("  Wrote ", out_csv)

fwrite(res_dt, file.path(OUT_DIR, "DA_luminal_vs_basal_MSstatsTMT_protein.csv"))

# Bridge file for subtype_cptac_ccle_benchmark.R / build_subtype_marker_panel.R (expects limma-like columns)
bridge_limma <- res_dt[, .(
  UniProtID = Protein,
  logFC = log2FC,
  adj.P.Val = adj.pvalue,
  P.Value = pvalue,
  Gene_symbol = Gene_symbol
)]
fwrite(bridge_limma, file.path(OUT_DIR, "DA_luminal_vs_basal_limma.csv"))
message("  Wrote DA_luminal_vs_basal_limma.csv (MSstatsTMT results; limma-compatible columns for downstream scripts)")

tab_present <- res_dt[order(adj.pvalue), .(
  Protein_accession = Protein,
  Gene_symbol = fifelse(is.na(Gene_symbol) | !nzchar(Gene_symbol), "", Gene_symbol),
  log2FC = round(log2FC, 4),
  FDR = round(adj.pvalue, 6),
  direction = fcase(
    log2FC > 0 & adj.pvalue < 0.05, "Luminal_up",
    log2FC < 0 & adj.pvalue < 0.05, "Basal_up",
    default = "NS"
)
)]
fwrite(tab_present, file.path(OUT_DIR, paste0("DA_MSstatsTMT_", contrast_name, "_gene_symbols.csv")))

fwrite(marker_dt[, .(Protein, Gene_symbol, log2FC, adj.pvalue, expected, observed, direction_ok)],
  file.path(OUT_DIR, paste0("DA_MSstatsTMT_", contrast_name, "_marker_sanity.csv"))
)

# Summary
n_sig <- res_dt[adj.pvalue < 0.05 & abs(log2FC) > 1, .N]
out_summary <- file.path(OUT_DIR, paste0("DA_MSstatsTMT_", contrast_name, "_summary.txt"))
sink(out_summary)
cat("CCLE subtype DA (MSstatsTMT) — ", contrast_name, "\n", sep = "")
cat("============================================\n\n")
cat("Design: 4 Luminal cell lines vs 4 Basal (8 BioReplicates).\n")
cat("  Luminal: ", paste(luminal_lines, collapse = ", "), "\n", sep = "")
cat("  Basal:   ", paste(basal_lines, collapse = ", "), "\n\n", sep = "")
cat("Method: MSstatsTMT::groupComparisonTMT on protein_summary.tsv (same settings as PDC000120).\n")
cat("Estimand: positive log2FC => higher in Luminal (Luminal - Basal).\n\n")
cat("Proteins tested: ", nrow(res_dt), "\n", sep = "")
cat("Significant (FDR < 0.05, |log2FC| > 1): ", n_sig, "\n\n", sep = "")
cat("Marker direction match: ", marker_dt[, sum(direction_ok)], " / ", nrow(marker_dt), "\n", sep = "")
sink(NULL)
message("  Wrote ", out_summary)

# Volcano
plot_dt <- copy(res_dt)
plot_dt[, neglog10FDR := -log10(pmax(adj.pvalue, 1e-300))]
plot_dt[, sig := abs(log2FC) > 1 & adj.pvalue < 0.05]
plot_dt[, rank := rank(adj.pvalue, ties.method = "first")]
plot_dt[, label := fifelse(sig & rank <= 20L, as.character(Protein), NA_character_)]
plot_dt[!is.na(Gene_symbol) & nzchar(Gene_symbol), label := fifelse(sig & rank <= 20L, Gene_symbol, label)]
all_markers <- c(markers_basal_up, markers_luminal_up)
plot_dt[Gene_symbol %in% all_markers, label := Gene_symbol]

out_pdf <- file.path(OUT_DIR, paste0("volcano_MSstatsTMT_", contrast_name, ".pdf"))
p <- ggplot(plot_dt, aes(x = log2FC, y = neglog10FDR)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.1) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55") +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3, box.padding = 0.4, min.segment.length = 0) +
  labs(
    title = paste0("MSstatsTMT ", contrast_name, " (CCLE corrected)"),
    subtitle = "4 Luminal vs 4 Basal lines; exploratory — not powered",
    x = "log2FC (Luminal - Basal)",
    y = "-log10(FDR)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(out_pdf, p, width = 8.5, height = 6.5, device = grDevices::pdf)
message("  Wrote ", out_pdf)

# =============================================================================
# Gene-matrix PCA (visualization only; same 8 columns as before)
# =============================================================================
gm_path <- file.path(CCLE_RES, "gene_matrix.csv")
if (file.exists(gm_path)) {
  message("PCA on gene_matrix (optional visualization)...")
  gm <- fread(gm_path, showProgress = FALSE)
  miss <- setdiff(keep_lines, names(gm))
  if (length(miss) == 0) {
    uid <- gm[[2]]
    M <- as.matrix(gm[, ..keep_lines])
    storage.mode(M) <- "numeric"
    M[!is.finite(M)] <- NA
    ok <- rowSums(is.finite(M[, luminal_lines])) >= 2L & rowSums(is.finite(M[, basal_lines])) >= 2L
    M <- M[ok, , drop = FALSE]
    Mp <- M
    for (j in seq_len(ncol(Mp))) {
      v <- Mp[, j]
      med <- median(v[is.finite(v)], na.rm = TRUE)
      v[!is.finite(v)] <- med
      Mp[, j] <- v
    }
    pc <- prcomp(t(Mp), center = TRUE, scale. = TRUE)
    group <- factor(rep(c("Luminal", "Basal"), c(4L, 4L)), levels = c("Luminal", "Basal"))
    pcs <- data.frame(
      sample = keep_lines,
      PC1 = pc$x[, 1],
      PC2 = pc$x[, 2],
      group = group
    )
    fwrite(pcs, file.path(OUT_DIR, "ccle_pca_scores.csv"))
    p_pca <- ggplot(pcs, aes(PC1, PC2, colour = group, label = sample)) +
      geom_point(size = 3, alpha = 0.85) +
      geom_text_repel(size = 2.5, max.overlaps = 20) +
      theme_bw() +
      labs(
        title = "PCA — CCLE gene matrix (8 lines, filtered genes)",
        subtitle = "MSstatsTMT DE is protein-level; PCA is gene-level for visualization"
      )
    ggsave(file.path(OUT_DIR, "ccle_luminal_basal_pca.pdf"), p_pca, width = 7, height = 5, device = grDevices::pdf)
    message("  Wrote ccle_luminal_basal_pca.pdf")
  }
}

# Canonical markers vs MSstatsTMT
marker_patterns <- c(
  ESR1 = "P03372", GATA3 = "P23771", FOXA1 = "P55317", KRT18 = "P05783", PGR = "P06454",
  KRT5 = "P13647", KRT14 = "P02533", KRT17 = "Q04695", EGFR = "P00533", FOXC1 = "Q12948"
)
mk_rows <- lapply(names(marker_patterns), function(nm) {
  pat <- marker_patterns[[nm]]
  hit <- res_dt[(Gene_symbol == nm) | grepl(pat, Protein, fixed = TRUE)]
  if (nrow(hit) == 0L) {
    return(data.table(gene = nm, pattern = pat, in_MSstatsTMT = FALSE))
  }
  h <- hit[order(adj.pvalue)][1L]
  data.table(
    gene = nm, pattern = pat, in_MSstatsTMT = TRUE,
    Protein = h$Protein, log2FC = h$log2FC, adj.pvalue = h$adj.pvalue
  )
})
mk <- rbindlist(mk_rows, fill = TRUE)
fwrite(mk, file.path(OUT_DIR, "canonical_markers_check.csv"))

writeLines(c(
  "CCLE Luminal vs Basal — MSstatsTMT (primary)",
  "=============================================",
  "Lines: 4 Luminal (MCF7, T-47D, CAMA-1, ZR-75-1) vs 4 Basal (HCC 1806, HCC1143, HCC70, MDA-MB-468).",
  "Inference: MSstatsTMT::groupComparisonTMT on protein_summary.tsv (same settings as CPTAC PDC000120).",
  "Gene-level limma (optional legacy): data/scripts/ccle_DA_luminal_basal_limma_gene_matrix.R",
  "",
  "Outputs: DA_MSstatsTMT_Luminal_vs_Basal*.csv, volcano_MSstatsTMT_*.pdf, PCA (if gene_matrix present)."
), file.path(OUT_DIR, "README.txt"))

message("Done. Outputs in ", OUT_DIR)
