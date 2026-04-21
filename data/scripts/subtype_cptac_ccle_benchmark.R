#!/usr/bin/env Rscript
# CPTAC vs CCLE subtype benchmark for slides (primary CCLE = Table S2 / ccle_sum — denser overlap vs CPTAC).
# - subtype_summary_table.csv (2 rows: CPTAC, CCLE Table S2)
# - subtype_benchmark_metadata.txt (thresholds, column names, sign convention, sample counts)
# - venn_luminal_up_cptac_vs_ccle.png + venn_basal_up_cptac_vs_ccle.png (diagram-only; wide canvas)
# - optional venn_*_with_gene_lists.png (diagram + gene text)
# - optional venn_*_cptac_vs_ccle_corrected.png if CCLE_corrected limma exists
# - subtype_marker_sanity.csv (markers: direction first, significance second)
# - subtype_slide_highlights.txt
#
# Run from repository root:
#   Rscript --no-init-file data/scripts/subtype_cptac_ccle_benchmark.R

suppressPackageStartupMessages({
  if (!requireNamespace("ggVennDiagram", quietly = TRUE))
    install.packages("ggVennDiagram", repos = "https://cloud.r-project.org", quiet = TRUE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", repos = "https://cloud.r-project.org", quiet = TRUE)
  if (!requireNamespace("patchwork", quietly = TRUE))
    install.packages("patchwork", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(ggVennDiagram)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

repo <- getwd()
if (!file.exists(file.path(repo, "data", "results", "PDC000120"))) {
  repo <- normalizePath(file.path(getwd(), ".."))
}
if (!file.exists(file.path(repo, "data", "results", "PDC000120"))) {
  repo <- normalizePath(file.path(getwd(), "..", ".."))
}
stopifnot(file.exists(file.path(repo, "reports")))

FDR_THR <- 0.05
MAX_GENES_SHOW_CPTAC_ONLY <- 45L
# Wider figures give the Venn circles more horizontal room (inches at 300 dpi).
VENN_DIAGRAM_WIDTH <- 13
VENN_DIAGRAM_HEIGHT <- 7.5
VENN_WITH_LISTS_WIDTH <- 16

path_cptac <- file.path(repo, "data", "results", "PDC000120", "DA_subtype_subset_runs", "DA_MSstatsTMT_Luminal_vs_Basal.csv")
path_cptac_sum <- file.path(repo, "data", "results", "PDC000120", "DA_subtype_subset_runs", "DA_MSstatsTMT_Luminal_vs_Basal_summary.txt")
path_ccle <- file.path(repo, "data", "results", "CCLE_corrected", "DA_luminal_vs_basal", "DA_luminal_vs_basal_limma.csv")
path_ccle_sum <- file.path(repo, "data", "results", "CCLE", "ccle_sum", "DA_luminal_vs_basal_table_s2", "DA_luminal_vs_basal_limma.csv")

if (!file.exists(path_cptac)) stop("Missing CPTAC file: ", path_cptac)
if (!file.exists(path_ccle_sum)) stop("Missing CCLE Table S2 file (primary CCLE for overlap): ", path_ccle_sum)
if (!file.exists(path_ccle)) warning("CCLE_corrected limma not found (optional 3-row table): ", path_ccle)

# ---- Helpers ----
extract_gene_uniprot <- function(uid) {
  uid <- as.character(uid)
  g <- sub("^[^|]*\\|[^|]*\\|([A-Z0-9]+)(?:-[0-9]+)?_HUMAN.*", "\\1", uid, perl = TRUE)
  ifelse(grepl("^[A-Z0-9]+$", g), g, NA_character_)
}

load_cptac_dedup <- function() {
  cpt <- fread(path_cptac, showProgress = FALSE)
  setDT(cpt)
  adj_col <- if ("adj.pvalue" %in% names(cpt)) "adj.pvalue" else "adj.P.Val"
  cpt[, gene := trimws(as.character(Gene_symbol))]
  cpt <- cpt[nzchar(gene) & !is.na(gene)]
  cpt <- cpt[order(get(adj_col))][, .SD[1], by = gene]
  cpt
}

load_ccle_limma_dedup <- function(path, id_col) {
  cc <- fread(path, showProgress = FALSE)
  setDT(cc)
  cc[, gene := extract_gene_uniprot(get(id_col))]
  cc <- cc[!is.na(gene)]
  cc <- cc[order(adj.P.Val)][, .SD[1], by = gene]
  cc
}

load_cptac_sets <- function(cpt) {
  adj_col <- if ("adj.pvalue" %in% names(cpt)) "adj.pvalue" else "adj.P.Val"
  cpt_sig <- cpt[get(adj_col) < FDR_THR]
  list(
    lum = sort(unique(cpt_sig[log2FC > 0, gene])),
    bas = sort(unique(cpt_sig[log2FC < 0, gene]))
  )
}

load_ccle_sets <- function(cc) {
  cc_sig <- cc[adj.P.Val < FDR_THR]
  list(
    lum = sort(unique(cc_sig[logFC > 0, gene])),
    bas = sort(unique(cc_sig[logFC < 0, gene]))
  )
}

cptac_sample_sizes <- function() {
  lum_n <- NA_integer_
  bas_n <- NA_integer_
  if (file.exists(path_cptac_sum)) {
    sl <- readLines(path_cptac_sum, warn = FALSE)
    i1 <- grep("Luminal:", sl)[1]
    if (!is.na(i1)) lum_n <- as.integer(sub(".*:\\s*([0-9]+).*", "\\1", sl[i1]))
    i2 <- grep("^\\s*Basal:", sl)[1]
    if (!is.na(i2)) bas_n <- as.integer(sub(".*:\\s*([0-9]+).*", "\\1", sl[i2]))
  }
  c(lum_n, bas_n)
}

jacc <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

fc_direction_label <- function(fc) {
  # Positive = higher in Luminal (Luminal - Basal)
  ifelse(!is.na(fc) & fc > 0, "Luminal-higher", ifelse(!is.na(fc) & fc < 0, "Basal-higher", NA_character_))
}

sig_label <- function(adj, thr = FDR_THR) {
  ifelse(!is.na(adj) & adj < thr, "yes", "no")
}

ccle_direction_display <- function(fc) {
  d <- fc_direction_label(fc)
  fifelse(is.na(fc), "not in CCLE matrix", d)
}

ccle_sig_display <- function(adj) {
  fifelse(is.na(adj), "n/a", sig_label(adj))
}

# ---- Load ----
cpt_full <- load_cptac_dedup()
# Primary CCLE for slides: Table S2 (more genes / overlap with CPTAC than corrected-only matrix).
cc_full <- load_ccle_limma_dedup(path_ccle_sum, "Protein_Id")
cc_corr_full <- if (file.exists(path_ccle)) {
  load_ccle_limma_dedup(path_ccle, "UniProtID")
} else {
  NULL
}

cpt <- load_cptac_sets(cpt_full)
cc <- load_ccle_sets(cc_full)
cpt_lum <- cpt$lum
cpt_bas <- cpt$bas
cc_lum <- cc$lum
cc_bas <- cc$bas

cc_corr <- if (!is.null(cc_corr_full)) load_ccle_sets(cc_corr_full) else list(lum = character(0), bas = character(0))
cc_corr_lum <- cc_corr$lum
cc_corr_bas <- cc_corr$bas

lum_bas_n <- cptac_sample_sizes()
LUM_CCLE <- 4L
BAS_CCLE <- 4L

adj_cpt <- if ("adj.pvalue" %in% names(cpt_full)) "adj.pvalue" else "adj.P.Val"

# ---- 1) Slide summary table: CPTAC + CCLE (Table S2) ----
summ2 <- data.table(
  Dataset = c("CPTAC", "CCLE_Table_S2"),
  Luminal_n = c(lum_bas_n[1], LUM_CCLE),
  Basal_n = c(lum_bas_n[2], BAS_CCLE),
  Luminal_up = c(length(cpt_lum), length(cc_lum)),
  Basal_up = c(length(cpt_bas), length(cc_bas)),
  Total_DA = c(length(cpt_lum) + length(cpt_bas), length(cc_lum) + length(cc_bas))
)
fwrite(summ2, file.path(repo, "reports", "subtype_summary_table.csv"))

# Extended table (corrected vs Table S2) when corrected limma exists
summ3 <- if (!is.null(cc_corr_full)) {
  data.table(
    Dataset = c("CPTAC", "CCLE_corrected", "CCLE_Table_S2"),
    Luminal_n = c(lum_bas_n[1], LUM_CCLE, LUM_CCLE),
    Basal_n = c(lum_bas_n[2], BAS_CCLE, BAS_CCLE),
    Luminal_up = c(length(cpt_lum), length(cc_corr_lum), length(cc_lum)),
    Basal_up = c(length(cpt_bas), length(cc_corr_bas), length(cc_bas)),
    Total_DA = c(
      length(cpt_lum) + length(cpt_bas),
      length(cc_corr_lum) + length(cc_corr_bas),
      length(cc_lum) + length(cc_bas)
    )
  )
} else {
  copy(summ2)
}
fwrite(summ3, file.path(repo, "reports", "subtype_cptac_ccle_summary.csv"))

# ---- Gene lists ----
writeLines(cpt_lum, file.path(repo, "reports", "cptac_luminal_up_genes.txt"))
writeLines(cpt_bas, file.path(repo, "reports", "cptac_basal_up_genes.txt"))
writeLines(cc_lum, file.path(repo, "reports", "ccle_luminal_up_genes.txt"))
writeLines(cc_bas, file.path(repo, "reports", "ccle_basal_up_genes.txt"))
writeLines(cc_lum, file.path(repo, "reports", "ccle_sum_luminal_up_genes.txt"))
writeLines(cc_bas, file.path(repo, "reports", "ccle_sum_basal_up_genes.txt"))
if (!is.null(cc_corr_full)) {
  writeLines(cc_corr_lum, file.path(repo, "reports", "ccle_corrected_luminal_up_genes.txt"))
  writeLines(cc_corr_bas, file.path(repo, "reports", "ccle_corrected_basal_up_genes.txt"))
}

# ---- 2) Metadata for slides ----
meta_txt <- c(
  "Subtype benchmark — definitions and file metadata",
  "===============================================",
  "",
  "Threshold (same for CPTAC and CCLE)",
  paste0("  FDR < ", FDR_THR, " on the adjusted p-value column."),
  "",
  "Sign convention (both cohorts)",
  "  Positive logFC / log2FC = Luminal minus Basal (higher abundance in Luminal when > 0).",
  "",
  "CPTAC",
  paste("  File:", path_cptac),
  paste("  Luminal samples (n):", lum_bas_n[1], "; Basal samples (n):", lum_bas_n[2]),
  "  Gene column: Gene_symbol (trimmed; one row per gene after dedup by min FDR)",
  "  Effect: log2FC",
  paste("  FDR column:", adj_cpt),
  "  Luminal-up set: FDR < threshold AND log2FC > 0",
  "  Basal-up set:   FDR < threshold AND log2FC < 0",
  "",
  "CCLE — primary for slides: Table S2 proteome (limma on Luminal vs Basal cell lines)",
  paste("  File:", path_ccle_sum),
  paste("  Luminal cell lines (n):", LUM_CCLE, "; Basal cell lines (n):", BAS_CCLE),
  "  Gene: parsed from Protein_Id (sp|ACCESSION|SYMBOL_HUMAN → SYMBOL)",
  "  Effect: logFC",
  "  FDR column: adj.P.Val",
  "  Luminal-up set: FDR < threshold AND logFC > 0",
  "  Basal-up set:   FDR < threshold AND logFC < 0",
  "",
  if (file.exists(path_ccle)) {
    c(
      "CCLE — optional comparison: corrected in-house matrix (smaller overlap with CPTAC)",
      paste("  File:", path_ccle),
      "  Gene column: UniProtID (same parsing rule)",
      ""
    )
  } else {
    character(0)
  },
  "Total_DA = Luminal_up + Basal_up (non-overlapping directional DE genes).",
  ""
)
writeLines(meta_txt, file.path(repo, "reports", "subtype_benchmark_metadata.txt"))

# ---- 3) Venn: diagram-only (two figures, by direction) ----
plot_venn_diagram_only <- function(s1, s2, lab1, lab2, main_title, path_png) {
  if (length(s1) == 0L && length(s2) == 0L) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Both DE sets empty", size = 5) +
      theme_void() +
      labs(title = main_title) +
      theme(plot.margin = margin(12, 28, 12, 28, "pt"))
    ggsave(path_png, p, width = VENN_DIAGRAM_WIDTH, height = VENN_DIAGRAM_HEIGHT, dpi = 300, bg = "white")
    return(invisible(NULL))
  }
  vl <- stats::setNames(list(s1, s2), c(lab1, lab2))
  p <- ggVennDiagram::ggVennDiagram(
    vl,
    label = "count",
    label_alpha = 0.9,
    edge_size = 0.55
  ) +
    scale_fill_gradient(low = "#f2f6ff", high = "#2d5a9e") +
    labs(
      title = main_title,
      subtitle = paste0("FDR < ", FDR_THR, " | counts = distinct gene symbols | CCLE = Table S2")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.margin = margin(14, 32, 14, 32, "pt"),
      legend.margin = margin(4, 8, 4, 8)
    )
  ggsave(path_png, p, width = VENN_DIAGRAM_WIDTH, height = VENN_DIAGRAM_HEIGHT, dpi = 300, bg = "white")
  invisible(NULL)
}

plot_venn_diagram_only(
  cpt_lum, cc_lum, "CPTAC", "CCLE (Table S2)",
  "Luminal-up genes (shared subtype direction)",
  file.path(repo, "reports", "venn_luminal_up_cptac_vs_ccle.png")
)
plot_venn_diagram_only(
  cpt_bas, cc_bas, "CPTAC", "CCLE (Table S2)",
  "Basal-up genes (shared subtype direction)",
  file.path(repo, "reports", "venn_basal_up_cptac_vs_ccle.png")
)

# ---- Optional: combined diagram + gene lists (not for minimal slide deck) ----
fmt_block <- function(title, genes, max_show = NULL) {
  if (length(genes) == 0) return(paste0(title, ": (none)"))
  g <- sort(unique(genes))
  if (!is.null(max_show) && length(g) > max_show) {
    paste0(
      title, " (n=", length(g), "): ",
      paste(g[seq_len(max_show)], collapse = ", "),
      "\n    ... +", length(g) - max_show, " more"
    )
  } else {
    paste0(title, " (n=", length(g), "): ", paste(g, collapse = ", "))
  }
}

plot_venn_with_lists <- function(s1, s2, lab1, lab2, main_title, path_png) {
  only1 <- setdiff(s1, s2)
  only2 <- setdiff(s2, s1)
  both <- intersect(s1, s2)
  txt <- paste(
    fmt_block("Intersection", both),
    "",
    fmt_block(paste0(lab1, " only"), only1, MAX_GENES_SHOW_CPTAC_ONLY),
    "",
    fmt_block(paste0(lab2, " only"), only2, if (length(only2) > 55L) 50L else NULL),
    sep = "\n"
  )
  txt_wrapped <- paste(strwrap(txt, width = 96), collapse = "\n")
  if (length(s1) == 0L && length(s2) == 0L) {
    plot_venn_diagram_only(s1, s2, lab1, lab2, main_title, path_png)
    return(invisible(NULL))
  }
  vl <- stats::setNames(list(s1, s2), c(lab1, lab2))
  p1 <- ggVennDiagram::ggVennDiagram(vl, label = "count", label_alpha = 0.88, edge_size = 0.55) +
    scale_fill_gradient(low = "#f0f4ff", high = "#3b6fb6") +
    labs(title = main_title, subtitle = paste0("FDR < ", FDR_THR)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 13))
  p2 <- ggplot() +
    annotate("text", x = 0, y = 1, label = txt_wrapped, hjust = 0, vjust = 1, size = 2.8, lineheight = 1.06) +
    xlim(0, 1) + ylim(0, 1.05) + theme_void() +
    labs(caption = "CPTAC-only lists truncated; see reports/cptac_*_genes.txt")
  ggsave(path_png, p1 / p2 + plot_layout(heights = c(0.45, 0.55)), width = VENN_WITH_LISTS_WIDTH, height = 13, dpi = 300, bg = "white")
}

plot_venn_with_lists(cpt_lum, cc_lum, "CPTAC", "CCLE (Table S2)", "Luminal-up (with gene lists)", file.path(repo, "reports", "venn_luminal_up_cptac_vs_ccle_with_gene_lists.png"))
plot_venn_with_lists(cpt_bas, cc_bas, "CPTAC", "CCLE (Table S2)", "Basal-up (with gene lists)", file.path(repo, "reports", "venn_basal_up_cptac_vs_ccle_with_gene_lists.png"))

if (!is.null(cc_corr_full) && (length(cc_corr_lum) + length(cc_corr_bas) > 0 || length(cpt_lum) + length(cpt_bas) > 0)) {
  plot_venn_diagram_only(
    cpt_lum, cc_corr_lum, "CPTAC", "CCLE (corrected)",
    "Luminal-up — CPTAC vs CCLE_corrected (smaller overlap)",
    file.path(repo, "reports", "venn_luminal_up_cptac_vs_ccle_corrected.png")
  )
  plot_venn_diagram_only(
    cpt_bas, cc_corr_bas, "CPTAC", "CCLE (corrected)",
    "Basal-up — CPTAC vs CCLE_corrected (smaller overlap)",
    file.path(repo, "reports", "venn_basal_up_cptac_vs_ccle_corrected.png")
  )
}

# ---- 4) Marker sanity table (direction first, significance second) ----
markers_luminal <- c("ESR1", "GATA3", "FOXA1", "KRT18", "PGR")
markers_basal <- c("KRT5", "KRT14", "KRT17", "EGFR", "FOXC1")
mk_genes <- c(markers_luminal, markers_basal)
expected <- c(rep("Luminal", 5L), rep("Basal", 5L))

cpt_key <- cpt_full[, .(gene, log2FC, adj = get(adj_cpt))]
setnames(cpt_key, "log2FC", "cpt_fc")
cc_key <- cc_full[, .(gene, logFC, adj.P.Val)]
setnames(cc_key, c("logFC", "adj.P.Val"), c("cc_fc", "cc_adj"))

ms <- data.table(Gene = mk_genes, Expected_subtype = expected)
ms <- merge(ms, cpt_key, by.x = "Gene", by.y = "gene", all.x = TRUE)
ms <- merge(ms, cc_key, by.x = "Gene", by.y = "gene", all.x = TRUE)

ms[, `:=`(
  CPTAC_direction = {
    d <- fc_direction_label(cpt_fc)
    fifelse(is.na(cpt_fc), "not in CPTAC table", d)
  },
  CPTAC_sig = fifelse(is.na(adj), "n/a", sig_label(adj)),
  CCLE_direction = ccle_direction_display(cc_fc),
  CCLE_sig = ccle_sig_display(cc_adj)
)]

ms[, interpretation := {
  ok_lum <- Expected_subtype == "Luminal" & CPTAC_direction == "Luminal-higher"
  ok_bas <- Expected_subtype == "Basal" & CPTAC_direction == "Basal-higher"
  fifelse(
    Expected_subtype == "Luminal",
    fifelse(ok_lum & CPTAC_sig == "yes", "CPTAC: strong (dir+sig)",
      fifelse(ok_lum, "CPTAC: direction OK, not sig",
        "CPTAC: wrong or missing")),
    fifelse(ok_bas & CPTAC_sig == "yes", "CPTAC: strong (dir+sig)",
      fifelse(ok_bas, "CPTAC: direction OK, not sig",
        "CPTAC: wrong or missing"))
  )
}]

# Simpler output columns for slides (user spec)
slide_mk <- ms[, .(
  Gene,
  Expected_subtype,
  CPTAC_direction,
  CPTAC_sig,
  CCLE_direction,
  CCLE_sig
)]
fwrite(slide_mk, file.path(repo, "reports", "subtype_marker_sanity.csv"))

# Narrative interpretation (CCLE column notes)
ms[, CCLE_note := fifelse(
  Expected_subtype == "Luminal" & CCLE_direction == "Luminal-higher" & CCLE_sig == "no",
  "direction supportive; low n → often not sig",
  fifelse(
    Expected_subtype == "Basal" & CCLE_direction == "Basal-higher" & CCLE_sig == "no",
    "direction supportive; low n → often not sig",
    fifelse(CCLE_sig == "yes" & (
      (Expected_subtype == "Luminal" & CCLE_direction == "Luminal-higher") |
      (Expected_subtype == "Basal" & CCLE_direction == "Basal-higher")
    ), "direction + sig", "")
  )
)]
fwrite(ms[, .(Gene, Expected_subtype, CPTAC_direction, CPTAC_sig, CCLE_direction, CCLE_sig, interpretation, CCLE_note)],
  file.path(repo, "reports", "subtype_marker_sanity_extended.csv")
)

# ---- 5) Slide highlights (talking points) ----
ov_lum <- intersect(cpt_lum, cc_lum)
ov_bas <- intersect(cpt_bas, cc_bas)
highlights <- c(
  "Slide highlights — subtype benchmark (CPTAC vs CCLE Table S2)",
  "=============================================================",
  "",
  "From the summary table",
  "- Both cohorts report many subtype-associated genes at FDR < 0.05 (CPTAC typically stronger due to n).",
  "- CCLE row uses Table S2 proteome (more overlap with CPTAC than the corrected-only matrix).",
  "- CCLE uses 4 Luminal + 4 Basal lines; interpret counts as exploratory.",
  "",
  "From the two Venns (Luminal-up separate from Basal-up)",
  "- Overlap is partial: shared subtype core exists but is not identical (different tissue/TMT context).",
  "- Direction is separated: do not mix Luminal-up and Basal-up in one Venn.",
  "",
  "From the marker sanity table",
  "- Read direction first, significance second: correct direction without FDR is still supportive with n=4.",
  "- Wrong direction (e.g. Luminal marker higher in Basal) deserves attention even if not significant.",
  "",
  paste("Overlap genes — Luminal-up:", paste(sort(ov_lum), collapse = ", ")),
  paste("Overlap genes — Basal-up:", paste(sort(ov_bas), collapse = ", ")),
  ""
)
writeLines(highlights, file.path(repo, "reports", "subtype_slide_highlights.txt"))

# ---- Legacy summary text ----
raw_ct <- {
  s <- cpt_full[get(adj_cpt) < FDR_THR]
  c(sum(s$log2FC > 0, na.rm = TRUE), sum(s$log2FC < 0, na.rm = TRUE))
}
txt <- c(
  "CPTAC vs CCLE — subtype benchmark (see subtype_summary_table.csv for 2-row slide table)",
  "",
  paste("FDR:", FDR_THR),
  paste("CPTAC file:", path_cptac),
  paste("CCLE Table S2 file:", path_ccle_sum),
  paste("CPTAC protein rows at FDR (before gene dedup): Lum+", raw_ct[1], "Bas+", raw_ct[2]),
  "",
  "--- 2-row slide table ---",
  paste(capture.output(print(summ2, row.names = FALSE)), collapse = "\n"),
  "",
  "Intersections (gene symbols):",
  paste("Luminal-up:", paste(sort(ov_lum), collapse = ", ")),
  paste("Basal-up:", paste(sort(ov_bas), collapse = ", ")),
  "",
  "Files: subtype_benchmark_metadata.txt, subtype_marker_sanity.csv, subtype_slide_highlights.txt",
  "Venns (diagram-only): venn_luminal_up_cptac_vs_ccle.png, venn_basal_up_cptac_vs_ccle.png",
  "Repro: Rscript --no-init-file data/scripts/subtype_cptac_ccle_benchmark.R"
)
writeLines(txt, file.path(repo, "reports", "subtype_cptac_ccle_summary.txt"))

ov_txt <- c(
  "Overlap (CPTAC vs CCLE Table S2 — primary)",
  "===========================================",
  paste("Luminal-up Jaccard:", round(jacc(cpt_lum, cc_lum), 5)),
  paste("Basal-up Jaccard:", round(jacc(cpt_bas, cc_bas), 5)),
  ""
)
if (!is.null(cc_corr_full)) {
  ov_txt <- c(
    ov_txt,
    "Overlap (CPTAC vs CCLE_corrected — optional)",
    "===========================================",
    paste("Luminal-up Jaccard:", round(jacc(cpt_lum, cc_corr_lum), 5)),
    paste("Basal-up Jaccard:", round(jacc(cpt_bas, cc_corr_bas), 5))
  )
}
writeLines(ov_txt, file.path(repo, "reports", "subtype_overlap_summary.txt"))

message("Wrote: reports/subtype_summary_table.csv, subtype_benchmark_metadata.txt, venn_* (diagram-only), subtype_marker_sanity*.csv, subtype_slide_highlights.txt")
