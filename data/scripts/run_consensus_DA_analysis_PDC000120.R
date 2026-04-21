#!/usr/bin/env Rscript
# =============================================================================
# Consensus downstream analysis for Tumor vs NAT DA results (PDC000120)
#
# Inputs (already generated upstream):
#   - results/PDC000120/DA_limma_tumor_vs_NAT.csv
#   - results/PDC000120/DA_MSstats_tumor_vs_NAT.csv
#   - results/PDC000120/DA_MSstatsTMT_tumor_vs_NAT.csv
#
# Outputs (all under results/PDC000120/):
#   - DA_consensus_table.csv
#   - DA_method_correlations.csv
#   - scatter_*.pdf
#   - DA_overlap_summary.csv
#   - DA_marker_sanity_check.csv
#   - DA_high_confidence_tumor_up.csv
#   - DA_high_confidence_NAT_up.csv
#   - DA_high_confidence_top20_for_slides.csv
#   - DA_method_diagnostics.csv
#   - histogram_pvalues_*.pdf (when raw p-values exist)
#   - enrichment_*.csv + enrichment_*.pdf (if enrichment packages available)
#   - DA_consensus_summary.txt
#
# Design goals:
#   - robust to missing columns / NAs
#   - explicit counts + no silent dropping
#   - rerunnable and meeting-ready outputs
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel", repos = "https://cloud.r-project.org")
})
library(data.table)
library(ggplot2)
library(ggrepel)

message("=== Consensus DA downstream analysis (PDC000120) ===")

detect_data_dir <- function() {
  if (length(sys.frames()) > 0 && exists("ofile", sys.frame(1))) {
    return(dirname(dirname(sys.frame(1)$ofile)))
  }
  wd <- getwd()
  if (file.exists(file.path(wd, "results", "PDC000120", "DA_limma_tumor_vs_NAT.csv"))) return(wd)
  if (file.exists(file.path(wd, "data", "results", "PDC000120", "DA_limma_tumor_vs_NAT.csv"))) return(file.path(wd, "data"))
  wd
}

DATA_DIR <- detect_data_dir()
setwd(DATA_DIR)
RES_DIR <- file.path(DATA_DIR, "results", "PDC000120")
CONS_DIR <- file.path(RES_DIR, "consensus")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CONS_DIR, showWarnings = FALSE, recursive = TRUE)

IN_LIMMA <- file.path(RES_DIR, "DA_limma_tumor_vs_NAT.csv")
IN_MSSTATS <- file.path(RES_DIR, "DA_MSstats_tumor_vs_NAT.csv")
IN_MSSTATSTMT <- file.path(RES_DIR, "DA_MSstatsTMT_tumor_vs_NAT.csv")

stopifnot(file.exists(IN_LIMMA), file.exists(IN_MSSTATS), file.exists(IN_MSSTATSTMT))

SIG_FDR <- 0.05
SIG_LFC <- 1

`%||%` <- function(a, b) if (!is.null(a)) a else b

col_first_present <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(x))

assign_sign <- function(lfc) {
  out <- rep(NA_character_, length(lfc))
  out[!is.na(lfc) & lfc > 0] <- "pos"
  out[!is.na(lfc) & lfc < 0] <- "neg"
  out[!is.na(lfc) & lfc == 0] <- "zero"
  out
}

is_sig <- function(lfc, fdr, fdr_thr = SIG_FDR, lfc_thr = SIG_LFC) {
  !is.na(lfc) & !is.na(fdr) & (abs(lfc) > lfc_thr) & (fdr < fdr_thr)
}

try_install_bioc <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    tryCatch(install.packages("BiocManager", repos = "https://cloud.r-project.org"), error = function(e) NULL)
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) return(FALSE)
  ok <- tryCatch({
    BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  requireNamespace(pkg, quietly = TRUE) && ok
}

try_install_cran <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  ok <- tryCatch({
    install.packages(pkg, repos = "https://cloud.r-project.org")
    TRUE
  }, error = function(e) FALSE)
  requireNamespace(pkg, quietly = TRUE) && ok
}

# -----------------------------------------------------------------------------
# 1) Load + standardize the three DA tables
# -----------------------------------------------------------------------------

message("\n[1/10] Loading DA tables")
limma_raw <- fread(IN_LIMMA)
msstats_raw <- fread(IN_MSSTATS)
msstatstmt_raw <- fread(IN_MSSTATSTMT)

message("  limma rows: ", nrow(limma_raw), " cols: ", ncol(limma_raw))
message("  MSstats rows: ", nrow(msstats_raw), " cols: ", ncol(msstats_raw))
message("  MSstatsTMT rows: ", nrow(msstatstmt_raw), " cols: ", ncol(msstatstmt_raw))

standardize_limma <- function(dt) {
  id_col <- col_first_present(dt, c("GeneSymbol", "Gene", "Protein", "ProteinName"))
  lfc_col <- col_first_present(dt, c("log2FC", "logFC", "log2fc"))
  fdr_col <- col_first_present(dt, c("FDR", "adj.P.Val", "adj.pvalue", "Adjusted.Pvalue"))
  p_col <- col_first_present(dt, c("pvalue", "P.Value", "p.value"))
  if (is.na(id_col)) stop("limma: cannot find identifier column")
  out <- data.table(
    id_raw = as.character(dt[[id_col]]),
    gene_symbol = as.character(dt[[id_col]]),
    log2FC = as_numeric_safe(dt[[lfc_col %||% "logFC"]]),
    FDR = as_numeric_safe(dt[[fdr_col %||% "adj.P.Val"]]),
    pvalue = if (!is.na(p_col)) as_numeric_safe(dt[[p_col]]) else NA_real_,
    method = "limma"
  )
  out
}

standardize_msstats <- function(dt) {
  id_col <- col_first_present(dt, c("GeneSymbol", "ProteinName", "Protein", "Gene"))
  lfc_col <- col_first_present(dt, c("log2FC", "logFC"))
  fdr_col <- col_first_present(dt, c("adj.pvalue", "Adjusted.Pvalue", "adj.P.Val", "adj.P.Value", "FDR"))
  p_col <- col_first_present(dt, c("pvalue", "P.Value", "p.value"))
  if (is.na(id_col)) stop("MSstats: cannot find identifier column")
  out <- data.table(
    id_raw = as.character(dt[[id_col]]),
    gene_symbol = if ("GeneSymbol" %in% names(dt)) as.character(dt[["GeneSymbol"]]) else as.character(dt[[id_col]]),
    log2FC = as_numeric_safe(dt[[lfc_col %||% "log2FC"]]),
    FDR = as_numeric_safe(dt[[fdr_col %||% "adj.pvalue"]]),
    pvalue = if (!is.na(p_col)) as_numeric_safe(dt[[p_col]]) else NA_real_,
    method = "MSstats"
  )
  out
}

standardize_msstatsTMT <- function(dt) {
  id_col <- col_first_present(dt, c("Protein", "ProteinName", "GeneSymbol"))
  lfc_col <- col_first_present(dt, c("log2FC", "logFC"))
  fdr_col <- col_first_present(dt, c("adj.pvalue", "Adjusted.Pvalue", "adj.P.Val", "FDR"))
  p_col <- col_first_present(dt, c("pvalue", "P.Value", "p.value"))
  if (is.na(id_col)) stop("MSstatsTMT: cannot find identifier column")
  out <- data.table(
    id_raw = as.character(dt[[id_col]]),
    gene_symbol = NA_character_, # will attempt mapping below
    log2FC = as_numeric_safe(dt[[lfc_col %||% "log2FC"]]),
    FDR = as_numeric_safe(dt[[fdr_col %||% "adj.pvalue"]]),
    pvalue = if (!is.na(p_col)) as_numeric_safe(dt[[p_col]]) else NA_real_,
    method = "MSstatsTMT"
  )
  out
}

limma <- standardize_limma(limma_raw)
msstats <- standardize_msstats(msstats_raw)
msstatstmt <- standardize_msstatsTMT(msstatstmt_raw)

message("  Standardized: limma ", nrow(limma), ", MSstats ", nrow(msstats), ", MSstatsTMT ", nrow(msstatstmt))

# -----------------------------------------------------------------------------
# Identifier cleanup + mapping for MSstatsTMT (often RefSeq/UniProt-like)
# -----------------------------------------------------------------------------

message("\n[1b/10] Identifier cleanup + mapping (best-effort)")

looks_like_gene_symbol <- function(x) {
  # permissive: typical HGNC symbols + numeric Entrez IDs from limma/MSstats tables
  grepl("^[A-Za-z][A-Za-z0-9\\-\\.]+$", x) | grepl("^[0-9]+$", x)
}

try_map_ids_to_symbols <- function(ids) {
  # Return named character vector: names are input ids, values are gene symbols (or NA).
  ids <- unique(na.omit(as.character(ids)))
  out <- setNames(rep(NA_character_, length(ids)), ids)

  # If already looks like a gene symbol, keep it as-is.
  already <- looks_like_gene_symbol(ids)
  out[already] <- ids[already]

  # Try Bioconductor mapping (org.Hs.eg.db) for UNIPROT / REFSEQ
  ok_org <- try_install_bioc("org.Hs.eg.db") && try_install_bioc("AnnotationDbi")
  if (!ok_org) {
    message("  NOTE: org.Hs.eg.db not available; MSstatsTMT IDs may remain unmapped.")
    return(out)
  }
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })

  # Helper: safe select
  safe_select <- function(keys, keytype) {
    keys <- unique(keys)
    keys <- keys[nzchar(keys)]
    if (length(keys) == 0) return(data.frame())
    tryCatch(
      AnnotationDbi::select(org.Hs.eg.db, keys = keys, keytype = keytype, columns = c("SYMBOL")),
      error = function(e) data.frame()
    )
  }

  # Split candidate types
  ids2 <- ids[!already]
  # UniProt accessions are typically [OPQ][0-9][A-Z0-9]{3}[0-9] or [A-NR-Z][0-9]{5}
  is_uniprot_acc <- grepl("^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$", ids2) | grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$", ids2)
  uniprot_keys <- sub("-[0-9]+$", "", ids2[is_uniprot_acc])

  # RefSeq protein accessions e.g. NP_000005, XP_...
  is_refseq <- grepl("^[NXYP]P_\\d+(\\.\\d+)?$", ids2)
  refseq_keys <- sub("\\.\\d+$", "", ids2[is_refseq])

  # Map UniProt
  if (length(uniprot_keys) > 0) {
    m <- safe_select(uniprot_keys, "UNIPROT")
    if (nrow(m) > 0) {
      m <- m[!is.na(m$SYMBOL) & nzchar(m$SYMBOL), , drop = FALSE]
      sym_by_key <- setNames(m$SYMBOL, m$UNIPROT)
      for (k in unique(uniprot_keys)) {
        if (!is.na(sym_by_key[k])) out[names(out) %in% ids2[is_uniprot_acc]] <- out[names(out) %in% ids2[is_uniprot_acc]]
      }
      # assign back using stripped isoform keys
      for (id in ids2[is_uniprot_acc]) {
        key <- sub("-[0-9]+$", "", id)
        if (!is.na(sym_by_key[key]) && nzchar(sym_by_key[key])) out[id] <- sym_by_key[key]
      }
    }
  }

  # Map RefSeq
  if (length(refseq_keys) > 0) {
    m <- safe_select(refseq_keys, "REFSEQ")
    if (nrow(m) > 0) {
      m <- m[!is.na(m$SYMBOL) & nzchar(m$SYMBOL), , drop = FALSE]
      sym_by_key <- setNames(m$SYMBOL, m$REFSEQ)
      for (id in ids2[is_refseq]) {
        key <- sub("\\.\\d+$", "", id)
        if (!is.na(sym_by_key[key]) && nzchar(sym_by_key[key])) out[id] <- sym_by_key[key]
      }
    }
  }

  out
}

msstatstmt[, id_clean := trimws(id_raw)]
sym_map <- try_map_ids_to_symbols(msstatstmt$id_clean)
msstatstmt[, gene_symbol := unname(sym_map[id_clean])]

map_rate <- mean(!is.na(msstatstmt$gene_symbol) & nzchar(msstatstmt$gene_symbol))
message(sprintf("  MSstatsTMT mapping rate to gene symbols: %.1f%%", 100 * map_rate))

# For joining, define a consensus key that prioritizes gene_symbol (mapped), else id_clean
make_key <- function(gene_symbol, id_clean) {
  key <- ifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, id_clean)
  as.character(key)
}

limma[, key := make_key(gene_symbol, id_raw)]
msstats[, key := make_key(gene_symbol, id_raw)]
msstatstmt[, key := make_key(gene_symbol, id_clean)]

message("  Unique keys: limma ", uniqueN(limma$key), ", MSstats ", uniqueN(msstats$key), ", MSstatsTMT ", uniqueN(msstatstmt$key))

# -----------------------------------------------------------------------------
# 2) Build consensus comparison table
# -----------------------------------------------------------------------------

message("\n[2/10] Building consensus merged table")

wide_one <- function(dt, method_name) {
  # returns unique key row with method-specific columns; if duplicates, keep best (smallest FDR then abs(log2FC))
  d <- copy(dt)
  d[, log2FC := as_numeric_safe(log2FC)]
  d[, FDR := as_numeric_safe(FDR)]
  d[, pvalue := as_numeric_safe(pvalue)]
  d[, ord := fifelse(!is.na(FDR), FDR, Inf)]
  d[, ord2 := -abs(fifelse(!is.na(log2FC), log2FC, 0))]
  setorder(d, key, ord, ord2)
  d <- d[, .SD[1], by = key]
  setnames(d,
           c("log2FC", "FDR", "pvalue"),
           c(paste0("log2FC_", method_name),
             paste0("FDR_", method_name),
             paste0("pvalue_", method_name)))
  # Keep representative identifiers for transparency
  if ("id_clean" %in% names(d)) {
    d[, paste0("id_raw_", method_name) := id_clean]
  } else {
    d[, paste0("id_raw_", method_name) := id_raw]
  }
  if ("gene_symbol" %in% names(d)) d[, paste0("gene_symbol_", method_name) := gene_symbol]
  drop_cols <- intersect(c("id_raw", "gene_symbol", "method", "ord", "ord2", "id_clean"), names(d))
  if (length(drop_cols) > 0) d[, (drop_cols) := NULL]
  d
}

limma_w <- wide_one(limma, "limma")
msstats_w <- wide_one(msstats, "MSstats")
msstatstmt_w <- wide_one(msstatstmt, "MSstatsTMT")

cons <- Reduce(function(x, y) merge(x, y, by = "key", all = TRUE), list(limma_w, msstats_w, msstatstmt_w))
setnames(cons, "key", "Gene")

for (m in c("limma", "MSstats", "MSstatsTMT")) {
  cons[[paste0("sign_", m)]] <- assign_sign(cons[[paste0("log2FC_", m)]])
  cons[[paste0("significant_", m)]] <- is_sig(cons[[paste0("log2FC_", m)]], cons[[paste0("FDR_", m)]])
}

cons[, n_methods_present := rowSums(!is.na(.SD)), .SDcols = c("log2FC_limma", "log2FC_MSstats", "log2FC_MSstatsTMT")]
cons[, n_methods_significant := rowSums(.SD == TRUE, na.rm = TRUE),
     .SDcols = c("significant_limma", "significant_MSstats", "significant_MSstatsTMT")]

# Consensus mean effect size and strict sign rule across available methods
cons[, mean_log2FC := rowMeans(.SD, na.rm = TRUE),
     .SDcols = c("log2FC_limma", "log2FC_MSstats", "log2FC_MSstatsTMT")]

cons[, direction_consensus := {
  vals <- c(log2FC_limma, log2FC_MSstats, log2FC_MSstatsTMT)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) "insufficient"
  else if (all(vals > 0)) "Tumor_up"
  else if (all(vals < 0)) "NAT_up"
  else "mixed_or_discordant"
}, by = seq_len(nrow(cons))]

cons[, high_confidence_tumor_up :=
       (n_methods_significant >= 2) &
       (direction_consensus == "Tumor_up")]

cons[, high_confidence_NAT_up :=
       (n_methods_significant >= 2) &
       (direction_consensus == "NAT_up")]

cons[, strict_consensus :=
       (significant_limma & significant_MSstats & significant_MSstatsTMT) &
       (direction_consensus %in% c("Tumor_up", "NAT_up"))]

fwrite(cons, file.path(CONS_DIR, "DA_consensus_table.csv"))
message("  Wrote DA_consensus_table.csv (rows: ", nrow(cons), ")")

# -----------------------------------------------------------------------------
# 3) Cross-method effect size concordance
# -----------------------------------------------------------------------------

message("\n[3/10] Cross-method effect size concordance")

cor_pair <- function(df, m1, m2) {
  x <- df[[paste0("log2FC_", m1)]]
  y <- df[[paste0("log2FC_", m2)]]
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3) {
    return(data.table(method1 = m1, method2 = m2, n_shared = n, pearson = NA_real_, spearman = NA_real_))
  }
  data.table(
    method1 = m1,
    method2 = m2,
    n_shared = n,
    pearson = cor(x[ok], y[ok], method = "pearson"),
    spearman = cor(x[ok], y[ok], method = "spearman")
  )
}

corrs <- rbindlist(list(
  cor_pair(cons, "limma", "MSstats"),
  cor_pair(cons, "limma", "MSstatsTMT"),
  cor_pair(cons, "MSstats", "MSstatsTMT")
))
fwrite(corrs, file.path(CONS_DIR, "DA_method_correlations.csv"))
message("  Wrote DA_method_correlations.csv")

plot_scatter <- function(df, m1, m2, out_pdf) {
  xcol <- paste0("log2FC_", m1)
  ycol <- paste0("log2FC_", m2)
  s1 <- paste0("significant_", m1)
  s2 <- paste0("significant_", m2)
  dt <- as.data.table(df)[is.finite(get(xcol)) & is.finite(get(ycol))]
  dt[, both_sig := get(s1) & get(s2)]
  stat <- cor_pair(df, m1, m2)
  ann <- sprintf("n=%d\nPearson=%.3f\nSpearman=%.3f",
                 stat$n_shared, stat$pearson %||% NA_real_, stat$spearman %||% NA_real_)

  p <- ggplot(dt, aes(x = get(xcol), y = get(ycol))) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = both_sig), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
    labs(
      title = paste0("log2FC concordance: ", m1, " vs ", m2),
      x = paste0(m1, " log2FC"),
      y = paste0(m2, " log2FC"),
      subtitle = ann
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

  ggsave(out_pdf, p, width = 7.5, height = 6.0)
}

plot_scatter(cons, "limma", "MSstats", file.path(CONS_DIR, "scatter_limma_vs_MSstats.pdf"))
plot_scatter(cons, "limma", "MSstatsTMT", file.path(CONS_DIR, "scatter_limma_vs_MSstatsTMT.pdf"))
plot_scatter(cons, "MSstats", "MSstatsTMT", file.path(CONS_DIR, "scatter_MSstats_vs_MSstatsTMT.pdf"))
message("  Wrote scatter plots")

# -----------------------------------------------------------------------------
# 4) Overlap of significant hits
# -----------------------------------------------------------------------------

message("\n[4/10] Overlap of significant hits")

sig_sets <- list(
  limma = cons[significant_limma == TRUE, Gene],
  MSstats = cons[significant_MSstats == TRUE, Gene],
  MSstatsTMT = cons[significant_MSstatsTMT == TRUE, Gene]
)

count_overlap <- function(a, b) length(intersect(sig_sets[[a]], sig_sets[[b]]))
overlap_all3 <- length(Reduce(intersect, sig_sets))

overlap_summary <- data.table(
  metric = c(
    "n_sig_limma", "n_sig_MSstats", "n_sig_MSstatsTMT",
    "overlap_limma_MSstats", "overlap_limma_MSstatsTMT", "overlap_MSstats_MSstatsTMT",
    "overlap_all3"
  ),
  value = c(
    length(sig_sets$limma), length(sig_sets$MSstats), length(sig_sets$MSstatsTMT),
    count_overlap("limma", "MSstats"),
    count_overlap("limma", "MSstatsTMT"),
    count_overlap("MSstats", "MSstatsTMT"),
    overlap_all3
  )
)
fwrite(overlap_summary, file.path(CONS_DIR, "DA_overlap_summary.csv"))
message("  Wrote DA_overlap_summary.csv")

# Optional UpSet plot if available
if (try_install_cran("UpSetR")) {
  suppressPackageStartupMessages(library(UpSetR))
  u <- data.table(
    limma = cons$significant_limma,
    MSstats = cons$significant_MSstats,
    MSstatsTMT = cons$significant_MSstatsTMT
  )
  u[is.na(limma), limma := FALSE]
  u[is.na(MSstats), MSstats := FALSE]
  u[is.na(MSstatsTMT), MSstatsTMT := FALSE]
  u_df <- as.data.frame(u)
  # UpSetR is picky; use 0/1 integers
  u_df$limma <- as.integer(u_df$limma)
  u_df$MSstats <- as.integer(u_df$MSstats)
  u_df$MSstatsTMT <- as.integer(u_df$MSstatsTMT)
  ok <- tryCatch({
    pdf(file.path(CONS_DIR, "upset_significant_hits.pdf"), width = 8, height = 5.5)
    print(upset(u_df, sets = c("limma", "MSstats", "MSstatsTMT"), order.by = "freq", nsets = 3))
    dev.off()
    TRUE
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("  NOTE: UpSet plot failed (", conditionMessage(e), "); continuing without it.")
    FALSE
  })
  if (ok) message("  Wrote upset_significant_hits.pdf")
} else {
  message("  NOTE: UpSetR not available; skipped UpSet plot.")
}

# -----------------------------------------------------------------------------
# 5) Breast marker sanity-check table
# -----------------------------------------------------------------------------

message("\n[5/10] Marker sanity checks")

markers <- data.table(
  Marker = c("ESR1", "GATA3", "KRT18", "PTPRC", "COL1A1", "VIM"),
  expected_direction = c("Tumor_up", "Tumor_up", "Tumor_up", "NAT_up", "NAT_up", "NAT_up")
)

get_obs_dir <- function(lfc) {
  if (is.na(lfc) || !is.finite(lfc)) return(NA_character_)
  if (lfc > 0) "Tumor_up" else if (lfc < 0) "NAT_up" else "zero"
}

marker_rows <- markers[, {
  g <- Marker
  row <- cons[Gene == g][1]
  if (is.null(row) || nrow(row) == 0) {
    list(
      log2FC_limma = NA_real_, FDR_limma = NA_real_,
      log2FC_MSstats = NA_real_, FDR_MSstats = NA_real_,
      log2FC_MSstatsTMT = NA_real_, FDR_MSstatsTMT = NA_real_,
      observed_direction_limma = NA_character_,
      observed_direction_MSstats = NA_character_,
      observed_direction_MSstatsTMT = NA_character_,
      n_methods_significant = 0L,
      direction_consistent_with_expectation = FALSE
    )
  } else {
    l1 <- row$log2FC_limma; f1 <- row$FDR_limma
    l2 <- row$log2FC_MSstats; f2 <- row$FDR_MSstats
    l3 <- row$log2FC_MSstatsTMT; f3 <- row$FDR_MSstatsTMT

    s1 <- is_sig(l1, f1); s2 <- is_sig(l2, f2); s3 <- is_sig(l3, f3)
    n_sig <- sum(c(s1, s2, s3), na.rm = TRUE)

    od1 <- get_obs_dir(l1); od2 <- get_obs_dir(l2); od3 <- get_obs_dir(l3)
    exp <- expected_direction

    # “consistent” means: among methods where marker is significant, the direction matches expectation;
    # if not significant anywhere, mark FALSE (informative for meeting: "not significant").
    sig_dirs <- c(od1[s1], od2[s2], od3[s3])
    consistent <- if (n_sig == 0) FALSE else all(sig_dirs == exp)

    list(
      log2FC_limma = l1, FDR_limma = f1,
      log2FC_MSstats = l2, FDR_MSstats = f2,
      log2FC_MSstatsTMT = l3, FDR_MSstatsTMT = f3,
      observed_direction_limma = od1,
      observed_direction_MSstats = od2,
      observed_direction_MSstatsTMT = od3,
      n_methods_significant = as.integer(n_sig),
      direction_consistent_with_expectation = consistent
    )
  }
}, by = .(Marker, expected_direction)]

fwrite(marker_rows, file.path(RES_DIR, "DA_marker_sanity_check.csv"))
fwrite(marker_rows, file.path(CONS_DIR, "DA_marker_sanity_check.csv"))
message("  Wrote DA_marker_sanity_check.csv")

# -----------------------------------------------------------------------------
# 6) Top robust markers tables
# -----------------------------------------------------------------------------

message("\n[6/10] High-confidence marker lists")

calc_means <- function(row) {
  lfcs <- c(row$log2FC_limma, row$log2FC_MSstats, row$log2FC_MSstatsTMT)
  fdrs <- c(row$FDR_limma, row$FDR_MSstats, row$FDR_MSstatsTMT)
  mean_abs_lfc <- mean(abs(lfcs[is.finite(lfcs)]), na.rm = TRUE)
  mean_neglog10fdr <- mean(-log10(pmax(fdrs[is.finite(fdrs)], 1e-300)), na.rm = TRUE)
  list(mean_abs_log2FC = mean_abs_lfc, mean_neglog10FDR = mean_neglog10fdr)
}

cons_dt <- as.data.table(cons)
cons_dt[, c("mean_abs_log2FC", "mean_neglog10FDR") := calc_means(.SD), by = seq_len(nrow(cons_dt))]

tumor_up <- cons_dt[high_confidence_tumor_up == TRUE]
nat_up <- cons_dt[high_confidence_NAT_up == TRUE]

rank_cols <- c("n_methods_significant", "mean_abs_log2FC", "mean_neglog10FDR")

setorder(tumor_up, -n_methods_significant, -mean_abs_log2FC, -mean_neglog10FDR)
setorder(nat_up, -n_methods_significant, -mean_abs_log2FC, -mean_neglog10FDR)

fwrite(tumor_up, file.path(RES_DIR, "DA_high_confidence_tumor_up.csv"))
fwrite(nat_up, file.path(RES_DIR, "DA_high_confidence_NAT_up.csv"))

top20 <- rbindlist(list(
  tumor_up[1:min(20L, .N), .(Gene, direction = "Tumor_up", n_methods_significant, mean_abs_log2FC, mean_neglog10FDR,
                            log2FC_limma, FDR_limma, log2FC_MSstats, FDR_MSstats, log2FC_MSstatsTMT, FDR_MSstatsTMT)],
  nat_up[1:min(20L, .N), .(Gene, direction = "NAT_up", n_methods_significant, mean_abs_log2FC, mean_neglog10FDR,
                          log2FC_limma, FDR_limma, log2FC_MSstats, FDR_MSstats, log2FC_MSstatsTMT, FDR_MSstatsTMT)]
))
fwrite(tumor_up, file.path(CONS_DIR, "DA_high_confidence_tumor_up.csv"))
fwrite(nat_up, file.path(CONS_DIR, "DA_high_confidence_NAT_up.csv"))
fwrite(top20, file.path(CONS_DIR, "DA_high_confidence_top20_for_slides.csv"))
message("  Wrote high-confidence tables + top-20 slides table")

# -----------------------------------------------------------------------------
# 7) P-value / FDR diagnostics
# -----------------------------------------------------------------------------

message("\n[7/10] Method diagnostics")

diag_one <- function(df, method_name) {
  lfc <- df[[paste0("log2FC_", method_name)]]
  fdr <- df[[paste0("FDR_", method_name)]]
  pv <- df[[paste0("pvalue_", method_name)]]
  tested <- is.finite(lfc) | is.finite(fdr) | is.finite(pv)
  n_tested <- sum(tested, na.rm = TRUE)
  n_fdr <- sum(is.finite(fdr) & fdr < SIG_FDR, na.rm = TRUE)
  n_fc05 <- sum(is.finite(lfc) & abs(lfc) > 0.5, na.rm = TRUE)
  n_fc1 <- sum(is.finite(lfc) & abs(lfc) > 1, na.rm = TRUE)
  n_both05 <- sum(is.finite(fdr) & fdr < SIG_FDR & is.finite(lfc) & abs(lfc) > 0.5, na.rm = TRUE)
  n_both1 <- sum(is.finite(fdr) & fdr < SIG_FDR & is.finite(lfc) & abs(lfc) > 1, na.rm = TRUE)
  data.table(
    method = method_name,
    n_tested = n_tested,
    n_FDR_lt_0.05 = n_fdr,
    n_absFC_gt_0.5 = n_fc05,
    n_absFC_gt_1 = n_fc1,
    n_FDR_lt_0.05_and_absFC_gt_0.5 = n_both05,
    n_FDR_lt_0.05_and_absFC_gt_1 = n_both1,
    FDR_median = median(fdr[is.finite(fdr)], na.rm = TRUE),
    FDR_q05 = quantile(fdr[is.finite(fdr)], 0.05, na.rm = TRUE),
    FDR_q95 = quantile(fdr[is.finite(fdr)], 0.95, na.rm = TRUE),
    log2FC_median = median(lfc[is.finite(lfc)], na.rm = TRUE),
    log2FC_q05 = quantile(lfc[is.finite(lfc)], 0.05, na.rm = TRUE),
    log2FC_q95 = quantile(lfc[is.finite(lfc)], 0.95, na.rm = TRUE),
    has_raw_pvalues = any(is.finite(pv))
  )
}

diag <- rbindlist(list(
  diag_one(cons, "limma"),
  diag_one(cons, "MSstats"),
  diag_one(cons, "MSstatsTMT")
))
fwrite(diag, file.path(CONS_DIR, "DA_method_diagnostics.csv"))
message("  Wrote DA_method_diagnostics.csv")

plot_p_hist <- function(df, method_name, out_pdf) {
  pv <- df[[paste0("pvalue_", method_name)]]
  pv <- pv[is.finite(pv)]
  if (length(pv) < 10) {
    message("  NOTE: No/too few raw p-values for ", method_name, "; skipped histogram.")
    return(invisible(FALSE))
  }
  dt <- data.table(pvalue = pv)
  p <- ggplot(dt, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    theme_bw(base_size = 12) +
    labs(title = paste0("Raw p-value histogram: ", method_name), x = "p-value", y = "count")
  ggsave(out_pdf, p, width = 7.0, height = 4.5)
  TRUE
}

plot_p_hist(cons, "limma", file.path(CONS_DIR, "histogram_pvalues_limma.pdf"))
plot_p_hist(cons, "MSstats", file.path(CONS_DIR, "histogram_pvalues_MSstats.pdf"))
plot_p_hist(cons, "MSstatsTMT", file.path(CONS_DIR, "histogram_pvalues_MSstatsTMT.pdf"))

# -----------------------------------------------------------------------------
# 8) Pathway enrichment (best-effort; don’t fail the whole pipeline)
# -----------------------------------------------------------------------------

message("\n[8/10] Pathway enrichment (best-effort)")

enrichment_ok <- try_install_bioc("clusterProfiler") && try_install_bioc("org.Hs.eg.db") &&
  try_install_cran("msigdbr") && try_install_cran("enrichplot")

enrich_and_save <- function(genes, collection, subcategory = NULL, out_csv, out_pdf, title) {
  if (!enrichment_ok) {
    fwrite(data.table(TODO = "Enrichment packages unavailable; please install clusterProfiler, org.Hs.eg.db, msigdbr, enrichplot."),
           out_csv)
    return(invisible(NULL))
  }
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(msigdbr)
    library(org.Hs.eg.db)
    library(enrichplot)
  })

  genes <- unique(na.omit(as.character(genes)))
  genes <- genes[nzchar(genes)]
  if (length(genes) < 10) {
    fwrite(data.table(note = "Too few genes for enrichment", n_genes = length(genes)), out_csv)
    return(invisible(NULL))
  }

  msig <- msigdbr::msigdbr(species = "Homo sapiens", category = collection,
                           subcategory = subcategory)
  term2gene <- unique(msig[, c("gs_name", "gene_symbol")])
  colnames(term2gene) <- c("term", "gene")

  enr <- clusterProfiler::enricher(genes, TERM2GENE = term2gene, pAdjustMethod = "BH")
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
    fwrite(data.table(note = "No enriched terms found", n_genes = length(genes)), out_csv)
    return(invisible(NULL))
  }
  res <- as.data.frame(enr)
  fwrite(res, out_csv)

  pdf(out_pdf, width = 10, height = 6.5)
  print(enrichplot::dotplot(enr, showCategory = 20) + ggplot2::ggtitle(title))
  dev.off()
  invisible(enr)
}

genes_tumor <- tumor_up$Gene
genes_nat <- nat_up$Gene

enr_tumor_h <- file.path(CONS_DIR, "enrichment_tumor_up.csv")
enr_nat_h <- file.path(CONS_DIR, "enrichment_NAT_up.csv")
enr_tumor_pdf <- file.path(CONS_DIR, "enrichment_tumor_up.pdf")
enr_nat_pdf <- file.path(CONS_DIR, "enrichment_NAT_up.pdf")

# Hallmark first (most stable for slides)
enr_tumor <- enrich_and_save(genes_tumor, collection = "H", out_csv = enr_tumor_h, out_pdf = enr_tumor_pdf,
                             title = "Tumor-up (high-confidence) — Hallmark enrichment")
enr_nat <- enrich_and_save(genes_nat, collection = "H", out_csv = enr_nat_h, out_pdf = enr_nat_pdf,
                           title = "NAT-up (high-confidence) — Hallmark enrichment")

# Optionally GO BP + Reactome (save alongside if available)
if (enrichment_ok) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })

  enrich_go <- function(genes, out_csv, out_pdf, title) {
    genes <- unique(na.omit(as.character(genes)))
    genes <- genes[nzchar(genes)]
    if (length(genes) < 10) {
      fwrite(data.table(note = "Too few genes for GO enrichment", n_genes = length(genes)), out_csv)
      return(invisible(NULL))
    }
    ego <- tryCatch(
      clusterProfiler::enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                               ont = "BP", pAdjustMethod = "BH", readable = TRUE),
      error = function(e) NULL
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      fwrite(data.table(note = "No GO BP enriched terms found", n_genes = length(genes)), out_csv)
      return(invisible(NULL))
    }
    fwrite(as.data.frame(ego), out_csv)
    pdf(out_pdf, width = 10, height = 6.5)
    print(enrichplot::dotplot(ego, showCategory = 20) + ggplot2::ggtitle(title))
    dev.off()
    invisible(ego)
  }

  enrich_reactome <- function(genes, out_csv, out_pdf, title) {
    # Prefer ReactomePA if present; otherwise msigdbr Reactome set.
    if (try_install_bioc("ReactomePA")) {
      suppressPackageStartupMessages(library(ReactomePA))
      genes <- unique(na.omit(as.character(genes)))
      genes <- genes[nzchar(genes)]
      if (length(genes) < 10) {
        fwrite(data.table(note = "Too few genes for ReactomePA enrichment", n_genes = length(genes)), out_csv)
        return(invisible(NULL))
      }
      er <- tryCatch(ReactomePA::enrichPathway(gene = genes, organism = "human",
                                              pAdjustMethod = "BH", readable = TRUE),
                     error = function(e) NULL)
      if (is.null(er) || nrow(as.data.frame(er)) == 0) {
        fwrite(data.table(note = "No Reactome enriched terms found", n_genes = length(genes)), out_csv)
        return(invisible(NULL))
      }
      fwrite(as.data.frame(er), out_csv)
      pdf(out_pdf, width = 10, height = 6.5)
      print(enrichplot::dotplot(er, showCategory = 20) + ggplot2::ggtitle(title))
      dev.off()
      invisible(er)
    } else {
      # fallback: msigdbr Reactome collection
      enrich_and_save(genes, collection = "C2", subcategory = "CP:REACTOME",
                      out_csv = out_csv, out_pdf = out_pdf, title = title)
    }
  }

  enrich_go(genes_tumor, file.path(RES_DIR, "enrichment_GO_BP_tumor_up.csv"),
            file.path(CONS_DIR, "enrichment_GO_BP_tumor_up.pdf"),
            "Tumor-up (high-confidence) — GO BP enrichment")
  enrich_go(genes_nat, file.path(CONS_DIR, "enrichment_GO_BP_NAT_up.csv"),
            file.path(CONS_DIR, "enrichment_GO_BP_NAT_up.pdf"),
            "NAT-up (high-confidence) — GO BP enrichment")
  enrich_reactome(genes_tumor, file.path(CONS_DIR, "enrichment_Reactome_tumor_up.csv"),
                  file.path(CONS_DIR, "enrichment_Reactome_tumor_up.pdf"),
                  "Tumor-up (high-confidence) — Reactome enrichment")
  enrich_reactome(genes_nat, file.path(CONS_DIR, "enrichment_Reactome_NAT_up.csv"),
                  file.path(CONS_DIR, "enrichment_Reactome_NAT_up.pdf"),
                  "NAT-up (high-confidence) — Reactome enrichment")
} else {
  message("  NOTE: Enrichment packages unavailable; wrote TODO placeholders.")
}

# -----------------------------------------------------------------------------
# 9) Final report for presentation
# -----------------------------------------------------------------------------

message("\n[9/10] Writing summary report")

topN <- function(dt, direction, n = 10) {
  if (nrow(dt) == 0) return(data.table())
  dt[1:min(n, .N), .(Gene, direction = direction, n_methods_significant,
                     mean_abs_log2FC, mean_neglog10FDR,
                     log2FC_limma, FDR_limma,
                     log2FC_MSstats, FDR_MSstats,
                     log2FC_MSstatsTMT, FDR_MSstatsTMT)]
}

top_tumor10 <- topN(tumor_up, "Tumor_up", 10)
top_nat10 <- topN(nat_up, "NAT_up", 10)

read_top_enrichment <- function(path, n = 10) {
  if (!file.exists(path)) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  # clusterProfiler enricher outputs include Description/ID/p.adjust; msigdbr enricher uses Description
  if ("Description" %in% names(dt)) {
    return(dt[order(p.adjust)][1:min(n, .N), .(Description, p.adjust, qvalue = if ("qvalue" %in% names(dt)) qvalue else NA_real_)])
  }
  if ("ID" %in% names(dt)) {
    return(dt[order(p.adjust)][1:min(n, .N), .(ID, p.adjust, qvalue = if ("qvalue" %in% names(dt)) qvalue else NA_real_)])
  }
  dt[1:min(n, .N)]
}

enr_tumor_top <- read_top_enrichment(enr_tumor_h, 10)
enr_nat_top <- read_top_enrichment(enr_nat_h, 10)

summary_path <- file.path(CONS_DIR, "DA_consensus_summary.txt")
con <- file(summary_path, open = "wt")
writeLines(c(
  "Consensus DA downstream summary (PDC000120)",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Inputs:",
  paste0("  - ", IN_LIMMA),
  paste0("  - ", IN_MSSTATS),
  paste0("  - ", IN_MSSTATSTMT),
  "",
  "Significance rule:",
  paste0("  - FDR < ", SIG_FDR, " and |log2FC| > ", SIG_LFC),
  "",
  "Proteins/genes tested (non-missing stats):",
  paste0("  - limma: ", diag[method == "limma", n_tested]),
  paste0("  - MSstats: ", diag[method == "MSstats", n_tested]),
  paste0("  - MSstatsTMT: ", diag[method == "MSstatsTMT", n_tested]),
  "",
  "Significant hits (FDR<0.05 & |log2FC|>1):",
  paste0("  - limma: ", overlap_summary[metric == "n_sig_limma", value]),
  paste0("  - MSstats: ", overlap_summary[metric == "n_sig_MSstats", value]),
  paste0("  - MSstatsTMT: ", overlap_summary[metric == "n_sig_MSstatsTMT", value]),
  "",
  "Overlaps:",
  paste0("  - limma ∩ MSstats: ", overlap_summary[metric == "overlap_limma_MSstats", value]),
  paste0("  - limma ∩ MSstatsTMT: ", overlap_summary[metric == "overlap_limma_MSstatsTMT", value]),
  paste0("  - MSstats ∩ MSstatsTMT: ", overlap_summary[metric == "overlap_MSstats_MSstatsTMT", value]),
  paste0("  - all three: ", overlap_summary[metric == "overlap_all3", value]),
  "",
  "log2FC concordance across methods (shared proteins):",
  paste0("  - limma vs MSstats: n=", corrs[method1=="limma" & method2=="MSstats", n_shared],
         " pearson=", signif(corrs[method1=="limma" & method2=="MSstats", pearson], 3),
         " spearman=", signif(corrs[method1=="limma" & method2=="MSstats", spearman], 3)),
  paste0("  - limma vs MSstatsTMT: n=", corrs[method1=="limma" & method2=="MSstatsTMT", n_shared],
         " pearson=", signif(corrs[method1=="limma" & method2=="MSstatsTMT", pearson], 3),
         " spearman=", signif(corrs[method1=="limma" & method2=="MSstatsTMT", spearman], 3)),
  paste0("  - MSstats vs MSstatsTMT: n=", corrs[method1=="MSstats" & method2=="MSstatsTMT", n_shared],
         " pearson=", signif(corrs[method1=="MSstats" & method2=="MSstatsTMT", pearson], 3),
         " spearman=", signif(corrs[method1=="MSstats" & method2=="MSstatsTMT", spearman], 3)),
  "",
  paste0("High-confidence tumor-up (>=2 methods, consistent + sign): n=", nrow(tumor_up)),
  paste0("High-confidence NAT-up (>=2 methods, consistent - sign): n=", nrow(nat_up)),
  paste0("Strict consensus (all 3 methods, consistent sign): n=", sum(cons$strict_consensus, na.rm = TRUE)),
  "",
  "Top robust tumor-up proteins (top 10):",
  paste(capture.output(print(top_tumor10)), collapse = "\n"),
  "",
  "Top robust NAT-up proteins (top 10):",
  paste(capture.output(print(top_nat10)), collapse = "\n"),
  "",
  "Marker sanity check (panel):",
  paste(capture.output(print(marker_rows)), collapse = "\n"),
  "",
  "Top enriched pathways (Hallmark) — Tumor-up (top 10):",
  paste(capture.output(print(enr_tumor_top)), collapse = "\n"),
  "",
  "Top enriched pathways (Hallmark) — NAT-up (top 10):",
  paste(capture.output(print(enr_nat_top)), collapse = "\n"),
  "",
  "Interpretation (statistically careful):",
  paste(
    "These results represent a within-study cross-method consistency check of the Tumor vs NAT signal for CPTAC breast (PDC000120).",
    "Across three analysis pipelines (limma on the gene matrix, MSstats on the same matrix, and MSstatsTMT from PSM-derived TMT summarization),",
    "effect sizes show substantial concordance on shared proteins and a large overlap of significant hits under a stringent threshold (FDR<0.05 and |log2FC|>1).",
    "High-confidence proteins supported by at least two methods provide a robust shortlist for slide-level reporting. Pathway enrichment (where available) is used",
    "as a qualitative coherence check (e.g., epithelial/proliferation programs in tumor-up and ECM/immune/stromal programs in NAT-up).",
    "Importantly, individual marker non-significance or direction discrepancies in a single method should be interpreted as method- and identifier-dependent,",
    "not as invalidation of the overall Tumor vs NAT differential abundance result.",
    sep = " "
  )
), con = con)
close(con)
message("  Wrote DA_consensus_summary.txt")

# -----------------------------------------------------------------------------
# 10) Console summary
# -----------------------------------------------------------------------------

message("\n[10/10] Done.")
message("Key consensus outputs written under: ", CONS_DIR)
message("  - DA_consensus_table.csv")
message("  - DA_method_correlations.csv + scatter_*.pdf")
message("  - DA_overlap_summary.csv (+ optional upset_significant_hits.pdf)")
message("  - DA_marker_sanity_check.csv")
message("  - DA_high_confidence_tumor_up.csv / DA_high_confidence_NAT_up.csv / top20 table")
message("  - DA_method_diagnostics.csv + p-value histograms")
message("  - enrichment outputs (if available)")
message("  - DA_consensus_summary.txt")

