#!/usr/bin/env Rscript
# =============================================================================
# Extract per-protein bridge channel summaries from msstats_input.tsv
# =============================================================================
# For each domain (CPTAC PDC000120, CCLE):
#   1. Read msstats_input.tsv, filter Condition == "Norm"
#   2. Aggregate PSMs to protein-level per plex (median log2 intensity)
#   3. Compute cross-plex summaries (median, MAD, mean, SD, IQR, n)
#   4. Map ProteinName to GeneSymbol
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
OUTDIR <- file.path(REPO, "reports/benchmark_master/methods/bridge_aware")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ─── Helper: extract bridge from one msstats_input.tsv ──────────────────
extract_bridge <- function(tsv_path, domain_label, gene_map = NULL) {
  cat("  Reading bridge rows from", basename(tsv_path), "...\n")

  # Stream-read: grep Norm rows then add header
  hdr <- "ProteinName\tPeptideSequence\tCharge\tMixture\tTechRepMixture\tRun\tChannel\tCondition\tBioReplicate\tIntensity\tPSM"
  cmd <- paste0("(echo '", hdr, "'; grep '\tNorm\t' '", tsv_path, "') ")
  bridge_raw <- fread(cmd = cmd, sep = "\t")
  cat("    Bridge PSM rows:", nrow(bridge_raw), "\n")

  # Convert intensity to numeric, take log2
  bridge_raw[, Intensity := as.numeric(Intensity)]
  bridge_raw <- bridge_raw[!is.na(Intensity) & Intensity > 0]
  bridge_raw[, log2int := log2(Intensity)]

  # Aggregate: median log2 intensity per (ProteinName, Mixture)
  plex_summary <- bridge_raw[, .(
    plex_median = median(log2int, na.rm = TRUE),
    plex_n_psm = .N
  ), by = .(ProteinName, Mixture)]
  cat("    Protein × plex combinations:", nrow(plex_summary), "\n")
  cat("    Unique proteins:", uniqueN(plex_summary$ProteinName), "\n")

  # Summarize across plexes
  protein_bridge <- plex_summary[, .(
    bridge_median = median(plex_median, na.rm = TRUE),
    bridge_mad = mad(plex_median, na.rm = TRUE),
    bridge_mean = mean(plex_median, na.rm = TRUE),
    bridge_sd = sd(plex_median, na.rm = TRUE),
    bridge_iqr = IQR(plex_median, na.rm = TRUE),
    bridge_n = .N,
    total_psms = sum(plex_n_psm)
  ), by = ProteinName]

  protein_bridge[, domain := domain_label]

  # Map to gene symbol if mapping provided
  if (!is.null(gene_map)) {
    protein_bridge <- merge(protein_bridge, gene_map, by = "ProteinName", all.x = TRUE)
  }

  cat("    Proteins with bridge summary:", nrow(protein_bridge), "\n")
  cat("    Median bridge_n (plexes per protein):", median(protein_bridge$bridge_n), "\n")
  protein_bridge
}

# ─── Build gene symbol mapping from gene_matrix ─────────────────────────
build_gene_map <- function(gm_path) {
  dt <- fread(gm_path, header = TRUE, nrows = 0)
  if ("UniProtID" %in% names(dt)) {
    dt <- fread(gm_path, header = TRUE, select = c("GeneSymbol", "UniProtID"))
    # Gene matrix first column is GeneSymbol, second may be UniProtID
  } else {
    dt <- fread(gm_path, header = TRUE, select = 1)
  }
  dt
}

# ─── CPTAC studies to extract bridge from ────────────────────────────────
CPTAC_STUDIES <- list(
  PDC000120 = list(
    tissue = "Breast",
    msstats = file.path(REPO, "data/results/PDC000120/msstats_input.tsv"),
    gene_matrix = file.path(REPO, "data/results/PDC000120/gene_matrix.csv")
  ),
  PDC000153 = list(
    tissue = "Lung",
    msstats = file.path(REPO, "data/results/PDC000153/msstats_input.tsv"),
    gene_matrix = file.path(REPO, "data/results/PDC000153/gene_matrix.csv")
  )
)

# ─── Extract bridge from each CPTAC study ────────────────────────────────
cptac_bridge_list <- list()

for (study_id in names(CPTAC_STUDIES)) {
  info <- CPTAC_STUDIES[[study_id]]
  if (!file.exists(info$msstats)) {
    cat("  Skipping", study_id, "— msstats_input.tsv not found\n")
    next
  }
  cat("\n=== CPTAC", study_id, "(", info$tissue, ") Bridge Extraction ===\n")
  bridge_dt <- extract_bridge(info$msstats, paste0("CPTAC_", study_id))

  # Map ProteinName → GeneSymbol using the study's gene_matrix
  gm_map <- fread(info$gene_matrix, header = TRUE, select = c("GeneSymbol", "UniProtID"))
  refseq_gene <- unique(gm_map[, .(ProteinName = UniProtID, GeneSymbol)])
  extra <- refseq_gene[grepl("^XXX_", ProteinName)]
  if (nrow(extra) > 0) {
    extra[, ProteinName := sub("^XXX_", "", ProteinName)]
    refseq_gene <- unique(rbind(refseq_gene, extra))
  }
  cat("  RefSeq → Gene mappings for", study_id, ":", nrow(refseq_gene), "\n")

  bridge_dt <- merge(bridge_dt, refseq_gene, by = "ProteinName", all.x = TRUE)
  bridge_dt[, study_id := study_id]
  cat("  Mapped", sum(!is.na(bridge_dt$GeneSymbol)), "of",
      nrow(bridge_dt), "proteins to gene symbols\n")

  # Save per-study bridge
  fwrite(bridge_dt, file.path(OUTDIR, paste0("bridge_summary_cptac_", study_id, ".tsv")), sep = "\t")
  cptac_bridge_list[[study_id]] <- bridge_dt
}

# Combine all CPTAC bridge summaries. For genes in multiple studies,
# keep the entry with the most plexes (highest bridge_n).
cptac_bridge_all <- rbindlist(cptac_bridge_list, fill = TRUE)
cptac_bridge_all[, domain := "CPTAC"]
cptac_bridge <- cptac_bridge_all[order(-bridge_n)]
cptac_bridge <- cptac_bridge[!duplicated(GeneSymbol) | is.na(GeneSymbol)]
cat("\n  Combined CPTAC bridge:", nrow(cptac_bridge), "proteins from",
    length(cptac_bridge_list), "studies\n")

# ─── Extract CCLE bridge ────────────────────────────────────────────────
cat("\n=== CCLE Bridge Extraction ===\n")
ccle_tsv <- file.path(REPO, "data/results/CCLE_corrected/msstats_input.tsv")
ccle_bridge <- extract_bridge(ccle_tsv, "CCLE")

# CCLE: ProteinName is full UniProt format "sp|ACCESSION|NAME_SPECIES"
gm_ccle_map <- fread(file.path(REPO, "data/results/CCLE_corrected/gene_matrix.csv"),
                      header = TRUE, select = c("GeneSymbol", "UniProtID"))
gm_ccle_map <- unique(gm_ccle_map[, .(ProteinName = trimws(UniProtID),
                                        GeneSymbol = as.character(GeneSymbol))])
ccle_bridge[, ProteinName := trimws(ProteinName)]
ccle_bridge <- merge(ccle_bridge, gm_ccle_map, by = "ProteinName", all.x = TRUE)
cat("  CCLE: mapped", sum(!is.na(ccle_bridge$GeneSymbol) & ccle_bridge$GeneSymbol != ""), "of",
    nrow(ccle_bridge), "proteins to gene symbols (via UniProtID lookup)\n")

# ─── Save bridge summaries ──────────────────────────────────────────────
cat("\n=== Saving Bridge Summaries ===\n")
fwrite(cptac_bridge, file.path(OUTDIR, "bridge_summary_cptac.tsv"), sep = "\t")
fwrite(ccle_bridge, file.path(OUTDIR, "bridge_summary_ccle.tsv"), sep = "\t")
fwrite(cptac_bridge_all, file.path(OUTDIR, "bridge_summary_cptac_all_studies.tsv"), sep = "\t")

# ─── Summary statistics ─────────────────────────────────────────────────
cat("\n=== Summary ===\n")
for (sid in names(cptac_bridge_list)) {
  b <- cptac_bridge_list[[sid]]
  cat("  CPTAC", sid, "bridge proteins:", nrow(b),
      "| with GeneSymbol:", sum(!is.na(b$GeneSymbol)),
      "| median plexes:", median(b$bridge_n), "\n")
}
cat("  Combined CPTAC bridge proteins:", nrow(cptac_bridge), "\n")
cat("  CCLE bridge proteins:", nrow(ccle_bridge), "\n")
cat("  CCLE with GeneSymbol:", sum(!is.na(ccle_bridge$GeneSymbol)), "\n")
cat("  CCLE median plexes per protein:", median(ccle_bridge$bridge_n), "\n")

shared <- intersect(
  cptac_bridge[!is.na(GeneSymbol), GeneSymbol],
  ccle_bridge[!is.na(GeneSymbol), GeneSymbol]
)
cat("  Shared genes with bridge in both:", length(shared), "\n")

cat("\nDone. Outputs in:", OUTDIR, "\n")
