#!/usr/bin/env Rscript
#
# PDC PSM -> protein matrix (gene symbols) using MSstats TMT
#
# Reads PDC .psm files, converts to MSstatsTMT format, runs protein summarization,
# maps RefSeq IDs to gene symbols, and writes a sample x gene matrix.
#
# Usage (run ONE command at a time; use --no-init-file if you see renv errors):
#
#   Rscript --no-init-file pdc_psm_to_protein_matrix.R
#
#   Rscript --no-init-file pdc_psm_to_protein_matrix.R --psm_dir pdc_psm --out_matrix my_matrix.csv
#
#   Rscript --no-init-file pdc_psm_to_protein_matrix.R --annotation samples_annotation.csv --out_matrix protein_matrix_gene_symbols.csv
#
#   Rscript --no-init-file pdc_psm_to_protein_matrix.R --species Mm
#
# Options:
#   --psm_dir      Directory containing .psm files (recursive). Default: pdc_psm/
#   --annotation   Optional CSV: Run, Channel, Condition, BioReplicate, Mixture (use Condition='Norm' for reference channel)
#   --out_matrix   Output matrix path. Default: protein_matrix_gene_symbols.csv
#   --species      org package for gene mapping: Hs (human), Mm (mouse). Default: Hs
#
# MSstatsTMT input (built from PSM): ProteinName, PeptideSequence, Charge, PSM, Mixture, TechRepMixture,
#   Run, Channel, Condition, BioReplicate, Intensity. Before summarization the script checks these and
#   sets reference_norm=FALSE if no 'Norm' channel is present.

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("MSstatsTMT", quietly = TRUE)) BiocManager::install("MSstatsTMT", update = FALSE, ask = FALSE)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
})

library(MSstatsTMT)
library(data.table)

# ------------------------------------------------------------------------------
# Parse PDC PSM file: Protein (RefSeq), PeptideSequence, Charge, TMT   intensities
# TMT columns are like "TMT10-126", values like "23820.3/0.48" or "0/?"
# ------------------------------------------------------------------------------
# Pick first matching column name (case-insensitive)
pick_col <- function(nms, candidates) {
  low <- tolower(nms)
  for (c in candidates) {
    j <- match(tolower(c), low)
    if (!is.na(j)) return(nms[j])
  }
  NULL
}

parse_pdc_psm_file <- function(path) {
  d <- data.table::fread(path, sep = "\t", header = TRUE, na.strings = c("", "NA", "0/?", "?"))
  nms <- names(d)
  # Strip BOM from first column name if present (rename by index to avoid duplicate-name errors)
  nms[1L] <- sub("^\ufeff", "", nms[1L])
  if (names(d)[1L] != nms[1L]) setnames(d, 1L, nms[1L])
  nms <- names(d)
  # Resolve columns by common names
  filecol <- pick_col(nms, c("FileName", "File Name", "file_name"))
  if (is.null(filecol)) filecol <- nms[1L]
  chargecol <- pick_col(nms, c("OriginalCharge", "QueryCharge", "Charge"))
  proteincol <- pick_col(nms, c("Protein", "ProteinName", "Protein.Accessions"))
  peptidecol <- pick_col(nms, c("PeptideSequence", "Annotated.Sequence", "Sequence"))
  if (is.null(chargecol) || is.null(proteincol) || is.null(peptidecol))
    return(data.table(ProteinName = character(), PeptideSequence = character(), Charge = integer(),
                      PSM = character(), Run = character(), Channel = character(), Intensity = numeric()))
  run_vals <- sub("\\.raw$", "", basename(as.character(d[[filecol]])))
  protein_vals <- sub("\\(pre=.*$", "", trimws(vapply(strsplit(as.character(d[[proteincol]]), ";"), `[`, character(1), 1)))
  peptide_vals <- as.character(d[[peptidecol]])
  charge_vals <- as.integer(d[[chargecol]])
  d[, Run := run_vals]
  d[, ProteinName := protein_vals]
  d[, PeptideSequence := peptide_vals]
  d[, Charge := charge_vals]
  d[, PSM := paste0(PeptideSequence, "_", Charge)]

  # TMT channel columns (e.g. TMT10-126 -> Channel 126)
  tmt_prefix <- "TMT10-"
  tmt_cols <- grep(paste0("^", tmt_prefix), names(d), value = TRUE)
  channel_names <- sub(tmt_prefix, "", tmt_cols)
  if (length(tmt_cols) == 0L) return(d[0L, .(ProteinName, PeptideSequence, Charge, PSM, Run, Channel = character(), Intensity = numeric())])

  # Parse intensity: "value/frac" or "0/?" -> numeric or NA
  parse_intensity <- function(x) {
    if (is.na(x)) return(NA_real_)
    x <- as.character(x)
    if (grepl("^0/?|^\\s*$", x)) return(NA_real_)
    num <- as.numeric(sub("^([0-9.Ee+-]+).*", "\\1", x))
    if (is.na(num) || num <= 0) return(NA_real_)
    num
  }

  # Long format: one row per (PSM, Run, Channel)
  out <- data.table(
    ProteinName = character(),
    PeptideSequence = character(),
    Charge = integer(),
    PSM = character(),
    Run = character(),
    Channel = character(),
    Intensity = numeric()
  )

  for (i in seq_along(tmt_cols)) {
    col <- tmt_cols[i]
    ch <- channel_names[i]
    intensities <- vapply(d[[col]], parse_intensity, numeric(1))
    valid <- !is.na(intensities)
    if (!any(valid)) next
    tmp <- d[valid, .(ProteinName, PeptideSequence, Charge, PSM, Run)]
    tmp[, Channel := ch]
    tmp[, Intensity := intensities[valid]]
    out <- rbind(out, tmp)
  }
  out
}

# ------------------------------------------------------------------------------
# Convert all PSM files in dir to one MSstatsTMT-style long-format data.table
# and build default annotation (Run, Channel -> Condition = Run_Channel, etc.)
# ------------------------------------------------------------------------------
pdc_psm_to_msstats_format <- function(psm_dir, annotation_file = NULL) {
  psm_files <- list.files(path = psm_dir, pattern = "\\.psm$", recursive = TRUE, full.names = TRUE)
  if (length(psm_files) == 0) stop("No .psm files found in ", psm_dir)

  message("Reading ", length(psm_files), " PSM files...")
  all_psm <- rbindlist(lapply(psm_files, parse_pdc_psm_file))

  # Standardize Channel for MSstats (no leading zeros; 126 not 126.0)
  all_psm[, Channel := as.character(Channel)]

  # Build annotation: Run, Channel, Condition, BioReplicate, Mixture, (Fraction, TechRepMixture)
  runs <- unique(all_psm$Run)
  channels <- unique(all_psm$Channel)
  if (is.null(annotation_file) || !file.exists(annotation_file)) {
    annotation <- expand.grid(
      Run = runs,
      Channel = channels,
      stringsAsFactors = FALSE
    )
    annotation$Fraction <- 1L
    annotation$TechRepMixture <- 1L
    annotation$Condition <- paste0(annotation$Run, "_", annotation$Channel)
    annotation$BioReplicate <- annotation$Condition
    annotation$Mixture <- "Mixture1"
  } else {
    annotation <- data.table::fread(annotation_file)
    need <- c("Run", "Channel", "Condition", "BioReplicate", "Mixture")
    if (!all(need %in% names(annotation))) stop("Annotation must have: ", paste(need, collapse = ", "))
    if (!"Fraction" %in% names(annotation)) annotation[, Fraction := 1L]
    if (!"TechRepMixture" %in% names(annotation)) annotation[, TechRepMixture := 1L]
  }

  # Merge annotation into PSM data
  all_psm <- merge(
    all_psm,
    annotation,
    by = c("Run", "Channel"),
    all.x = TRUE
  )
  # Drop rows with no annotation (optional)
  all_psm <- all_psm[!is.na(Condition)]

  # Remove shared peptides (one protein per PSM; already one protein per row)
  # Keep unique peptide-protein pairs: drop if same peptide maps to multiple proteins
  n_per_pep <- all_psm[, .N, by = .(PeptideSequence, Charge)]
  multi_prot <- all_psm[, .(n_prot = uniqueN(ProteinName)), by = .(PeptideSequence, Charge)][n_prot > 1]
  if (nrow(multi_prot) > 0) {
    all_psm <- all_psm[!multi_prot, on = c("PeptideSequence", "Charge")]
  }

  list(psm = all_psm, annotation = annotation)
}

# ------------------------------------------------------------------------------
# RefSeq (NP_ / XP_ ) to gene symbol via org.Hs.eg.db or org.Mm.eg.db
# ------------------------------------------------------------------------------
refseq_to_gene <- function(refseq_ids, species = "Hs") {
  pkg <- switch(species,
    Hs = "org.Hs.eg.db",
    Mm = "org.Mm.eg.db",
    stop("species must be Hs or Mm")
  )
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  suppressPackageStartupMessages(do.call(library, list(pkg, character.only = TRUE)))
  db <- if (species == "Hs") org.Hs.eg.db else org.Mm.eg.db
  ids <- unique(na.omit(refseq_ids))
  ids_no_ver <- sub("\\.[0-9]+$", "", ids)
  out <- tryCatch(
    AnnotationDbi::select(db, keys = unique(ids_no_ver), columns = "SYMBOL", keytype = "REFSEQ"),
    error = function(e) data.frame(REFSEQ = character(), SYMBOL = character())
  )
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame(ProteinName = refseq_ids, GeneSymbol = NA_character_, stringsAsFactors = FALSE))
  }
  out <- out[!is.na(out$SYMBOL) & out$SYMBOL != "", ]
  base_to_symbol <- setNames(out$SYMBOL, sub("\\.[0-9]+$", "", out$REFSEQ))
  full_ids <- unique(refseq_ids)
  full_ids <- full_ids[!is.na(full_ids)]
  data.frame(
    ProteinName = full_ids,
    GeneSymbol = ifelse(
      sub("\\.[0-9]+$", "", full_ids) %in% names(base_to_symbol),
      base_to_symbol[sub("\\.[0-9]+$", "", full_ids)],
      NA_character_
    ),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  psm_dir <- "pdc_psm"
  annotation_file <- NULL
  out_matrix <- "protein_matrix_gene_symbols.csv"
  species <- "Hs"

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--psm_dir" && i < length(args)) {
      psm_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--annotation" && i < length(args)) {
      annotation_file <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--out_matrix" && i < length(args)) {
      out_matrix <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--species" && i < length(args)) {
      species <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  # 1) PDC PSM -> long format + annotation
  conv <- pdc_psm_to_msstats_format(psm_dir, annotation_file)
  input_dt <- conv$psm

  # MSstatsTMT expects columns: ProteinName, PeptideSequence, Charge, PSM, Mixture, TechRepMixture, Run, Channel, BioReplicate, Condition, Intensity
  input_dt <- input_dt[, .(
    ProteinName, PeptideSequence, Charge, PSM,
    Mixture, TechRepMixture, Run, Channel, BioReplicate, Condition, Intensity
  )]

  # Remove rows with NA intensity
  input_dt <- input_dt[!is.na(Intensity) & Intensity > 0]

  if (nrow(input_dt) == 0) stop("No valid intensities after parsing. Check PSM format and TMT columns.")

  message("PSM-level rows for MSstatsTMT: ", nrow(input_dt))

  # ---------- Data check: what MSstatsTMT expects vs what we have ----------
  # MSstatsTMT expects: ProteinName, PeptideSequence, Charge, PSM, Mixture, TechRepMixture,
  #                     Run, Channel, Condition, BioReplicate, Intensity
  required_cols <- c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", "TechRepMixture",
                     "Run", "Channel", "Condition", "BioReplicate", "Intensity")
  missing <- setdiff(required_cols, names(input_dt))
  if (length(missing) > 0) stop("Input missing columns required by MSstatsTMT: ", paste(missing, collapse = ", "))

  n_runs <- uniqueN(input_dt$Run)
  n_channels <- uniqueN(input_dt$Channel)
  n_conditions <- uniqueN(input_dt$Condition)
  n_mixtures <- uniqueN(input_dt$Mixture)
  has_norm <- any(tolower(input_dt$Condition) == "norm")
  reference_norm <- has_norm && n_runs > 1L
  if (!has_norm && n_runs > 1L) message("No 'Norm' channel found; using reference_norm = FALSE (run-to-run normalization skipped).")
  if (has_norm) message("'Norm' channel found; using reference_norm = TRUE for run-to-run normalization.")

  message("Data summary for MSstatsTMT:")
  message("  Runs: ", n_runs, "  Channels: ", n_channels, "  Conditions: ", n_conditions, "  Mixtures: ", n_mixtures)
  message("  Intensity range: [", min(input_dt$Intensity, na.rm = TRUE), ", ", max(input_dt$Intensity, na.rm = TRUE), "]")
  message("  Proteins: ", uniqueN(input_dt$ProteinName), "  Unique PSMs: ", uniqueN(input_dt$PSM))

  # Optional: save a small sample for inspection (first 50k rows)
  sample_path <- "pdc_msstats_input_sample.csv"
  if (nrow(input_dt) > 50000L) {
    fwrite(head(input_dt, 50000L), sample_path)
    message("  Sample of 50,000 rows written to ", sample_path, " for inspection.")
  }

  # 2) Protein summarization (MSstats TMT)
  message("Running MSstatsTMT proteinSummarization...")
  quant <- proteinSummarization(
    input_dt,
    method = "msstats",
    global_norm = TRUE,
    reference_norm = reference_norm,
    remove_norm_channel = FALSE,
    remove_empty_channel = TRUE,
    MBimpute = TRUE
  )

  # 3) Map protein IDs to gene symbols
  message("Mapping RefSeq to gene symbols...")
  prot_col <- if ("Protein" %in% names(quant)) "Protein" else "ProteinName"
  prot_ids <- unique(quant[[prot_col]])
  gene_map <- refseq_to_gene(prot_ids, species = species)
  setnames(gene_map, "ProteinName", prot_col)
  quant <- merge(quant, gene_map, by = prot_col, all.x = TRUE)
  setDT(quant)
  if (!"GeneSymbol" %in% names(quant)) quant[, GeneSymbol := get(prot_col)]
  quant[, GeneSymbol := fifelse(is.na(GeneSymbol) | GeneSymbol == "", get(prot_col), GeneSymbol)]

  # 4) Build sample x gene matrix (one row per Run_Channel, one column per gene)
  # Use Abundance from quant
  sample_id <- paste0(quant$Run, "_", quant$Channel)
  mat_dt <- quant[, .(sample = sample_id, GeneSymbol, Abundance)]
  # If multiple proteins map to same gene in same sample, take median (or mean)
  mat_wide <- dcast(mat_dt, GeneSymbol ~ sample, value.var = "Abundance", fun.aggregate = median, na.rm = TRUE)
  mat <- as.matrix(mat_wide[, -1, with = FALSE])
  rownames(mat) <- mat_wide$GeneSymbol

  # 5) Write matrix
  data.table::fwrite(
    data.table(GeneSymbol = rownames(mat), mat),
    out_matrix,
    sep = ",",
    quote = FALSE
  )
  message("Wrote ", nrow(mat), " genes x ", ncol(mat), " samples to ", out_matrix)
  invisible(list(quant = quant, matrix = mat, path = out_matrix))
}

main()
