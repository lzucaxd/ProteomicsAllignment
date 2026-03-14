#!/usr/bin/env Rscript
# =============================================================================
# PDC / CPTAC PSM → MSstatsTMT → Gene Matrix (Robust + Paper-Ready)
# =============================================================================
#
# Production pipeline: raw PDC .psm files → MSstatsTMT long-format → protein
# summarization (with bridge/reference channel) → gene symbol matrix.
#
# Requirements: annotation CSV with Run, Channel, Condition, BioReplicate,
# Mixture, Fraction, TechRepMixture. Bridge channel must have Condition = "Norm".
# If annotation is missing, a template is written and the script exits.
#
# Usage:
#   Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R \
#     --psm_dir pdc_psm --annotation annotation.csv --outdir results
#
# Arguments:
#   --psm_dir           Directory containing .psm files (recursive). Required.
#   --annotation        CSV: Run, Channel, Condition, BioReplicate, Mixture, Fraction, TechRepMixture. Optional if --reference_channel set.
#   --reference_channel Bridge channel label (e.g. 131 for TMT10, 126C for TMT11). If set and no annotation, auto-fill Norm/BioReplicate and run.
#   --outdir            Output directory. Default: .
#   --species           Hs (human) or Mm (mouse). Default: Hs
#   --MBimpute          TRUE/FALSE. Default: FALSE (no imputation).
#   --max_runs          Use only first N runs (subsample for fast testing). Default: all.
#   --force_parse       Re-parse PSM even if parsed_psm_long.tsv exists.
#   --force_summarize   Re-run proteinSummarization even if protein_summary.tsv exists.
#   --sample_txt         CPTAC study design: *.sample.txt (FileNameRegEx, AnalyticalSample, channel cols, POOL=bridge). If set, build/audit/rebuild annotation from it.
#   --replace_annotation When rebuilding from sample.txt, overwrite annotation_filled.csv with corrected annotation.
#
# Outputs:
#   parsed_psm_long.tsv   Long-format PSM (one row per PSM x Run x Channel)
#   msstats_input.tsv     Input passed to proteinSummarization
#   protein_summary.tsv   Protein-level abundances from MSstatsTMT
#   gene_matrix.csv       Final sample x gene matrix (Norm channel removed)
#   qc_summary.txt        QC counts and intensity distributions
#   annotation_audit.txt   If --sample_txt: audit of current vs sample.txt (plex, channel, bridge, fraction).
#   annotation_filled_corrected.csv  Corrected annotation from sample.txt (or built when no annotation exists).
#   normalization_audit.txt  If --sample_txt: validation of Norm/POOL; script stops if validation fails.
#
# Fractionated TMT annotation (when using --sample_txt):
#   Mixture = AnalyticalSample (plex id from sample.txt); same across fractions of one plex.
#   BioReplicate = value in sample.txt column(Channel) (actual sample ID); same channel across fractions = same BioReplicate.
#   Condition = "Norm" only when that value is POOL; else "Sample".
#   Run and Fraction preserved from parsed filenames.
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  if (!requireNamespace("MSstatsTMT", quietly = TRUE))
    BiocManager::install("MSstatsTMT", update = FALSE, ask = FALSE)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
})
library(MSstatsTMT)
library(data.table)

# ------------------------------------------------------------------------------
# 1) Column name resolution (avoid hardcoding; PDC/CPTAC column names vary)
# ------------------------------------------------------------------------------
pick_col <- function(nms, candidates) {
  low <- tolower(nms)
  for (c in candidates) {
    j <- match(tolower(c), low)
    if (!is.na(j)) return(nms[j])
  }
  NULL
}

# ------------------------------------------------------------------------------
# 2) Auto-detect TMT format (TMT10-, TMT11-, TMT16-, etc.) and channel names
# ------------------------------------------------------------------------------
detect_tmt_columns <- function(nms) {
  # Match TMTnn- or TMTn- prefix (e.g. TMT10-126, TMT11-131C); exclude TotalAb, FractionOfTotalAb, Flags
  tmt_idx <- grep("^TMT[0-9]+-[0-9]", nms, ignore.case = TRUE)
  if (length(tmt_idx) == 0L) tmt_idx <- grep("^TMT[0-9]+-", nms, ignore.case = TRUE)
  tmt_idx <- tmt_idx[!grepl("TotalAb|FractionOfTotalAb|Flags$", nms[tmt_idx], ignore.case = TRUE)]
  if (length(tmt_idx) == 0L) return(list(prefix = NULL, cols = character(), channels = character()))
  tmt_cols <- nms[tmt_idx]
  prefix <- sub("^(TMT[0-9]+)-.*$", "\\1-", tmt_cols[1L])
  channels <- sub(paste0("^", gsub("-", "\\\\-", prefix)), "", tmt_cols)
  list(prefix = prefix, cols = tmt_cols, channels = channels)
}

# Parse intensity: "value/frac" or "0/?" → numeric; zero or missing → NA
parse_intensity <- function(x) {
  if (is.na(x)) return(NA_real_)
  x <- as.character(x)
  if (grepl("^0/?|^\\s*$", x)) return(NA_real_)
  num <- as.numeric(sub("^([0-9.Ee+-]+).*", "\\1", x))
  if (is.na(num) || num <= 0) return(NA_real_)
  num
}

# ------------------------------------------------------------------------------
# 3) Run / Plex / Fraction from filename (CPTAC/PDC style)
# Biological reasoning: Run = plex (one TMT multiplex); each file = one fraction.
# ------------------------------------------------------------------------------
# Biological reasoning: Run = one MS run (one file); Mixture = plex (group of fractions);
# Fraction = fraction number. Do not treat each fraction as independent biological run.
parse_run_plex_fraction <- function(filename) {
  base <- sub("\\.raw$", "", basename(as.character(filename)))
  # Run = file identifier (one Run per file for MSstatsTMT)
  run_id <- base
  # Mixture = plex/sample group (fractions of same plex share Mixture)
  plex_match <- regmatches(base, regexpr("plex_[0-9]+", base, ignore.case = TRUE))
  if (length(plex_match)) {
    mixture_id <- plex_match[1L]
  } else {
    # CPTAC: 05CPTAC_..._BL_f24 → Mixture = 05CPTAC_..._BL (strip _f24)
    if (length(regmatches(base, regexpr("_f[0-9]+$", base, ignore.case = TRUE))))
      mixture_id <- sub("_f[0-9]+$", "", base, ignore.case = TRUE)
    else
      mixture_id <- base
  }
  # Fraction number: f12, f24, _f12, etc.
  frac_match <- regmatches(base, regexpr("_?f([0-9]+)(_|$)", base, ignore.case = TRUE))
  frac_num <- if (length(frac_match)) as.integer(sub("_?f([0-9]+).*", "\\1", frac_match[1L])) else 1L
  list(Run = run_id, Mixture = mixture_id, Fraction = frac_num, TechRepMixture = 1L)
}

# ------------------------------------------------------------------------------
# 4) Protein ID normalization: RefSeq (NP_, XP_), UniProt (sp|P12345|), strip version, first of ";"
# ------------------------------------------------------------------------------
normalize_protein_id <- function(x) {
  x <- as.character(x)
  # First accession if multiple (semicolon-separated)
  first <- trimws(vapply(strsplit(x, ";"), `[`, character(1), 1))
  # Strip (pre=...,post=...) if present
  first <- sub("\\(pre=.*$", "", first)
  # UniProt: sp|P12345|HUMAN → P12345
  first <- sub("^[a-z]*\\|([A-Z0-9]+)\\|.*$", "\\1", first, ignore.case = TRUE)
  # Strip version .1, .2
  first <- sub("\\.[0-9]+$", "", first)
  first
}

# ------------------------------------------------------------------------------
# 5) Parse a single PDC/CPTAC PSM file → long format (one row per PSM x Channel)
# Uses auto-detected TMT columns; preserves Fraction, Run, Mixture.
# ------------------------------------------------------------------------------
parse_one_psm_file <- function(path) {
  d <- fread(path, sep = "\t", header = TRUE, fill = TRUE, na.strings = c("", "NA", "0/?", "?"))
  nms <- trimws(names(d))
  nms[1L] <- sub("^\ufeff", "", nms[1L])
  if (any(duplicated(nms))) nms <- make.unique(nms, sep = "_")
  setnames(d, seq_len(ncol(d)), nms)  # rename by index to avoid duplicate 'old' names
  nms <- names(d)

  filecol <- pick_col(nms, c("FileName", "File Name", "file_name"))
  if (is.null(filecol)) filecol <- nms[1L]
  chargecol <- pick_col(nms, c("OriginalCharge", "QueryCharge", "Charge"))
  proteincol <- pick_col(nms, c("Protein", "ProteinName", "Protein.Accessions"))
  peptidecol <- pick_col(nms, c("PeptideSequence", "Annotated.Sequence", "Sequence"))
  if (is.null(chargecol) || is.null(proteincol) || is.null(peptidecol))
    return(NULL)

  tmt <- detect_tmt_columns(nms)
  if (length(tmt$cols) == 0L) return(NULL)

  run_meta <- parse_run_plex_fraction(d[[filecol]][1L])
  protein_vals <- normalize_protein_id(d[[proteincol]])
  peptide_vals <- as.character(d[[peptidecol]])
  charge_vals <- as.integer(d[[chargecol]])

  # Build long format without rbind loop: one data.table per channel, then rbindlist once
  out_list <- vector("list", length(tmt$cols))
  for (i in seq_along(tmt$cols)) {
    intens <- vapply(d[[tmt$cols[i]]], parse_intensity, numeric(1))
    valid <- !is.na(intens) & intens > 0
    if (!any(valid)) next
    out_list[[i]] <- data.table(
      ProteinName = protein_vals[valid],
      PeptideSequence = peptide_vals[valid],
      Charge = charge_vals[valid],
      PSM = paste0(peptide_vals[valid], "_", charge_vals[valid]),
      Run = run_meta$Run,
      Mixture = run_meta$Mixture,
      Fraction = run_meta$Fraction,
      TechRepMixture = run_meta$TechRepMixture,
      Channel = tmt$channels[i],
      Intensity = intens[valid]
    )
  }
  out_list <- out_list[!sapply(out_list, is.null)]
  if (length(out_list) == 0L) return(NULL)
  rbindlist(out_list)
}

# ------------------------------------------------------------------------------
# 6) Load all PSM files (data.table; list of DTs then rbindlist)
# ------------------------------------------------------------------------------
load_all_psm <- function(psm_dir) {
  files <- list.files(psm_dir, pattern = "\\.psm$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0L) stop("No .psm files found in ", psm_dir)
  message("Reading ", length(files), " PSM files...")
  out <- rbindlist(lapply(files, parse_one_psm_file), use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0L) stop("No valid PSM rows after parsing. Check TMT columns and format.")
  out
}

# ------------------------------------------------------------------------------
# 7) Annotation: required columns; bridge = Norm. If missing, build template;
#    if reference_channel given, auto-fill (Norm for bridge, BioReplicate=Run_Channel) and continue.
# ------------------------------------------------------------------------------
require_annotation <- function(annotation_path, parsed_psm, outdir, reference_channel = NULL) {
  if (!is.null(annotation_path) && file.exists(annotation_path)) {
    ann <- fread(annotation_path)
    need <- c("Run", "Channel", "Condition", "BioReplicate", "Mixture", "Fraction", "TechRepMixture")
    missing <- setdiff(need, names(ann))
    if (length(missing) > 0L)
      stop("Annotation must contain: ", paste(need, collapse = ", "), ". Missing: ", paste(missing, collapse = ", "))
    return(ann)
  }
  # Build template: one row per (Run, Channel, Fraction) so merge matches all parsed rows
  u <- unique(parsed_psm[, .(Run, Channel, Fraction)])
  template <- copy(u)
  template$Mixture <- template$Run
  template$TechRepMixture <- 1L
  # Auto-fill if reference_channel provided (e.g. "131" or "126C")
  if (!is.null(reference_channel) && nzchar(trimws(reference_channel))) {
    ref <- trimws(as.character(reference_channel))
    template[, Condition := fifelse(Channel == ref, "Norm", "Sample")]
    template[, BioReplicate := paste0(Run, "_", Channel)]
    template <- template[, c("Run", "Channel", "Condition", "BioReplicate", "Mixture", "Fraction", "TechRepMixture")]
    filled_path <- file.path(outdir, "annotation_filled.csv")
    fwrite(template, filled_path)
    message("Auto-filled annotation (reference channel ", ref, " = Norm) -> ", filled_path)
    return(template)
  }
  template$Condition <- ""
  template$BioReplicate <- ""
  template <- template[, c("Run", "Channel", "Condition", "BioReplicate", "Mixture", "Fraction", "TechRepMixture")]
  template_path <- file.path(outdir, "annotation_template.csv")
  fwrite(template, template_path)
  stop(
    "Annotation file is required. A template was written to: ", template_path, "\n",
    "Either fill it (Condition='Norm' for bridge channel, BioReplicate=sample ID) and re-run with --annotation ", template_path,
    " or re-run with --reference_channel 131 (or your bridge channel) to auto-fill and run."
  )
}

# ------------------------------------------------------------------------------
# 7b) Load CPTAC study design (*.sample.txt): FileNameRegEx, AnalyticalSample, channel columns, POOL = bridge
# ------------------------------------------------------------------------------
# Returns list(design_dt, channel_cols, bridge_channel_by_row).
# design_dt: one row per plex; columns FileNameRegEx, AnalyticalSample, and channel names (126, 127N, ...).
# bridge_channel_by_row: character vector, same length as nrow(design_dt), channel label where value == "POOL".
# Detect TMT channel columns dynamically: exclude reserved headers; include columns matching TMT channel pattern (e.g. 126, 127N, 131C).
load_sample_txt <- function(path) {
  d <- fread(path, sep = "\t", header = TRUE, fill = TRUE)
  nms <- trimws(names(d))
  nms[1L] <- sub("^\ufeff", "", nms[1L])
  setnames(d, seq_len(ncol(d)), nms)
  nms <- names(d)
  filecol <- pick_col(nms, c("FileNameRegEx", "FileNameRegex", "File Name RegEx"))
  analcol <- pick_col(nms, c("AnalyticalSample", "Analytical Sample"))
  if (is.null(filecol) || is.null(analcol))
    stop("sample.txt must contain FileNameRegEx and AnalyticalSample. Found: ", paste(nms, collapse = ", "))
  reserved <- tolower(c("FileNameRegEx", "FileNameRegex", "AnalyticalSample", "LabelReagent", "Ratios", filecol, analcol))
  # TMT channel pattern: 3 digits optionally followed by C or N (e.g. 126, 127N, 127C, 131, 131C, 134N)
  is_channel <- function(x) {
    if (tolower(x) %in% reserved) return(FALSE)
    grepl("^[0-9]{3}[CN]?$", x, ignore.case = TRUE)
  }
  channel_cols <- nms[vapply(nms, is_channel, logical(1))]
  if (length(channel_cols) == 0L)
    stop("sample.txt must contain at least one TMT channel column (e.g. 126, 127N, ... 131). Found: ", paste(nms, collapse = ", "))
  # Identify bridge: column whose value is "POOL" (case-insensitive)
  bridge_channel_by_row <- character(nrow(d))
  for (i in seq_len(nrow(d))) {
    for (ch in channel_cols) {
      val <- trimws(as.character(d[[ch]][i]))
      if (toupper(val) == "POOL") {
        bridge_channel_by_row[i] <- ch
        break
      }
    }
  }
  list(design_dt = d, channel_cols = channel_cols, bridge_channel_by_row = bridge_channel_by_row,
       filecol = filecol, analcol = analcol)
}

# ------------------------------------------------------------------------------
# 7c) Match Run to plex: which FileNameRegEx matches this Run?
# ------------------------------------------------------------------------------
# Returns integer vector (plex row index per run); NA if no match or multiple matches.
match_runs_to_plex <- function(runs, design_dt, filecol) {
  out <- rep(NA_integer_, length(runs))
  for (i in seq_along(runs)) {
    run <- runs[i]
    matches <- which(vapply(design_dt[[filecol]], function(pat) {
      grepl(pat, run, ignore.case = TRUE)
    }, logical(1)))
    if (length(matches) == 1L) out[i] <- matches
    # if 0 or >1, leave NA
  }
  out
}

# ------------------------------------------------------------------------------
# 7d) Audit current annotation vs sample.txt design; return list(audit_ok, lines, run_to_plex, expected_bio)
# ------------------------------------------------------------------------------
audit_annotation <- function(annotation_current, sample_design, parsed_psm, outdir) {
  design_dt <- sample_design$design_dt
  channel_cols <- sample_design$channel_cols
  bridge_by_row <- sample_design$bridge_channel_by_row
  filecol <- sample_design$filecol
  analcol <- sample_design$analcol
  lines <- c("=== Annotation audit (current vs CPTAC sample.txt) ===", "")

  runs <- unique(as.character(parsed_psm$Run))
  run_to_plex <- match_runs_to_plex(runs, design_dt, filecol)
  names(run_to_plex) <- runs

  # Build expected: for each (Run, Channel), expected Mixture and BioReplicate from sample.txt
  expected_bio <- list()
  for (run in runs) {
    plex_idx <- run_to_plex[run]
    if (is.na(plex_idx)) next
    for (ch in channel_cols) {
      biorep_exp <- trimws(as.character(design_dt[[ch]][plex_idx]))
      mixture_exp <- as.character(design_dt[[analcol]][plex_idx])
      expected_bio[[length(expected_bio) + 1L]] <- data.table(Run = run, Channel = ch, BioReplicate_expected = biorep_exp, Mixture_expected = mixture_exp)
    }
  }
  expected_dt <- if (length(expected_bio) > 0L) rbindlist(expected_bio) else data.table(Run = character(), Channel = character(), BioReplicate_expected = character(), Mixture_expected = character())

  audit_ok <- TRUE

  # --- 1) Runs matching zero plex rows ---
  no_match <- runs[is.na(run_to_plex[runs])]
  if (length(no_match) > 0L) {
    lines <- c(lines, "1) RUNS MATCHING ZERO PLEX ROWS (FileNameRegEx):", paste("   Count:", length(no_match)), paste("   Runs:", paste(head(no_match, 20), collapse = ", "), if (length(no_match) > 20) "..." else ""), "")
    audit_ok <- FALSE
  } else {
    lines <- c(lines, "1) RUNS MATCHING ZERO PLEX ROWS: None (all runs match exactly one plex).", "")
  }

  # --- 2) Runs matching multiple plex rows ---
  multi_match <- runs[vapply(runs, function(r) {
    if (!is.na(run_to_plex[r])) return(FALSE)
    sum(vapply(design_dt[[filecol]], function(pat) grepl(pat, r, ignore.case = TRUE), logical(1))) > 1L
  }, logical(1))]
  if (length(multi_match) > 0L) {
    lines <- c(lines, "2) RUNS MATCHING MULTIPLE PLEX ROWS:", paste("   Count:", length(multi_match)), paste("   Runs:", paste(head(multi_match, 20), collapse = ", "), if (length(multi_match) > 20) "..." else ""), "")
    audit_ok <- FALSE
  } else {
    lines <- c(lines, "2) RUNS MATCHING MULTIPLE PLEX ROWS: None.", "")
  }

  ann <- as.data.table(annotation_current)
  setkey(ann, Run, Channel)
  setkey(expected_dt, Run, Channel)
  merged_audit <- merge(ann[, .(Run, Channel, Condition, BioReplicate, Mixture, Fraction)], expected_dt, by = c("Run", "Channel"), all = TRUE)

  # --- 3) Current Mixture equals AnalyticalSample? ---
  mixture_mismatch <- merged_audit[!is.na(Mixture_expected) & as.character(Mixture) != as.character(Mixture_expected)]
  if (nrow(mixture_mismatch) > 0L) {
    lines <- c(lines, "3) CURRENT MIXTURE EQUALS ANALYTICALSAMPLE (sample.txt): No.", paste("   Mismatch rows:", nrow(mixture_mismatch)), "")
    audit_ok <- FALSE
  } else {
    lines <- c(lines, "3) CURRENT MIXTURE EQUALS ANALYTICALSAMPLE: Yes (or N/A).", "")
  }

  # --- 4) Current BioReplicate equals channel-mapped sample ID? ---
  bio_mismatch <- merged_audit[!is.na(BioReplicate_expected) & as.character(BioReplicate) != as.character(BioReplicate_expected)]
  if (nrow(bio_mismatch) > 0L) {
    lines <- c(lines, "4) CURRENT BIOREPLICATE EQUALS CHANNEL-MAPPED SAMPLE ID (sample.txt): No.", paste("   Mismatch rows:", nrow(bio_mismatch)), "")
    audit_ok <- FALSE
  } else {
    lines <- c(lines, "4) CURRENT BIOREPLICATE EQUALS CHANNEL-MAPPED SAMPLE ID: Yes (or N/A).", "")
  }

  # --- 5) Current Condition matches POOL (Norm only for POOL)? ---
  norm_rows <- merged_audit[tolower(Condition) == "norm"]
  pool_rows <- merged_audit[toupper(as.character(BioReplicate_expected)) == "POOL"]
  norm_ok <- (nrow(norm_rows) == 0L && nrow(pool_rows) == 0L) || (nrow(norm_rows) > 0L && nrow(pool_rows) > 0L && nrow(merge(norm_rows[, .(Run, Channel)], pool_rows[, .(Run, Channel)], by = c("Run", "Channel"))) == nrow(norm_rows))
  if (!norm_ok) {
    lines <- c(lines, "5) CURRENT CONDITION MATCHES POOL (Norm only where sample.txt = POOL): No.", "")
    audit_ok <- FALSE
  } else {
    lines <- c(lines, "5) CURRENT CONDITION MATCHES POOL: Yes.", "")
  }
  norm_per_plex <- norm_rows[, .(n_norm = .N), by = Mixture][n_norm > 1L]
  if (nrow(norm_per_plex) > 0L) {
    lines <- c(lines, "   FAIL: Multiple Norm channels in same plex (Mixture):", paste(norm_per_plex$Mixture, collapse = ", "), "")
    audit_ok <- FALSE
  }

  # --- 6) (Mixture, Channel) stable across fractions (same BioReplicate)? ---
  frac_check <- merged_audit[!is.na(Mixture_expected), .(n_bio = uniqueN(BioReplicate)), by = .(Mixture_expected, Channel)][n_bio > 1L]
  if (nrow(frac_check) > 0L) {
    lines <- c(lines, "6) (MIXTURE, CHANNEL) STABLE ACROSS FRACTIONS (same BioReplicate): No.", paste("   (Mixture, Channel) pairs with differing BioReplicate:", nrow(frac_check)), "")
    audit_ok <- FALSE
  } else {
    lines <- c(lines, "6) (MIXTURE, CHANNEL) STABLE ACROSS FRACTIONS: Yes (or N/A).", "")
  }

  # --- 7) Appears to use incorrect old logic (Mixture=Run, BioReplicate=Run_Channel)? ---
  mixture_equals_run <- merged_audit[!is.na(Mixture) & as.character(Mixture) == as.character(Run)]
  run_channel_pattern <- merged_audit[, paste0(as.character(Run), "_", as.character(Channel))]
  biorep_equals_run_channel <- merged_audit[!is.na(BioReplicate) & as.character(BioReplicate) == run_channel_pattern]
  old_logic_mixture <- nrow(mixture_equals_run) > 0L
  old_logic_biorep <- nrow(biorep_equals_run_channel) > 0L
  if (old_logic_mixture || old_logic_biorep) {
    lines <- c(lines, "7) APPEARS TO USE INCORRECT OLD LOGIC (Mixture=Run and/or BioReplicate=Run_Channel): Yes.",
      paste("   Mixture == Run rows:", nrow(mixture_equals_run)),
      paste("   BioReplicate == Run_Channel rows:", nrow(biorep_equals_run_channel)), "")
    audit_ok <- FALSE
  } else {
    lines <- c(lines, "7) APPEARS TO USE INCORRECT OLD LOGIC: No.", "")
  }

  n_plex <- nrow(design_dt)
  n_chan <- length(channel_cols)
  expected_samples <- n_plex * (n_chan - 1L)
  actual_samples <- uniqueN(ann[tolower(Condition) != "norm", BioReplicate])
  lines <- c(lines, "SAMPLE COUNT:", paste("   Expected non-POOL (from sample.txt):", expected_samples), paste("   Actual unique BioReplicate (non-Norm):", actual_samples), "")
  if (actual_samples != expected_samples) lines <- c(lines, "   (Mismatch may be acceptable if fractions differ per plex.)", "")

  lines <- c(lines, if (audit_ok) "Result: PASS (annotation matches sample.txt)." else "Result: FAIL (annotation will be rebuilt from sample.txt).")
  audit_path <- file.path(outdir, "annotation_audit.txt")
  writeLines(lines, audit_path)
  message("Annotation audit written to ", audit_path)

  list(audit_ok = audit_ok, audit_lines = lines, run_to_plex = run_to_plex, expected_dt = expected_dt,
       sample_design = sample_design, design_dt = design_dt)
}

# ------------------------------------------------------------------------------
# 7e) Rebuild annotation from sample.txt and parsed PSM (unique Run, Channel, Fraction)
# ------------------------------------------------------------------------------
rebuild_annotation_from_sample_txt <- function(parsed_psm, sample_design) {
  design_dt <- sample_design$design_dt
  channel_cols <- sample_design$channel_cols
  filecol <- sample_design$filecol
  analcol <- sample_design$analcol

  u <- unique(parsed_psm[, .(Run, Channel, Fraction)])
  runs <- unique(u$Run)
  run_to_plex <- match_runs_to_plex(runs, design_dt, filecol)
  names(run_to_plex) <- runs

  out_list <- vector("list", nrow(u))
  for (i in seq_len(nrow(u))) {
    run <- as.character(u$Run[i])
    ch <- as.character(u$Channel[i])
    frac <- as.integer(u$Fraction[i])
    plex_idx <- run_to_plex[run]
    if (is.na(plex_idx) || !(ch %in% channel_cols)) next
    biorep <- trimws(as.character(design_dt[[ch]][plex_idx]))
    mixture <- as.character(design_dt[[analcol]][plex_idx])
    cond <- if (toupper(biorep) == "POOL") "Norm" else "Sample"
    out_list[[i]] <- data.table(Run = run, Channel = ch, Condition = cond, BioReplicate = biorep, Mixture = mixture, Fraction = frac, TechRepMixture = 1L)
  }
  out <- rbindlist(out_list)
  out <- out[!is.na(Run)]
  setkey(out, Run, Channel, Fraction)
  out <- unique(out)
  unmatched_runs <- setdiff(unique(as.character(parsed_psm$Run)), unique(as.character(out$Run)))
  if (length(unmatched_runs) > 0L)
    message("Warning: ", length(unmatched_runs), " run(s) did not match any FileNameRegEx in sample.txt and will be excluded from annotation: ", paste(head(unmatched_runs, 5), collapse = ", "), if (length(unmatched_runs) > 5) " ..." else "")
  out
}

# ------------------------------------------------------------------------------
# 7f) Normalization validation (before MSstatsTMT): one Norm per plex, Norm == POOL, Ratios consistent. Stops on failure.
# ------------------------------------------------------------------------------
validate_normalization_and_audit <- function(annotation_dt, sample_design, outdir) {
  design_dt <- sample_design$design_dt
  bridge_by_row <- sample_design$bridge_channel_by_row
  analcol <- sample_design$analcol
  nms <- names(design_dt)
  ratioscol <- pick_col(nms, c("Ratios", "Ratio"))
  lines <- c("=== Normalization audit (MSstatsTMT reference channel) ===", "")
  validation_ok <- TRUE

  n_plexes <- uniqueN(annotation_dt$Mixture)
  n_fractions_per_plex <- annotation_dt[, .(n_frac = uniqueN(paste(Run, Fraction))), by = Mixture]$n_frac
  n_samples_per_plex <- annotation_dt[tolower(Condition) != "norm", .(n_samp = uniqueN(BioReplicate)), by = Mixture]$n_samp
  lines <- c(lines,
    "COUNTS:",
    paste("   Number of plexes (Mixture):", n_plexes),
    paste("   Fractions per plex (min-max):", paste(min(n_fractions_per_plex), max(n_fractions_per_plex), sep = "-")),
    paste("   Samples per plex (non-Norm, min-max):", paste(min(n_samples_per_plex), max(n_samples_per_plex), sep = "-")),
    "")

  # Fractionated TMT: one *unique* Norm channel label per plex (same channel repeated per fraction).
  # Do not require one Norm *row* per plex — we have one Norm row per (Mixture, Fraction).
  norm_per_plex <- annotation_dt[tolower(Condition) == "norm", .(
    n_unique_norm_channels = uniqueN(Channel),
    norm_channel = paste(unique(Channel), collapse = ","),
    n_norm_rows = .N
  ), by = Mixture]
  if (any(norm_per_plex$n_unique_norm_channels != 1L)) {
    lines <- c(lines, "1) ONE UNIQUE NORM CHANNEL LABEL PER PLEX: FAIL.", paste("   Plexes with != 1 unique Norm channel:", paste(norm_per_plex[n_unique_norm_channels != 1L]$Mixture, collapse = ", ")), "")
    validation_ok <- FALSE
  } else {
    lines <- c(lines, "1) ONE UNIQUE NORM CHANNEL LABEL PER PLEX: PASS.", "")
  }
  lines <- c(lines, "   (Norm channel is repeated across fractions — one row per fraction; row count per plex = number of fractions.)", "")
  for (i in seq_len(nrow(norm_per_plex))) {
    lines <- c(lines, paste("   ", norm_per_plex$Mixture[i], ": unique Norm channel = ", norm_per_plex$norm_channel[i], ", Norm rows (fractions) = ", norm_per_plex$n_norm_rows[i], sep = ""))
  }
  lines <- c(lines, "")

  norm_rows <- annotation_dt[tolower(Condition) == "norm"]
  if (nrow(norm_rows) > 0L) {
    lines <- c(lines, "2) NORM CHANNEL MAPS TO POOL (BioReplicate = POOL) AND IS STABLE ACROSS FRACTIONS:", "")
    all_pool <- all(toupper(as.character(norm_rows$BioReplicate)) == "POOL")
    if (!all_pool) {
      lines <- c(lines, "   FAIL: Some Norm rows have BioReplicate != POOL.", "")
      validation_ok <- FALSE
    } else {
      lines <- c(lines, "   PASS: All Norm rows have BioReplicate = POOL (stable across fractions).", "")
    }
    lines <- c(lines, paste("   (Checked", nrow(norm_rows), "Norm rows across all fractions.)"), "")
  }
  non_norm <- annotation_dt[tolower(Condition) != "norm"]
  non_pool_are_sample <- nrow(non_norm) == 0L || all(tolower(non_norm$Condition) == "sample")
  if (!non_pool_are_sample) {
    lines <- c(lines, "3) ALL NON-POOL CHANNELS ARE SAMPLE: FAIL (some non-Norm rows have Condition != Sample).", "")
    validation_ok <- FALSE
  } else {
    lines <- c(lines, "3) ALL NON-POOL CHANNELS ARE SAMPLE: PASS.", "")
  }
  if (!is.null(ratioscol)) {
    ratios_vals <- trimws(unique(as.character(design_dt[[ratioscol]])))
    lines <- c(lines, "4) RATIOS COLUMN (consistent with bridge channel):", paste("   Column:", ratioscol), paste("   Values in design:", paste(ratios_vals, collapse = ", ")), "   (Ratios should indicate normalization against the bridge/POOL channel.)", "")
  } else {
    lines <- c(lines, "4) RATIOS COLUMN: Not present in sample.txt.", "")
  }
  if (!validation_ok) {
    lines <- c(lines, "Normalization design is invalid. Fix annotation (e.g. use annotation_filled_corrected.csv from sample.txt) and re-run.")
    n_audit_path <- file.path(outdir, "normalization_audit.txt")
    writeLines(lines, n_audit_path)
    stop("Normalization validation failed. See ", n_audit_path)
  }
  lines <- c(lines, "Reference normalization will use POOL channel as bridge. MSstatsTMT reference_norm=TRUE is valid.", "")
  n_audit_path <- file.path(outdir, "normalization_audit.txt")
  writeLines(lines, n_audit_path)
  message("Normalization audit written to ", n_audit_path)
}

# ------------------------------------------------------------------------------
# 8) Merge annotation into PSM long; drop rows without matching annotation
# ------------------------------------------------------------------------------
merge_annotation <- function(parsed_psm, annotation) {
  ann <- copy(annotation)
  setDT(ann)
  # Drop annotation columns from parsed_psm so merge doesn't duplicate (Condition, BioReplicate, TechRepMixture, Mixture)
  drop_from_psm <- c("Condition", "BioReplicate", "TechRepMixture", "Mixture")
  drop_from_psm <- intersect(drop_from_psm, names(parsed_psm))
  if (length(drop_from_psm)) parsed_psm <- parsed_psm[, setdiff(names(parsed_psm), drop_from_psm), with = FALSE]
  cols_key <- c("Run", "Channel", "Fraction")
  out <- merge(parsed_psm, ann, by = cols_key, all.x = FALSE)
  if (nrow(out) == 0L) stop("No rows after merging with annotation. Check Run/Channel/Fraction match.")
  out
}

# ------------------------------------------------------------------------------
# 9) Remove shared peptides (peptide mapped to >1 protein)
# ------------------------------------------------------------------------------
remove_shared_peptides <- function(d) {
  multi <- d[, .(n = uniqueN(ProteinName)), by = .(PeptideSequence, Charge)][n > 1L]
  if (nrow(multi) == 0L) return(d)
  d[!multi, on = c("PeptideSequence", "Charge")]
}

# ------------------------------------------------------------------------------
# 10) Gene symbol mapping (RefSeq / UniProt → gene symbol via org.Hs.eg.db / org.Mm.eg.db)
# ------------------------------------------------------------------------------
protein_to_gene <- function(protein_ids, species = "Hs") {
  pkg <- if (species == "Hs") "org.Hs.eg.db" else "org.Mm.eg.db"
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
  suppressPackageStartupMessages(do.call(library, list(pkg, character.only = TRUE)))
  db <- if (species == "Hs") org.Hs.eg.db else org.Mm.eg.db

  ids <- unique(na.omit(protein_ids))
  # Try with and without version (RefSeq: NP_123.1 → NP_123)
  ids_flat <- unique(c(ids, sub("\\.[0-9]+$", "", ids)))
  # RefSeq (NP_, XP_): try REFSEQ keytype
  out_refseq <- tryCatch(
    AnnotationDbi::select(db, keys = ids_flat, columns = "SYMBOL", keytype = "REFSEQ"),
    error = function(e) data.frame(REFSEQ = character(), SYMBOL = character())
  )
  # UniProt (e.g. P12345): try UNIPROT keytype if available
  keytypes <- AnnotationDbi::keytypes(db)
  out_uniprot <- data.frame(UNIPROT = character(), SYMBOL = character())
  if ("UNIPROT" %in% keytypes) {
    out_uniprot <- tryCatch(
      AnnotationDbi::select(db, keys = ids_flat, columns = "SYMBOL", keytype = "UNIPROT"),
      error = function(e) out_uniprot
    )
  }
  # ALIAS often has RefSeq/UniProt
  out_alias <- data.frame(ALIAS = character(), SYMBOL = character())
  if ("ALIAS" %in% keytypes) {
    out_alias <- tryCatch(
      AnnotationDbi::select(db, keys = ids_flat, columns = "SYMBOL", keytype = "ALIAS"),
      error = function(e) out_alias
    )
  }

  symbol_by_id <- character(length(ids))
  names(symbol_by_id) <- ids
  if (nrow(out_refseq) > 0L && "SYMBOL" %in% names(out_refseq))
    for (i in seq_len(nrow(out_refseq)))
      if (!is.na(out_refseq$SYMBOL[i]) && nzchar(out_refseq$SYMBOL[i]))
        symbol_by_id[out_refseq$REFSEQ[i]] <- out_refseq$SYMBOL[i]
  if (nrow(out_uniprot) > 0L && "SYMBOL" %in% names(out_uniprot))
    for (i in seq_len(nrow(out_uniprot)))
      if (!is.na(out_uniprot$SYMBOL[i]) && nzchar(out_uniprot$SYMBOL[i]))
        symbol_by_id[out_uniprot$UNIPROT[i]] <- out_uniprot$SYMBOL[i]
  if (nrow(out_alias) > 0L && "SYMBOL" %in% names(out_alias))
    for (i in seq_len(nrow(out_alias)))
      if (!is.na(out_alias$SYMBOL[i]) && nzchar(out_alias$SYMBOL[i])) {
        a <- out_alias$ALIAS[i]
        if (a %in% names(symbol_by_id) && !nzchar(symbol_by_id[a]))
          symbol_by_id[a] <- out_alias$SYMBOL[i]
      }

  # Map versioned IDs to symbol if flat form was found (e.g. NP_123.1 -> symbol from NP_123)
  for (id in ids) {
    if (nzchar(symbol_by_id[id])) next
    flat <- sub("\\.[0-9]+$", "", id)
    if (flat != id && nzchar(symbol_by_id[flat]))
      symbol_by_id[id] <- symbol_by_id[flat]
  }

  # Retain original ID when mapping fails (paper requirement)
  data.frame(
    ProteinName = ids,
    GeneSymbol = ifelse(nzchar(symbol_by_id[ids]), symbol_by_id[ids], ids),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 11) QC summary and intensity distribution (before/after summarization)
# ------------------------------------------------------------------------------
qc_summary <- function(parsed_psm, msstats_input, protein_summary, gene_matrix, outdir, annotation_audit_info = NULL) {
  prot_n <- if ("Protein" %in% names(protein_summary)) uniqueN(protein_summary$Protein) else uniqueN(protein_summary$ProteinName)
  total_runs <- uniqueN(parsed_psm$Run)
  total_fractions <- uniqueN(parsed_psm[, .(Run, Fraction)])
  lines <- c(
    "=== PDC/CPTAC PSM → MSstatsTMT → Gene Matrix QC ===",
    paste("Total runs parsed (fraction files):", total_runs),
    paste("Total fractions (Run x Fraction):", total_fractions),
    paste("PSM rows (long format):", nrow(parsed_psm)),
    paste("Peptides retained (unique PeptideSequence_Charge):", uniqueN(msstats_input[, .(PeptideSequence, Charge)])),
    paste("Proteins summarized:", prot_n),
    paste("Genes in final matrix:", nrow(gene_matrix)),
    "",
    "Intensity (PSM-level, before normalization):",
    paste("  Min:", min(parsed_psm$Intensity, na.rm = TRUE)),
    paste("  Max:", max(parsed_psm$Intensity, na.rm = TRUE)),
    paste("  Median:", median(parsed_psm$Intensity, na.rm = TRUE)),
    paste("  Mean:", round(mean(parsed_psm$Intensity, na.rm = TRUE), 2))
  )
  if ("Abundance" %in% names(protein_summary)) {
    lines <- c(lines,
      "",
      "Abundance (protein-level, after MSstatsTMT):",
      paste("  Min:", min(protein_summary$Abundance, na.rm = TRUE)),
      paste("  Max:", max(protein_summary$Abundance, na.rm = TRUE)),
      paste("  Median:", median(protein_summary$Abundance, na.rm = TRUE)),
      paste("  Mean:", round(mean(protein_summary$Abundance, na.rm = TRUE), 2))
    )
  }
  if (!is.null(annotation_audit_info)) {
    lines <- c(lines,
      "",
      "--- Annotation / design (from sample.txt) ---",
      paste("Total runs parsed:", annotation_audit_info$total_runs),
      paste("Total plexes detected:", annotation_audit_info$total_plexes),
      paste("Total unique mixtures:", annotation_audit_info$total_mixtures),
      paste("Total unique BioReplicates:", annotation_audit_info$total_bioreplicates),
      paste("Expected non-POOL sample count (from sample.txt):", annotation_audit_info$expected_non_pool_count),
      paste("Actual non-Norm BioReplicate count:", annotation_audit_info$actual_non_norm_count),
      paste("Fractions per mixture (min-max):", annotation_audit_info$fractions_per_plex),
      paste("BioReplicate stable across fractions:", annotation_audit_info$biorep_stable_across_fractions),
      paste("Bridge channels:", annotation_audit_info$bridge_channels),
      paste("Annotation audit result:", annotation_audit_info$annotation_mismatches)
    )
  }
  qc_path <- file.path(outdir, "qc_summary.txt")
  writeLines(lines, qc_path)
  message("QC summary written to ", qc_path)
}

# ------------------------------------------------------------------------------
# 12) Main
# ------------------------------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  psm_dir <- NULL
  annotation_path <- NULL
  reference_channel <- NULL
  outdir <- "."
  species <- "Hs"
  MBimpute <- FALSE
  max_runs <- NULL
  force_parse <- FALSE
  force_summarize <- FALSE
  sample_txt_path <- NULL
  replace_annotation <- FALSE

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--psm_dir" && i < length(args)) { psm_dir <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--annotation" && i < length(args)) { annotation_path <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--reference_channel" && i < length(args)) { reference_channel <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--outdir" && i < length(args)) { outdir <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--species" && i < length(args)) { species <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--MBimpute" && i < length(args)) { MBimpute <- as.logical(args[i + 1]); i <- i + 2 }
    else if (args[i] == "--max_runs" && i < length(args)) { max_runs <- as.integer(args[i + 1]); i <- i + 2 }
    else if (args[i] == "--force_parse") { force_parse <- TRUE; i <- i + 1 }
    else if (args[i] == "--force_summarize") { force_summarize <- TRUE; i <- i + 1 }
    else if (args[i] == "--sample_txt" && i < length(args)) { sample_txt_path <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--replace_annotation") { replace_annotation <- TRUE; i <- i + 1 }
    else i <- i + 1
  }

  if (is.null(psm_dir) || !dir.exists(psm_dir))
    stop("--psm_dir must point to a directory containing .psm files.")

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  parsed_path <- file.path(outdir, "parsed_psm_long.tsv")

  # Step 1: Parse PSM (or load from checkpoint)
  if (!force_parse && file.exists(parsed_path)) {
    message("Loading parsed PSM from checkpoint: ", parsed_path)
    parsed_psm <- fread(parsed_path, sep = "\t", header = TRUE)
  } else {
    parsed_psm <- load_all_psm(psm_dir)
    fwrite(parsed_psm, parsed_path, sep = "\t")
    message("Saved ", nrow(parsed_psm), " rows to ", parsed_path)
  }

  if (!is.null(max_runs) && max_runs > 0L) {
    runs_keep <- unique(parsed_psm$Run)[seq_len(min(max_runs, uniqueN(parsed_psm$Run)))]
    parsed_psm <- parsed_psm[Run %in% runs_keep]
    message("Subsampled to ", length(runs_keep), " runs (--max_runs ", max_runs, ")")
  }

  # Prefer annotation in outdir if no --annotation given
  if (is.null(annotation_path)) {
    ann_in_outdir <- file.path(outdir, "annotation_filled.csv")
    if (file.exists(ann_in_outdir)) annotation_path <- ann_in_outdir
  }

  # Step 2: Annotation — from CPTAC sample.txt (build/audit/rebuild) or require as before
  annotation_audit_info <- NULL
  if (!is.null(sample_txt_path) && file.exists(sample_txt_path)) {
    sample_design <- load_sample_txt(sample_txt_path)
    design_dt <- sample_design$design_dt
    n_plex_design <- nrow(design_dt)
    n_chan <- length(sample_design$channel_cols)
    expected_non_pool <- n_plex_design * (n_chan - 1L)

    if (is.null(annotation_path) || !file.exists(annotation_path)) {
      message("No existing annotation; building annotation from sample.txt and parsed PSM...")
      annotation <- rebuild_annotation_from_sample_txt(parsed_psm, sample_design)
      corrected_path <- file.path(outdir, "annotation_filled_corrected.csv")
      fwrite(annotation, corrected_path, sep = ",")
      message("Annotation written to ", corrected_path)
      if (replace_annotation) {
        fill_path <- file.path(outdir, "annotation_filled.csv")
        fwrite(annotation, fill_path, sep = ",")
        message("Overwritten ", fill_path, " (--replace_annotation)")
      }
      audit_result <- list(audit_ok = TRUE)
      # Write a short audit note (no "current" to compare)
      writeLines(c("=== Annotation built from sample.txt (no prior annotation) ===", "", "Result: Annotation created from CPTAC sample.txt.", paste("Rows:", nrow(annotation))), file.path(outdir, "annotation_audit.txt"))
    } else {
      message("Loading current annotation and study design (sample.txt) for audit...")
      annotation_current <- fread(annotation_path, header = TRUE)
      need <- c("Run", "Channel", "Condition", "BioReplicate", "Mixture", "Fraction", "TechRepMixture")
      if (length(setdiff(need, names(annotation_current))) > 0L)
        stop("When using --sample_txt, existing annotation must have columns: ", paste(need, collapse = ", "))
      audit_result <- audit_annotation(annotation_current, sample_design, parsed_psm, outdir)
      if (!audit_result$audit_ok) {
        message("Annotation audit failed; rebuilding from sample.txt...")
        annotation <- rebuild_annotation_from_sample_txt(parsed_psm, sample_design)
        corrected_path <- file.path(outdir, "annotation_filled_corrected.csv")
        fwrite(annotation, corrected_path, sep = ",")
        message("Corrected annotation written to ", corrected_path)
        if (replace_annotation) {
          fwrite(annotation, annotation_path, sep = ",")
          message("Overwritten ", annotation_path, " (--replace_annotation)")
        }
      } else {
        annotation <- as.data.table(annotation_current)
      }
    }

    validate_normalization_and_audit(annotation, sample_design, outdir)

    n_runs <- uniqueN(parsed_psm$Run)
    n_plex <- uniqueN(annotation$Mixture)
    n_mixtures <- uniqueN(annotation$Mixture)
    n_biorep <- uniqueN(annotation$BioReplicate)
    actual_non_norm <- uniqueN(annotation[tolower(Condition) != "norm", BioReplicate])
    bridge_dt <- annotation[tolower(Condition) == "norm", .(Mixture, Channel)]
    bridge <- if (nrow(bridge_dt) > 0L) paste(unique(paste(bridge_dt$Mixture, bridge_dt$Channel, sep = ":")), collapse = "; ") else "none"
    fpp <- annotation[, .(n = uniqueN(Fraction)), by = Mixture]$n
    frac_per_plex <- if (length(fpp) > 0L) paste(min(fpp), max(fpp), sep = "-") else "1"
    stable <- annotation[, .(n_bio = uniqueN(BioReplicate)), by = .(Mixture, Channel)][n_bio > 1L]
    biorep_stable <- if (nrow(stable) == 0L) "Yes" else "No"
    annotation_audit_info <- list(
      total_runs = n_runs,
      total_plexes = n_plex,
      total_mixtures = n_mixtures,
      total_bioreplicates = n_biorep,
      expected_non_pool_count = expected_non_pool,
      actual_non_norm_count = actual_non_norm,
      fractions_per_plex = frac_per_plex,
      biorep_stable_across_fractions = biorep_stable,
      bridge_channels = bridge,
      annotation_mismatches = if (!audit_result$audit_ok) "FAIL (corrected from sample.txt)" else "PASS"
    )
  } else {
    annotation <- require_annotation(annotation_path, parsed_psm, outdir, reference_channel)
  }

  # Step 3: Merge annotation; remove shared peptides
  merged <- merge_annotation(parsed_psm, annotation)
  merged <- remove_shared_peptides(merged)

  # MSstatsTMT required columns; aggregate to peptide level
  input_dt <- merged[, .(
    ProteinName, PeptideSequence, Charge, PSM,
    Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate, Fraction, Intensity
  )]
  input_dt <- input_dt[!is.na(Intensity) & Intensity > 0]
  input_dt <- input_dt[, .(
    Intensity = as.double(median(Intensity, na.rm = TRUE)),
    PSM = PSM[1L]
  ), by = c("ProteinName", "PeptideSequence", "Charge", "Mixture", "TechRepMixture", "Run", "Channel", "Condition", "BioReplicate", "Fraction")]
  if (nrow(input_dt) == 0L) stop("No rows after peptide-level aggregation.")
  input_dt[, Fraction := NULL]

  msstats_path <- file.path(outdir, "msstats_input.tsv")
  fwrite(input_dt, msstats_path, sep = "\t")
  message("Saved MSstatsTMT input to ", msstats_path)

  has_norm <- any(tolower(input_dt$Condition) == "norm")
  if (!has_norm) warning("No channel with Condition = 'Norm' (bridge). reference_norm will be FALSE; consider using --sample_txt for correct annotation.")

  prot_path <- file.path(outdir, "protein_summary.tsv")

  # Step 4: Protein summarization (or load from checkpoint). reference_norm only when valid bridge present.
  if (!force_summarize && file.exists(prot_path)) {
    message("Loading protein summary from checkpoint: ", prot_path)
    quant <- fread(prot_path, sep = "\t", header = TRUE)
    setDT(quant)
  } else {
    message("Running MSstatsTMT proteinSummarization...")
    quant <- proteinSummarization(
      input_dt,
      method = "msstats",
      global_norm = TRUE,
      reference_norm = has_norm,
      remove_norm_channel = TRUE,
      remove_empty_channel = TRUE,
      MBimpute = MBimpute
    )
    if (is.list(quant) && !inherits(quant, "data.frame")) {
      quant <- rbindlist(lapply(quant, as.data.table), use.names = TRUE, fill = TRUE)
    } else {
      quant <- as.data.table(quant)
    }
    if (nrow(quant) == 0L) stop("proteinSummarization returned no rows.")
    fwrite(quant, prot_path, sep = "\t")
    message("Saved protein summary to ", prot_path)
  }

  # Step 5: Map protein → gene symbol
  message("Mapping proteins to gene symbols...")
  prot_col <- if ("Protein" %in% names(quant)) "Protein" else "ProteinName"
  gene_map <- protein_to_gene(unique(quant[[prot_col]]), species)
  setnames(gene_map, "ProteinName", prot_col)
  quant <- merge(quant, gene_map, by = prot_col, all.x = TRUE)
  # Ensure both branches are character (avoid 'yes' integer / 'no' character error from select())
  quant[, GeneSymbol := fifelse(is.na(GeneSymbol) | GeneSymbol == "", as.character(get(prot_col)), as.character(GeneSymbol))]

  # Step 6: Build sample × gene matrix (BioReplicate = sample; exclude Norm channel)
  quant_norm_removed <- quant[tolower(Condition) != "norm"]
  sample_id <- quant_norm_removed$BioReplicate
  mat_dt <- quant_norm_removed[, .(sample = sample_id, GeneSymbol, Abundance)]
  mat_wide <- dcast(mat_dt, GeneSymbol ~ sample, value.var = "Abundance", fun.aggregate = median, na.rm = TRUE)
  mat <- as.matrix(mat_wide[, -1, with = FALSE])
  rownames(mat) <- mat_wide$GeneSymbol

  gene_path <- file.path(outdir, "gene_matrix.csv")
  fwrite(data.table(GeneSymbol = rownames(mat), mat), gene_path)
  message("Saved gene matrix (", nrow(mat), " x ", ncol(mat), ") to ", gene_path)

  # Step 7: QC (use full quant for abundance distribution; gene count from matrix)
  qc_summary(parsed_psm, input_dt, quant, data.table(GeneSymbol = rownames(mat)), outdir, annotation_audit_info)
}

tryCatch(main(), error = function(e) {
  message("Error: ", conditionMessage(e))
  quit(status = 1, save = "no")
})
