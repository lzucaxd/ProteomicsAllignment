#!/usr/bin/env Rscript
# =============================================================================
# Benchmark subset strategies
# =============================================================================
# Deterministic, documented subsetting for both benchmark tasks.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Map annotation sample_id to gene-matrix column names
# ---------------------------------------------------------------------------
# Prefer exact / case-insensitive equality, then a single regex hit (hyphens →
# "." for grep). Used so union metadata aligns with CCLE matrix colnames.
# ---------------------------------------------------------------------------
map_sample_id_to_matrix_col <- function(sid, matrix_cols) {
  if (!nzchar(sid)) return(NA_character_)
  if (sid %in% matrix_cols) return(sid)
  eq <- matrix_cols[tolower(matrix_cols) == tolower(sid)]
  if (length(eq) >= 1L) return(eq[1L])
  pat <- gsub("-", ".", sid, fixed = TRUE)
  m <- grep(pat, matrix_cols, ignore.case = TRUE, value = TRUE)
  if (length(m) == 1L) return(m[1L])
  if (length(m) > 1L) {
    ex <- matrix_cols[tolower(matrix_cols) == tolower(sid)]
    if (length(ex) == 1L) return(ex[1L])
    return(NA_character_)
  }
  NA_character_
}

map_sample_ids_to_matrix_cols <- function(ids, matrix_cols) {
  vapply(ids, map_sample_id_to_matrix_col, character(1L), matrix_cols = matrix_cols)
}

# ---------------------------------------------------------------------------
# Task A: Breast subtype — mixture-balanced Basal vs Luminal subset (CPTAC)
# ---------------------------------------------------------------------------
build_subtype_subset_cptac <- function(subtype_mapping_path, outdir = NULL) {
  sm <- fread(subtype_mapping_path)

  required <- c("matrix_sample_id", "pam50", "mixture", "sample_type")
  missing <- setdiff(required, names(sm))
  if (length(missing)) stop("Subtype mapping missing columns: ", paste(missing, collapse = ", "))

  tumors <- sm[tolower(sample_type) == "tumor" &
               tolower(pam50) %in% c("basal", "luma", "lumb")]
  tumors[, subtype := ifelse(tolower(pam50) == "basal", "Basal", "Luminal")]

  mix_counts <- tumors[, .(n_Basal = sum(subtype == "Basal"),
                           n_Luminal = sum(subtype == "Luminal"),
                           total = .N), by = mixture]
  # Rule: drop if either subtype absent
  mix_counts[, keep := (n_Basal >= 1) & (n_Luminal >= 1)]
  # Rule: drop near-absent minority
  mix_counts[keep == TRUE, keep := !(pmin(n_Basal, n_Luminal) == 1 & total >= 6)]

  kept_mixtures <- mix_counts[keep == TRUE, mixture]
  subset_dt <- tumors[mixture %in% kept_mixtures]

  summary_lines <- c(
    "=== Breast Subtype Subset (CPTAC PDC000120) ===",
    "",
    paste("Full tumor Basal:", tumors[subtype == "Basal", .N]),
    paste("Full tumor Luminal:", tumors[subtype == "Luminal", .N]),
    paste("Mixtures (full):", uniqueN(tumors$mixture)),
    "",
    paste("Subset Basal:", subset_dt[subtype == "Basal", .N]),
    paste("Subset Luminal:", subset_dt[subtype == "Luminal", .N]),
    paste("Mixtures (subset):", uniqueN(subset_dt$mixture)),
    paste("Mixtures dropped:", uniqueN(tumors$mixture) - uniqueN(subset_dt$mixture)),
    "",
    "Rules applied:",
    "  1. Tumor-only with PAM50 in {Basal, LumA, LumB}",
    "  2. LumA + LumB collapsed to Luminal",
    "  3. Drop mixture if n_Basal == 0 OR n_Luminal == 0",
    "  4. Drop mixture if min(n_Basal, n_Luminal) == 1 AND total >= 6"
  )

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    fwrite(mix_counts, file.path(outdir, "subtype_subset_summary.tsv"), sep = "\t")
    writeLines(summary_lines, file.path(outdir, "subtype_subset_notes.txt"))
    fwrite(subset_dt, file.path(outdir, "subtype_subset_samples.csv"))
    message("Subtype subset saved to ", outdir)
  }

  list(
    subset = subset_dt,
    mixture_summary = mix_counts,
    notes = summary_lines
  )
}

# ---------------------------------------------------------------------------
# Task A: Breast subtype — CCLE from union metadata (full panel, default)
# ---------------------------------------------------------------------------
# Reads data/processed/union/sample_meta_breast_subtype.csv (or alternate path),
# keeps Basal/Luminal CCLE rows, maps sample_id → actual matrix column names.
# sample_id in the returned table is the matrix column name (for intersect()).
#
# Legacy 8-line panel: pass union_meta_path = NULL and ccle_sample_info_path.
# ---------------------------------------------------------------------------
build_subtype_subset_ccle <- function(ccle_matrix_path,
                                      union_meta_path = NULL,
                                      ccle_sample_info_path = NULL,
                                      outdir = NULL) {
  cols <- setdiff(names(fread(ccle_matrix_path, nrows = 0L)),
                  c("GeneSymbol", "UniProtID", "Gene"))

  if (!is.null(union_meta_path) && nzchar(union_meta_path) && file.exists(union_meta_path)) {
    um <- fread(union_meta_path)
    need <- c("sample_id", "domain", "condition")
    miss <- setdiff(need, names(um))
    if (length(miss)) stop("Union subtype meta missing columns: ", paste(miss, collapse = ", "))

    ccle <- um[toupper(domain) == "CCLE" &
                 tolower(condition) %in% c("basal", "luminal")]
    if (nrow(ccle) < 4L) stop("Union meta has fewer than 4 CCLE Basal/Luminal rows")

    ccle[, subtype := fifelse(tolower(condition) == "basal", "Basal", "Luminal")]
    ccle[, matrix_col := map_sample_ids_to_matrix_cols(sample_id, cols)]

    dropped <- ccle[is.na(matrix_col)]
    if (nrow(dropped) > 0L) {
      message("  CCLE subtype: dropped ", nrow(dropped),
              " union row(s) with no matrix column match (e.g. ",
              paste(head(dropped$sample_id, 3L), collapse = ", "), ")")
    }
    ccle <- ccle[!is.na(matrix_col)]
    ccle <- unique(ccle, by = "matrix_col")
    subset_dt <- ccle[, .(sample_id = matrix_col, subtype, domain = "CCLE")]

    summary_lines <- c(
      "=== Breast Subtype Subset (CCLE) — union metadata ===",
      "",
      paste("Union meta path:", union_meta_path),
      paste("Basal lines (in matrix):", subset_dt[subtype == "Basal", .N]),
      paste("Luminal lines (in matrix):", subset_dt[subtype == "Luminal", .N]),
      paste("Total:", nrow(subset_dt))
    )
  } else {
    if (is.null(ccle_sample_info_path) || !file.exists(ccle_sample_info_path)) {
      stop("Provide union_meta_path to sample_meta_breast_subtype.csv, or a valid ",
           "ccle_sample_info_path for the legacy 8-line panel.")
    }
    basal_lines <- c("HCC70", "HCC 1806", "HCC1143", "MDA-MB-468")
    luminal_lines <- c("CAMA-1", "MCF7", "T-47D", "ZR-75-1")
    bmc <- map_sample_ids_to_matrix_cols(basal_lines, cols)
    lmc <- map_sample_ids_to_matrix_cols(luminal_lines, cols)
    subset_dt <- rbind(
      data.table(sample_id = bmc[!is.na(bmc)], subtype = "Basal", domain = "CCLE"),
      data.table(sample_id = lmc[!is.na(lmc)], subtype = "Luminal", domain = "CCLE")
    )
    summary_lines <- c(
      "=== Breast Subtype Subset (CCLE) — legacy 8-line panel ===",
      "",
      paste("Basal:", subset_dt[subtype == "Basal", .N]),
      paste("Luminal:", subset_dt[subtype == "Luminal", .N])
    )
  }

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    fwrite(subset_dt, file.path(outdir, "subtype_subset_ccle_samples.csv"))
    writeLines(summary_lines, file.path(outdir, "subtype_subset_ccle_notes.txt"))
    message("CCLE subtype subset saved to ", outdir)
  }

  subset_dt
}

# ---------------------------------------------------------------------------
# Task B: Breast vs Lung — CPTAC (all tumors from PDC000120 + PDC000153)
# ---------------------------------------------------------------------------
build_breast_vs_lung_subset_cptac <- function(breast_matrix_path,
                                                lung_matrix_path,
                                                breast_subtype_mapping_path = NULL,
                                                outdir = NULL) {
  breast_dt <- fread(breast_matrix_path, nrows = 0)
  lung_dt <- fread(lung_matrix_path, nrows = 0)

  breast_samples <- setdiff(names(breast_dt), c("GeneSymbol", "UniProtID"))
  lung_samples <- setdiff(names(lung_dt), c("GeneSymbol", "UniProtID"))

  subset_dt <- data.table(
    sample_id = c(breast_samples, lung_samples),
    cancer_type = c(rep("Breast", length(breast_samples)),
                    rep("Lung", length(lung_samples))),
    study_id = c(rep("PDC000120", length(breast_samples)),
                 rep("PDC000153", length(lung_samples))),
    domain = "CPTAC"
  )

  summary_lines <- c(
    "=== Breast vs Lung Subset (CPTAC) ===",
    "",
    paste("Breast samples (PDC000120):", length(breast_samples)),
    paste("Lung samples (PDC000153):", length(lung_samples)),
    "",
    "Note: Cancer type is perfectly confounded with study.",
    "All representations share this confound equally."
  )

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    fwrite(subset_dt, file.path(outdir, "breast_vs_lung_subset_samples.csv"))
    fwrite(data.table(
      domain = "CPTAC",
      breast_n = length(breast_samples),
      lung_n = length(lung_samples),
      breast_study = "PDC000120",
      lung_study = "PDC000153",
      confound_note = "cancer_type == study_id"
    ), file.path(outdir, "breast_vs_lung_subset_summary.tsv"), sep = "\t")
    writeLines(summary_lines, file.path(outdir, "breast_vs_lung_subset_notes.txt"))
    message("Breast vs Lung subset saved to ", outdir)
  }

  list(subset = subset_dt, notes = summary_lines)
}

# ---------------------------------------------------------------------------
# Task B: Breast vs Lung — CCLE
# ---------------------------------------------------------------------------
build_breast_vs_lung_subset_ccle <- function(ccle_sample_info_path, outdir = NULL) {
  info <- fread(ccle_sample_info_path)
  name_col <- if ("Cell Line" %in% names(info)) "Cell Line" else names(info)[1]
  tissue_col <- if ("Tissue of Origin" %in% names(info)) "Tissue of Origin" else "tissue"

  breast <- info[tolower(get(tissue_col)) == "breast", get(name_col)]
  lung <- info[tolower(get(tissue_col)) == "lung", get(name_col)]

  subset_dt <- data.table(
    sample_id = c(breast, lung),
    cancer_type = c(rep("Breast", length(breast)), rep("Lung", length(lung))),
    domain = "CCLE"
  )

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    fwrite(subset_dt, file.path(outdir, "breast_vs_lung_ccle_samples.csv"))
  }

  list(subset = subset_dt,
       notes = c(paste("CCLE Breast:", length(breast)),
                 paste("CCLE Lung:", length(lung))))
}
