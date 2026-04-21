#!/usr/bin/env Rscript
# =============================================================================
# Method 0 — Raw Representation Adapter
# =============================================================================
# Returns the shared gene matrix unchanged (inner-join of CPTAC and CCLE genes,
# filtered by minimum observation fraction). No transformation applied.
#
# Usage:
#   source("scripts/methods/method_interface.R")
#   source("scripts/methods/run_raw_representation.R")
#   result <- run_raw_representation(cptac_mat, ccle_mat, cptac_meta, ccle_meta, outdir)
# =============================================================================

run_raw_representation <- function(cptac_mat, ccle_mat,
                                    cptac_meta, ccle_meta,
                                    outdir = "reports/benchmark_master/methods/raw",
                                    min_obs_frac = 0.1) {

  shared <- intersect_features(cptac_mat, ccle_mat, min_obs_frac = min_obs_frac)
  combined <- combine_domains(shared$mat_a, shared$mat_b, cptac_meta, ccle_meta)

  feature_meta <- data.frame(
    gene = c(combined$genes, shared$genes_dropped),
    included = c(rep(TRUE, length(combined$genes)), rep(FALSE, length(shared$genes_dropped))),
    exclusion_reason = c(rep(NA_character_, length(combined$genes)), shared$drop_reason),
    stringsAsFactors = FALSE
  )

  notes <- c(
    "Method: raw",
    paste0("Date: ", Sys.time()),
    "",
    "Description:",
    "  No transformation applied. The output matrix is the inner join of CPTAC",
    "  and CCLE gene matrices after filtering genes with >=10% observed values in",
    "  both domains. Values are log2 protein abundances from MSstatsTMT.",
    "",
    paste0("Genes included: ", length(combined$genes)),
    paste0("Genes excluded (low observation): ", length(shared$genes_dropped)),
    paste0("CPTAC samples: ", sum(combined$sample_meta$domain == "CPTAC")),
    paste0("CCLE samples:  ", sum(combined$sample_meta$domain == "CCLE")),
    "",
    "What 'raw' means:",
    "  - MSstatsTMT protein summarization with reference normalization (bridge channel)",
    "  - Gene symbol mapping and median collapse per gene",
    "  - No cross-domain alignment or batch correction",
    "  - Each domain retains its own scale, centering, and spread",
    "",
    "Limitations:",
    "  - Systematic domain differences (instrument, protocol, sample type) are present",
    "  - Direct cross-domain fold-change comparison is confounded by batch effects",
    "  - This method serves as a baseline to quantify how much correction is needed"
  )

  result <- make_method_result(
    matrix       = combined$matrix,
    sample_meta  = combined$sample_meta,
    feature_meta = feature_meta,
    method_name  = "raw",
    method_notes = notes,
    qc_paths     = c(notes_file = file.path(outdir, "method_notes.txt"))
  )
  save_method_result(result, outdir)
  result
}
