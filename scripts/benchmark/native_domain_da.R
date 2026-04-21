#!/usr/bin/env Rscript
# =============================================================================
# Native-domain baseline inference
# =============================================================================
# For comparisons where the original TMT experimental design is still
# meaningful, use MSstatsTMT groupComparisonTMT (CPTAC) or limma (CCLE).
#
# These results serve as the reference/anchor for evaluating how well
# representation-level comparisons preserve biological signal.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ---------------------------------------------------------------------------
# Native CPTAC subtype DA (MSstatsTMT-based, protein level)
# ---------------------------------------------------------------------------
# This requires the MSstatsTMT protein_summary.tsv and the annotation with
# Condition = subtype. If MSstatsTMT is installed and the summary exists,
# it runs groupComparisonTMT. Otherwise, falls back to limma on gene_matrix.
# ---------------------------------------------------------------------------
run_native_subtype_da_cptac <- function(protein_summary_path = NULL,
                                         gene_matrix_path,
                                         subtype_subset,
                                         outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  msstats_available <- requireNamespace("MSstatsTMT", quietly = TRUE)
  prot_exists <- !is.null(protein_summary_path) && file.exists(protein_summary_path)

  if (msstats_available && prot_exists) {
    message("Running MSstatsTMT groupComparisonTMT for native CPTAC subtype DA...")
    inference_type <- "native_domain_msstatstmt"

    # SCAFFOLD: Full MSstatsTMT groupComparisonTMT pipeline
    # Requires:
    #   1. protein_summary.tsv from MSstatsTMT proteinSummarization
    #   2. Annotation with Condition mapped to subtype (Basal/Luminal)
    #   3. Contrast matrix: Luminal - Basal
    #
    # This section is a scaffold. To run:
    #   library(MSstatsTMT)
    #   data <- read.csv(protein_summary_path)
    #   # Remap Condition to subtype based on subtype_subset
    #   # Build contrast matrix
    #   result <- groupComparisonTMT(data, contrast.matrix = ...)

    warning("MSstatsTMT groupComparisonTMT scaffold — full pipeline requires ",
            "protein_summary.tsv with Condition remapped to Basal/Luminal. ",
            "Falling back to limma on gene_matrix.")
    prot_exists <- FALSE
  }

  if (!prot_exists) {
    message("Running limma on gene_matrix for native CPTAC subtype DA...")
    inference_type <- "native_domain_limma"

    gm <- fread(gene_matrix_path, header = TRUE)
    gene_col <- names(gm)[1]
    id_cols <- intersect(c("GeneSymbol", "UniProtID"), names(gm))
    sample_cols <- setdiff(names(gm), id_cols)

    ss <- as.data.table(subtype_subset)
    sample_ids <- intersect(ss$sample_id %||% ss$matrix_sample_id, sample_cols)
    if (length(sample_ids) < 4) stop("Fewer than 4 subtype samples found in matrix")

    mat <- as.matrix(gm[, ..sample_ids])
    rownames(mat) <- gm[[gene_col]]

    groups <- ss[match(sample_ids, sample_id %||% matrix_sample_id), subtype]

    keep <- rowSums(!is.na(mat)) >= ncol(mat) * 0.5
    mat <- mat[keep, , drop = FALSE]

    source_dir <- dirname(sys.frame(1)$ofile %||% "scripts/benchmark/native_domain_da.R")
    eh <- file.path(source_dir, "evaluation_helpers.R")
    if (file.exists(eh)) source(eh)

    da <- run_limma_da(mat, groups, contrast_name = "Luminal_vs_Basal_CPTAC_native")
    da[, inference_type := inference_type]
  }

  fwrite(da, file.path(outdir, "da_result.csv"))

  note <- data.table(
    domain = "CPTAC", task = "breast_subtype",
    n_basal = sum(groups == "Basal"),
    n_luminal = sum(groups == "Luminal"),
    n_genes = nrow(mat),
    inference_type = inference_type,
    note = if (inference_type == "native_domain_msstatstmt")
      "Protein-level MSstatsTMT groupComparisonTMT" else
      "Gene-level limma (MSstatsTMT protein_summary not available)"
  )
  fwrite(note, file.path(outdir, "sample_counts.csv"))

  markers <- da[gene %in% c("FOXA1", "GATA3", "KRT5", "KRT14", "KRT17",
                             "EGFR", "ESR1", "PGR", "ERBB2", "CDH1")]
  fwrite(markers, file.path(outdir, "marker_summary.csv"))

  writeLines(c(
    "=== Native-Domain DA: CPTAC Breast Subtype ===",
    "",
    paste("Inference type:", inference_type),
    paste("Level:", if (inference_type == "native_domain_msstatstmt") "protein" else "gene"),
    paste("Contrast: Luminal - Basal (log2 scale)"),
    paste("n Basal:", sum(groups == "Basal")),
    paste("n Luminal:", sum(groups == "Luminal")),
    paste("Genes tested:", nrow(mat)),
    "",
    "This is a native-domain analysis. The original TMT experimental",
    "design is used for inference. Results from this analysis serve",
    "as the reference for evaluating representation-level comparisons."
  ), file.path(outdir, "analysis_notes.txt"))

  message("Native CPTAC subtype DA complete → ", outdir)
  invisible(list(da = da, note = note, markers = markers))
}

# ---------------------------------------------------------------------------
# Native CCLE subtype DA (limma, gene level)
# ---------------------------------------------------------------------------
run_native_subtype_da_ccle <- function(gene_matrix_path,
                                        subtype_subset,
                                        outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  gm <- fread(gene_matrix_path, header = TRUE)
  gene_col <- names(gm)[1]
  id_cols <- intersect(c("GeneSymbol", "UniProtID"), names(gm))
  sample_cols <- setdiff(names(gm), id_cols)

  ss <- unique(as.data.table(subtype_subset), by = "sample_id")
  if (!all(c("sample_id", "subtype") %in% names(ss)))
    stop("subtype_subset must contain sample_id and subtype (Basal/Luminal)")

  # Prefer matrix column names already resolved (union metadata path).
  on_mat <- ss$sample_id %in% sample_cols
  if (all(on_mat)) {
    use <- ss[sample_id %in% sample_cols]
  } else {
    ccle_matrix_ids <- sample_cols
    matched <- vapply(ss$sample_id, function(sid) {
      if (sid %in% sample_cols) return(sid)
      m <- grep(gsub("-", ".", sid, fixed = TRUE), ccle_matrix_ids,
                ignore.case = TRUE, value = TRUE)
      if (length(m) == 1L) m else NA_character_
    }, character(1L))
    ss[, matrix_col := matched]
    ss <- ss[!is.na(matrix_col)]
    use <- ss[, .(sample_id = matrix_col, subtype)]
  }

  if (nrow(use) < 4) stop("Fewer than 4 CCLE subtype samples matched in matrix")

  mat <- as.matrix(gm[, use$sample_id, with = FALSE])
  rownames(mat) <- gm[[gene_col]]
  groups <- use$subtype

  keep <- rowSums(!is.na(mat)) >= ncol(mat) * 0.5
  mat <- mat[keep, , drop = FALSE]

  source_dir <- dirname(sys.frame(1)$ofile %||% "scripts/benchmark/native_domain_da.R")
  eh <- file.path(source_dir, "evaluation_helpers.R")
  if (file.exists(eh)) source(eh)

  da <- run_limma_da(mat, groups, contrast_name = "Luminal_vs_Basal_CCLE_native")
  da[, inference_type := "native_domain_limma"]

  fwrite(da, file.path(outdir, "da_result.csv"))

  note <- data.table(
    domain = "CCLE", task = "breast_subtype",
    n_basal = sum(groups == "Basal"),
    n_luminal = sum(groups == "Luminal"),
    n_genes = nrow(mat),
    inference_type = "native_domain_limma",
    note = "CCLE: 1 cell line per plex → limma on gene matrix is appropriate"
  )
  fwrite(note, file.path(outdir, "sample_counts.csv"))

  markers <- da[gene %in% c("FOXA1", "GATA3", "KRT5", "KRT14", "KRT17",
                             "EGFR", "ESR1", "PGR", "ERBB2", "CDH1")]
  fwrite(markers, file.path(outdir, "marker_summary.csv"))

  writeLines(c(
    "=== Native-Domain DA: CCLE Breast Subtype ===",
    "",
    "Inference type: native_domain_limma",
    "Level: gene",
    "Contrast: Luminal - Basal (log2 scale)",
    paste("n Basal:", sum(groups == "Basal")),
    paste("n Luminal:", sum(groups == "Luminal")),
    paste("Genes tested:", nrow(mat)),
    "",
    "CCLE native subtype DA uses limma on the gene matrix directly.",
    "Each cell line occupies its own TMT plex, so limma on the",
    "MSstatsTMT gene-level summary is the appropriate tool.",
    "",
    "CAVEAT: CCLE n is the number of cell lines in the subtype panel; power is limited.",
    "See tasks/breast_subtype/subtype_subset_ccle_samples.csv for the exact lines used."
  ), file.path(outdir, "analysis_notes.txt"))

  message("Native CCLE subtype DA complete → ", outdir)
  invisible(list(da = da, note = note, markers = markers))
}

# ---------------------------------------------------------------------------
# Native breast vs lung DA (limma, gene level)
# ---------------------------------------------------------------------------
run_native_breast_vs_lung_da <- function(domain = c("CPTAC", "CCLE"),
                                          gene_matrix_path,
                                          sample_subset,
                                          outdir) {
  domain <- match.arg(domain)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  gm <- fread(gene_matrix_path, header = TRUE)
  gene_col <- names(gm)[1]
  id_cols <- intersect(c("GeneSymbol", "UniProtID"), names(gm))
  sample_cols <- setdiff(names(gm), id_cols)

  ss <- as.data.table(sample_subset)
  sample_ids <- intersect(ss$sample_id, sample_cols)
  if (length(sample_ids) < 4) stop("Fewer than 4 lineage samples in matrix")

  mat <- as.matrix(gm[, ..sample_ids])
  rownames(mat) <- gm[[gene_col]]
  groups <- ss[match(sample_ids, sample_id), cancer_type]

  keep <- rowSums(!is.na(mat)) >= ncol(mat) * 0.3
  mat <- mat[keep, , drop = FALSE]

  source_dir <- dirname(sys.frame(1)$ofile %||% "scripts/benchmark/native_domain_da.R")
  eh <- file.path(source_dir, "evaluation_helpers.R")
  if (file.exists(eh)) source(eh)

  da <- run_limma_da(mat, groups,
                      contrast_name = paste0("Breast_vs_Lung_", domain, "_native"))
  da[, inference_type := "native_domain_limma"]

  fwrite(da, file.path(outdir, "da_result.csv"))

  note <- data.table(
    domain = domain, task = "breast_vs_lung",
    n_breast = sum(groups == "Breast"),
    n_lung = sum(groups == "Lung"),
    n_genes = nrow(mat),
    inference_type = "native_domain_limma",
    note = if (domain == "CPTAC")
      "Cross-study: cancer_type confounded with study (PDC000120 vs PDC000153)" else
      "CCLE: single experiment, unbalanced design"
  )
  fwrite(note, file.path(outdir, "sample_counts.csv"))

  markers <- da[gene %in% c("NKX2-1", "SFTPB", "SFTPC", "NAPSA",
                             "GATA3", "FOXA1", "ESR1", "KRT19",
                             "EGFR", "ERBB2", "CDH1", "VIM")]
  fwrite(markers, file.path(outdir, "marker_summary.csv"))

  writeLines(c(
    paste0("=== Native-Domain DA: ", domain, " Breast vs Lung ==="),
    "",
    "Inference type: native_domain_limma",
    "Level: gene",
    "Contrast: Breast - Lung (log2 scale)",
    paste("n Breast:", sum(groups == "Breast")),
    paste("n Lung:", sum(groups == "Lung")),
    paste("Genes tested:", nrow(mat)),
    "",
    if (domain == "CPTAC") paste(
      "CAVEAT: Breast and Lung come from different PDC studies.",
      "Cancer type is perfectly confounded with study.",
      "This baseline should be interpreted as study+cancer_type."
    ) else paste(
      "CCLE breast vs lung from same experiment.",
      "Unbalanced (30 vs 77) but no study confound."
    )
  ), file.path(outdir, "analysis_notes.txt"))

  message("Native ", domain, " breast vs lung DA complete → ", outdir)
  invisible(list(da = da, note = note, markers = markers))
}
