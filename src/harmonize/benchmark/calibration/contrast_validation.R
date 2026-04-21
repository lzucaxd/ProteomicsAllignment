#!/usr/bin/env Rscript
# =============================================================================
# Contrast Validation — check marker logFC signs against expected directions
# =============================================================================
# Reusable function: source this file and call validate_contrast_direction().
# Depends on: evaluation_helpers.R (for run_limma_da)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

validate_contrast_direction <- function(matrix, groups, expected_directions,
                                         contrast_label = "auto") {
  # ---------------------------------------------------------------------------
  # Args:
  #   matrix          : numeric matrix, genes x samples
  #   groups          : character/factor vector of condition labels (length = ncol)
  #   expected_directions : data.table/data.frame with columns:
  #                         gene (char), expected_sign (+1 or -1)
  #   contrast_label  : human-readable string (e.g. "Luminal_vs_Basal")
  #
  # Returns list:
  #   $marker_table   : per-gene data.table
  #   $summary        : aggregate stats
  #   $likely_flipped : TRUE if fraction_correct < 0.5
  # ---------------------------------------------------------------------------

  if (length(groups) != ncol(matrix))
    stop("groups length (", length(groups), ") != matrix columns (", ncol(matrix), ")")

  expected_directions <- as.data.table(expected_directions)
  stopifnot(all(c("gene", "expected_sign") %in% names(expected_directions)))

  group_factor <- factor(groups)
  lvls <- levels(group_factor)
  if (length(lvls) != 2)
    stop("Exactly 2 group levels required, got: ", paste(lvls, collapse = ", "))

  # Explicit contrast: lvls[2] - lvls[1] (alphabetical)
  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- lvls
  contrast_str <- paste0(lvls[2], " - ", lvls[1])

  if (contrast_label == "auto")
    contrast_label <- paste0(lvls[2], "_vs_", lvls[1])

  fit <- lmFit(matrix, design)
  cm <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  tt <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"))
  tt$gene <- rownames(topTable(fit2, number = Inf, sort.by = "none"))

  # Merge with expected markers
  marker_table <- merge(
    tt[, .(gene, logFC, t, P.Value, adj.P.Val)],
    expected_directions,
    by = "gene"
  )

  if (nrow(marker_table) == 0) {
    return(list(
      marker_table = marker_table,
      summary = list(
        contrast_label = contrast_label,
        contrast_direction = contrast_str,
        n_markers_tested = 0L,
        n_correct = 0L,
        fraction_correct = NA_real_
      ),
      likely_flipped = NA
    ))
  }

  marker_table[, observed_sign := sign(logFC)]
  marker_table[, correct := (observed_sign == expected_sign)]

  n_tested <- nrow(marker_table)
  n_correct <- sum(marker_table$correct)
  frac <- n_correct / n_tested

  summary_list <- list(
    contrast_label = contrast_label,
    contrast_direction = contrast_str,
    n_markers_tested = n_tested,
    n_correct = n_correct,
    fraction_correct = frac,
    levels_order = paste(lvls, collapse = " < ")
  )

  likely_flipped <- frac < 0.5

  list(
    marker_table = marker_table,
    summary = summary_list,
    likely_flipped = likely_flipped
  )
}
