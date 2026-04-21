#!/usr/bin/env Rscript
# =============================================================================
# Within-Domain Concordance Ceiling
# =============================================================================
# Estimates the maximum achievable FC agreement *within* a single domain by
# randomly splitting samples into two halves, running limma on each half, and
# measuring cross-half FC correlation / same-direction fraction.
#
# For small-n domains (e.g. CCLE, <=12 per condition), a leave-one-out
# jackknife variant is used instead of random splits.
#
# Usage (CLI):
#   Rscript concordance_ceiling.R \
#     --matrix input.csv --meta meta.csv \
#     --domain CPTAC --contrast-a Basal --contrast-b Luminal \
#     --n-splits 100 --seed 42 --outdir /path
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

compute_concordance_ceiling <- function(matrix, sample_meta, domain,
                                         contrast_a, contrast_b,
                                         n_splits = 100, seed = 42,
                                         jackknife_threshold = 12,
                                         intersection_genes = NULL,
                                         force_split_half = FALSE) {
  set.seed(seed)

  meta <- as.data.table(sample_meta)
  setnames(meta, "domain", "domain_col", skip_absent = TRUE)
  dom_meta <- meta[toupper(domain_col) == toupper(domain) &
                   condition %in% c(contrast_a, contrast_b)]
  dom_cols <- intersect(dom_meta$sample_id, colnames(matrix))
  dom_meta <- dom_meta[sample_id %in% dom_cols]

  if (nrow(dom_meta) < 8)
    stop("Too few samples in ", domain, " (", nrow(dom_meta), ") for ceiling estimation")

  dom_mat <- matrix[, dom_meta$sample_id, drop = FALSE]

  # Filter genes
  for (g in c(contrast_a, contrast_b)) {
    g_cols <- which(dom_meta$condition == g)
    if (length(g_cols) > 0) {
      na_frac <- rowMeans(is.na(dom_mat[, g_cols, drop = FALSE]))
      dom_mat[na_frac > 0.5, g_cols] <- NA
    }
  }
  keep <- rowSums(!is.na(dom_mat)) >= ncol(dom_mat) * 0.5
  dom_mat <- dom_mat[keep, , drop = FALSE]

  n_per_group <- min(table(dom_meta$condition))
  use_jackknife <- !isTRUE(force_split_half) && (n_per_group <= jackknife_threshold)

  run_half_da <- function(sids) {
    if (length(sids) < 4) return(NULL)
    half_mat <- dom_mat[, sids, drop = FALSE]
    half_groups <- dom_meta[match(sids, sample_id), condition]
    if (length(unique(half_groups)) < 2) return(NULL)

    gf <- factor(half_groups, levels = c(contrast_a, contrast_b))
    design <- model.matrix(~ 0 + gf)
    colnames(design) <- levels(gf)
    fit <- lmFit(half_mat, design)
    cm <- makeContrasts(contrasts = paste0(contrast_b, " - ", contrast_a), levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cm))
    tt <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"))
    tt$gene <- rownames(topTable(fit2, number = Inf, sort.by = "none"))
    tt
  }

  compare_halves <- function(da1, da2) {
    shared <- merge(da1[, .(gene, fc1 = logFC)],
                    da2[, .(gene, fc2 = logFC)], by = "gene")
    shared <- shared[is.finite(fc1) & is.finite(fc2)]
    if (!is.null(intersection_genes) && length(intersection_genes) > 0)
      shared <- shared[gene %in% intersection_genes]
    if (nrow(shared) < 10) return(NULL)

    data.table(
      fc_correlation = cor(shared$fc1, shared$fc2),
      same_direction_frac = mean(sign(shared$fc1) == sign(shared$fc2)),
      n_genes = nrow(shared)
    )
  }

  split_results <- vector("list", n_splits)

  if (use_jackknife) {
    cat("  Using leave-one-out jackknife for", domain,
        "(", n_per_group, "per condition)\n")

    sids_a <- dom_meta[condition == contrast_a, sample_id]
    sids_b <- dom_meta[condition == contrast_b, sample_id]
    all_sids <- c(sids_a, sids_b)
    n_loo <- length(all_sids)

    for (i in seq_len(min(n_loo, n_splits))) {
      left_out <- all_sids[i]
      remaining <- setdiff(all_sids, left_out)

      # Split remaining into two halves (balanced)
      rem_a <- intersect(remaining, sids_a)
      rem_b <- intersect(remaining, sids_b)
      half_size_a <- floor(length(rem_a) / 2)
      half_size_b <- floor(length(rem_b) / 2)

      sel_a1 <- sample(rem_a, half_size_a)
      sel_b1 <- sample(rem_b, half_size_b)
      sel_a2 <- setdiff(rem_a, sel_a1)
      sel_b2 <- setdiff(rem_b, sel_b1)

      da1 <- run_half_da(c(sel_a1, sel_b1))
      da2 <- run_half_da(c(sel_a2, sel_b2))

      if (!is.null(da1) && !is.null(da2)) {
        split_results[[i]] <- compare_halves(da1, da2)
        if (!is.null(split_results[[i]])) split_results[[i]]$split <- i
      }
    }
  } else {
    cat("  Using random split-half for", domain,
        "(", n_per_group, "per condition)\n")

    sids_a <- dom_meta[condition == contrast_a, sample_id]
    sids_b <- dom_meta[condition == contrast_b, sample_id]

    for (i in seq_len(n_splits)) {
      half_size_a <- floor(length(sids_a) / 2)
      half_size_b <- floor(length(sids_b) / 2)

      sel_a1 <- sample(sids_a, half_size_a)
      sel_b1 <- sample(sids_b, half_size_b)
      sel_a2 <- setdiff(sids_a, sel_a1)
      sel_b2 <- setdiff(sids_b, sel_b1)

      da1 <- run_half_da(c(sel_a1, sel_b1))
      da2 <- run_half_da(c(sel_a2, sel_b2))

      if (!is.null(da1) && !is.null(da2)) {
        split_results[[i]] <- compare_halves(da1, da2)
        if (!is.null(split_results[[i]])) split_results[[i]]$split <- i
      }

      if (i %% 25 == 0) cat("  Split", i, "/", n_splits, "\n")
    }
  }

  split_dt <- rbindlist(split_results[!sapply(split_results, is.null)])

  if (nrow(split_dt) == 0)
    stop("No valid split-half results for ", domain)

  ceiling_summary <- data.table(
    domain = domain,
    method_type = if (use_jackknife) "jackknife" else "split_half",
    n_valid_splits = nrow(split_dt),
    ceiling_fc_correlation = mean(split_dt$fc_correlation, na.rm = TRUE),
    ceiling_fc_corr_sd = sd(split_dt$fc_correlation, na.rm = TRUE),
    ceiling_same_dir_frac = mean(split_dt$same_direction_frac, na.rm = TRUE),
    ceiling_same_dir_sd = sd(split_dt$same_direction_frac, na.rm = TRUE),
    median_n_genes = median(split_dt$n_genes)
  )

  list(splits = split_dt, summary = ceiling_summary)
}


# ── CLI entry point ──────────────────────────────────────────────────────────
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  parse_arg <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) == 0) return(NULL)
    args[idx + 1]
  }

  matrix_path <- parse_arg("--matrix")
  meta_path   <- parse_arg("--meta")
  domain      <- parse_arg("--domain")
  contrast_a  <- parse_arg("--contrast-a")
  contrast_b  <- parse_arg("--contrast-b")
  outdir      <- parse_arg("--outdir")
  n_splits    <- as.integer(parse_arg("--n-splits") %||% "100")
  seed        <- as.integer(parse_arg("--seed") %||% "42")
  inter_file  <- parse_arg("--intersection-genes-file")
  force_sh    <- "--force-split-half" %in% args

  if (is.null(matrix_path) || is.null(meta_path) || is.null(domain) || is.null(outdir))
    stop("Required: --matrix, --meta, --domain, --outdir")
  if (is.null(contrast_a) || is.null(contrast_b))
    stop("Required: --contrast-a and --contrast-b")

  mat  <- as.matrix(fread(matrix_path), rownames = 1)
  meta <- fread(meta_path)

  inter_genes <- NULL
  if (!is.null(inter_file) && file.exists(inter_file)) {
    inter_genes <- trimws(readLines(inter_file, warn = FALSE))
    inter_genes <- inter_genes[nzchar(inter_genes)]
    cat("Intersection gene filter:", length(inter_genes), "genes\n")
  }

  cat("Computing concordance ceiling for", domain, "\n")
  res <- compute_concordance_ceiling(mat, meta, domain, contrast_a, contrast_b,
                                      n_splits = n_splits, seed = seed,
                                      intersection_genes = inter_genes,
                                      force_split_half = force_sh)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(res$splits, file.path(outdir, paste0("concordance_ceiling_", tolower(domain), ".csv")))
  fwrite(res$summary, file.path(outdir, paste0("ceiling_summary_", tolower(domain), ".csv")))
  cat("Saved to:", outdir, "\n")
}
