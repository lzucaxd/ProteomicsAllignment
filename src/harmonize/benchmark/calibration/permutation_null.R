#!/usr/bin/env Rscript
# =============================================================================
# Permutation Null Calibration
# =============================================================================
# Generates null distributions for cross-domain agreement metrics by shuffling
# condition labels *within* each domain independently.
#
# Usage (CLI):
#   Rscript permutation_null.R \
#     --matrix input.csv --meta meta.csv \
#     --contrast-a Basal --contrast-b Luminal \
#     --n-perm 200 --seed 42 \
#     --markers "ESR1,PGR,GATA3,FOXA1,EGFR,KRT5,KRT17,FOXC1" \
#     --expected-signs "1,1,1,1,-1,-1,-1,-1" \
#     --outdir /path/to/output
#
# Or source() and call compute_permutation_null() directly.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

compute_permutation_null <- function(matrix, sample_meta, contrast_a, contrast_b,
                                      n_perm = 200, seed = 42,
                                      marker_genes = NULL, expected_directions = NULL,
                                      intersection_genes = NULL) {
  # ---------------------------------------------------------------------------
  # Shuffle condition labels within each domain, run limma, compute cross-domain
  # agreement metrics. Returns list($null_distribution, $observed, $summary).
  # ---------------------------------------------------------------------------
  set.seed(seed)

  meta <- as.data.table(sample_meta)
  stopifnot(all(c("sample_id", "domain", "condition") %in% names(meta)))

  # Only keep samples in the matrix and with the two contrast levels
  meta <- meta[sample_id %in% colnames(matrix) &
               condition %in% c(contrast_a, contrast_b)]

  # Observed DA per domain
  run_domain_da <- function(mat, meta_dt, domain, ca, cb) {
    dom_sids <- meta_dt[toupper(domain_col) == toupper(domain), sample_id]
    dom_cols <- intersect(dom_sids, colnames(mat))
    if (length(dom_cols) < 4) return(NULL)

    dom_mat    <- mat[, dom_cols, drop = FALSE]
    dom_groups <- meta_dt[match(dom_cols, sample_id), condition]
    if (!all(c(ca, cb) %in% unique(dom_groups))) return(NULL)

    # Filter genes
    for (g in c(ca, cb)) {
      g_cols <- which(dom_groups == g)
      if (length(g_cols) > 0) {
        na_frac <- rowMeans(is.na(dom_mat[, g_cols, drop = FALSE]))
        dom_mat[na_frac > 0.5, g_cols] <- NA
      }
    }
    keep <- rowSums(!is.na(dom_mat)) >= ncol(dom_mat) * 0.5
    dom_mat <- dom_mat[keep, , drop = FALSE]
    if (nrow(dom_mat) < 10) return(NULL)

    gf <- factor(dom_groups, levels = c(ca, cb))
    design <- model.matrix(~ 0 + gf)
    colnames(design) <- levels(gf)
    fit <- lmFit(dom_mat, design)
    cm <- makeContrasts(contrasts = paste0(cb, " - ", ca), levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cm))
    tt <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"))
    tt$gene <- rownames(topTable(fit2, number = Inf, sort.by = "none"))
    tt
  }

  compute_agreement_metrics <- function(da_cptac, da_ccle, mk_genes, exp_dirs, gene_keep = NULL) {
    shared <- merge(da_cptac[, .(gene, fc_cptac = logFC)],
                    da_ccle[, .(gene, fc_ccle = logFC)], by = "gene")
    shared <- shared[is.finite(fc_cptac) & is.finite(fc_ccle)]
    if (!is.null(gene_keep) && length(gene_keep) > 0)
      shared <- shared[gene %in% gene_keep]
    if (nrow(shared) < 3) return(NULL)

    fc_corr <- cor(shared$fc_cptac, shared$fc_ccle, method = "pearson")
    same_dir <- mean(sign(shared$fc_cptac) == sign(shared$fc_ccle))
    med_fc_diff <- median(abs(shared$fc_cptac - shared$fc_ccle))

    marker_conc <- NA_real_
    if (!is.null(mk_genes) && !is.null(exp_dirs)) {
      mk_data <- merge(shared, exp_dirs, by = "gene")
      if (nrow(mk_data) > 0) {
        cptac_ok <- mean(sign(mk_data$fc_cptac) == mk_data$expected_sign)
        ccle_ok  <- mean(sign(mk_data$fc_ccle) == mk_data$expected_sign)
        marker_conc <- (cptac_ok + ccle_ok) / 2
      }
    }

    data.table(fc_correlation = fc_corr, same_direction_frac = same_dir,
               median_abs_fc_diff = med_fc_diff, marker_concordance = marker_conc)
  }

  setnames(meta, "domain", "domain_col")

  # Observed
  da_cptac_obs <- run_domain_da(matrix, meta, "CPTAC", contrast_a, contrast_b)
  da_ccle_obs  <- run_domain_da(matrix, meta, "CCLE", contrast_a, contrast_b)

  if (is.null(da_cptac_obs) || is.null(da_ccle_obs))
    stop("Cannot compute observed agreement: insufficient data in one or both domains")

  observed <- compute_agreement_metrics(da_cptac_obs, da_ccle_obs,
                                         marker_genes, expected_directions,
                                         gene_keep = intersection_genes)

  # Permutation null
  null_rows <- vector("list", n_perm)

  for (i in seq_len(n_perm)) {
    meta_perm <- copy(meta)

    # Shuffle condition within each domain independently
    for (dom in unique(meta_perm$domain_col)) {
      idx <- which(meta_perm$domain_col == dom)
      meta_perm$condition[idx] <- sample(meta_perm$condition[idx])
    }

    da_cptac_perm <- run_domain_da(matrix, meta_perm, "CPTAC", contrast_a, contrast_b)
    da_ccle_perm  <- run_domain_da(matrix, meta_perm, "CCLE", contrast_a, contrast_b)

    if (!is.null(da_cptac_perm) && !is.null(da_ccle_perm)) {
      null_rows[[i]] <- compute_agreement_metrics(da_cptac_perm, da_ccle_perm,
                                                    marker_genes, expected_directions,
                                                    gene_keep = intersection_genes)
      null_rows[[i]]$perm <- i
    }

    if (i %% 50 == 0) cat("  Permutation", i, "/", n_perm, "\n")
  }

  null_dist <- rbindlist(null_rows[!sapply(null_rows, is.null)])

  # Summary: observed vs null
  summary_rows <- list()
  for (metric in c("fc_correlation", "same_direction_frac",
                    "median_abs_fc_diff", "marker_concordance")) {
    obs_val  <- observed[[metric]]
    null_vals <- null_dist[[metric]]
    null_vals <- null_vals[is.finite(null_vals)]

    if (length(null_vals) < 5 || !is.finite(obs_val)) next

    null_mean <- mean(null_vals)
    null_sd   <- sd(null_vals)
    z_score   <- if (null_sd > 0) (obs_val - null_mean) / null_sd else NA_real_

    if (metric == "median_abs_fc_diff") {
      p_val <- mean(null_vals <= obs_val)
    } else {
      p_val <- mean(null_vals >= obs_val)
    }

    summary_rows[[metric]] <- data.table(
      metric = metric, observed = obs_val,
      null_mean = null_mean, null_sd = null_sd,
      z_score = z_score, p_value = p_val, n_perm = length(null_vals)
    )
  }

  summary_dt <- rbindlist(summary_rows)

  list(null_distribution = null_dist, observed = observed, summary = summary_dt)
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
  contrast_a  <- parse_arg("--contrast-a")
  contrast_b  <- parse_arg("--contrast-b")
  outdir      <- parse_arg("--outdir")
  n_perm      <- as.integer(parse_arg("--n-perm") %||% "200")
  seed        <- as.integer(parse_arg("--seed") %||% "42")
  mk_str      <- parse_arg("--markers")
  es_str      <- parse_arg("--expected-signs")
  inter_file  <- parse_arg("--intersection-genes-file")

  if (is.null(matrix_path) || is.null(meta_path) || is.null(outdir))
    stop("Required: --matrix, --meta, --outdir")
  if (is.null(contrast_a) || is.null(contrast_b))
    stop("Required: --contrast-a and --contrast-b")

  mat  <- as.matrix(fread(matrix_path), rownames = 1)
  meta <- fread(meta_path)

  mk_genes <- NULL
  exp_dirs <- NULL
  if (!is.null(mk_str) && !is.null(es_str)) {
    mk_genes <- trimws(strsplit(mk_str, ",")[[1]])
    exp_signs <- as.integer(trimws(strsplit(es_str, ",")[[1]]))
    exp_dirs <- data.table(gene = mk_genes, expected_sign = exp_signs)
  }

  inter_genes <- NULL
  if (!is.null(inter_file) && file.exists(inter_file)) {
    inter_genes <- trimws(readLines(inter_file, warn = FALSE))
    inter_genes <- inter_genes[nzchar(inter_genes)]
    cat("Intersection gene filter:", length(inter_genes), "genes\n")
  }

  cat("Running permutation null calibration (", n_perm, " permutations)\n")

  res <- compute_permutation_null(mat, meta, contrast_a, contrast_b,
                                   n_perm = n_perm, seed = seed,
                                   marker_genes = mk_genes,
                                   expected_directions = exp_dirs,
                                   intersection_genes = inter_genes)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(res$null_distribution, file.path(outdir, "null_distribution.csv"))
  fwrite(res$observed, file.path(outdir, "observed_metrics.csv"))
  fwrite(res$summary, file.path(outdir, "observed_vs_null_summary.csv"))
  cat("Saved to:", outdir, "\n")
}
