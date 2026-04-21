#!/usr/bin/env Rscript
# Step 4: Union + intersection cross-domain FC metrics; writes fc_agreement.csv + cross_domain_metrics.csv

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
results_root <- "reports/benchmark_master/benchmark_results"
inter_dir <- "data/processed"
for (i in seq_along(args)) {
  if (args[i] == "--repo-root" && i < length(args)) repo_root <- args[i + 1]
  if (args[i] == "--results-root" && i < length(args)) results_root <- args[i + 1]
  if (args[i] == "--intersection-dir" && i < length(args)) inter_dir <- args[i + 1]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
tasks <- c("breast_subtype", "breast_vs_lung")

compute_block <- function(da_c, da_e, inter_genes) {
  merged <- merge(
    da_c[, .(gene, logFC_cptac = logFC, adj.P.Val_cptac = adj.P.Val)],
    da_e[, .(gene, logFC_ccle = logFC, adj.P.Val_ccle = adj.P.Val)],
    by = "gene", all = FALSE
  )
  merged <- merged[is.finite(logFC_cptac) & is.finite(logFC_ccle)]
  merged[, same_direction := sign(logFC_cptac) == sign(logFC_ccle)]
  merged[, fc_diff := logFC_cptac - logFC_ccle]

  met <- function(df, label) {
    if (nrow(df) < 3) {
      return(data.table(
        gene_set = label, n_genes = nrow(df),
        fc_correlation = NA_real_, same_dir_fraction = NA_real_,
        median_abs_fc_diff = NA_real_
      ))
    }
    data.table(
      gene_set = label,
      n_genes = nrow(df),
      fc_correlation = cor(df$logFC_cptac, df$logFC_ccle, use = "complete.obs"),
      same_dir_fraction = mean(df$same_direction, na.rm = TRUE),
      median_abs_fc_diff = median(abs(df$logFC_cptac - df$logFC_ccle), na.rm = TRUE)
    )
  }

  u <- met(merged, "union")
  sub <- merged[gene %in% inter_genes]
  i <- met(sub, "intersection")
  list(merged = merged, metrics = rbind(u, i))
}

all_rows <- list()

for (m in methods) {
  for (task in tasks) {
    inter_file <- file.path(repo_root, inter_dir, paste0("intersection_genes_", task, ".txt"))
    if (!file.exists(inter_file))
      inter_file <- file.path(repo_root, "data/processed_union", paste0("intersection_genes_", task, ".txt"))
    inter_genes <- if (file.exists(inter_file)) {
      trimws(readLines(inter_file, warn = FALSE))
    } else {
      character(0)
    }
    inter_genes <- inter_genes[nzchar(inter_genes)]

    base <- file.path(repo_root, results_root, m, task, "representation_da")
    f1 <- file.path(base, "cptac", "da_limma_result.csv")
    f2 <- file.path(base, "ccle", "da_limma_result.csv")
    if (!file.exists(f1) || !file.exists(f2)) {
      warning("Missing DA for ", m, " ", task)
      next
    }
    da_c <- fread(f1)
    da_e <- fread(f2)
    blk <- compute_block(da_c, da_e, inter_genes)

    fwrite(blk$merged, file.path(base, "fc_agreement.csv"))
    met <- copy(blk$metrics)
    met[, `:=`(method = m, task = task)]
    rest <- setdiff(names(met), c("method", "task"))
    setcolorder(met, c("method", "task", rest))
    fwrite(met, file.path(base, "cross_domain_metrics.csv"))

    cat("\n=== Cross-domain FC:", m, task, "===\n")
    print(met)

    all_rows[[length(all_rows) + 1]] <- met
  }
}

if (length(all_rows)) {
  out_all <- file.path(repo_root, results_root, "cross_domain_metrics_all.csv")
  fwrite(rbindlist(all_rows), out_all)
  cat("\nWrote combined:", out_all, "\n")
}
cat("Step 4 complete.\n")
