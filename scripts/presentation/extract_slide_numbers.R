#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})
ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1]))) else normalizePath(file.path(getwd(), "scripts", "presentation"))
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

summary_path <- file.path(REPO, "reports/benchmark_master/benchmark_results/comparison_summary.csv")
if (!file.exists(summary_path)) {
  stop("Missing ", summary_path, " — run the benchmark first.")
}
summary <- fread(summary_path)

cat("\n=== SLIDE 3: SAMPLE COUNTS ===\n")
for (i in seq_len(nrow(summary))) {
  row <- summary[i, ]
  cat(sprintf(
    "  %s / %s: struct_n_samples=%s, n_ccle_samples=%s\n",
    row$method, row$task,
    row$struct_n_samples, row$n_ccle_samples
  ))
}

cat("\n=== SLIDE 4: GENE COVERAGE ===\n")
for (task in unique(summary$task)) {
  audit_file <- file.path(REPO, sprintf(
    "reports/benchmark_master/diagnostics/gene_coverage_audit_%s.csv", task
  ))
  if (file.exists(audit_file)) {
    audit <- fread(audit_file)
    cat(sprintf("  %s:\n", task))
    print(table(audit$category))
    cat(sprintf(
      "  Total: %d, both_domains: %d\n\n",
      nrow(audit), sum(audit$category == "both_domains", na.rm = TRUE)
    ))
  } else {
    cat(sprintf("  %s: audit file missing\n", task))
  }
}

cat("\n=== SLIDE 11: DA GENE COUNTS + EFFECT SIZES (raw only) ===\n")
cat("(Canonical DA: representation_da/<cptac|ccle>/da_limma_result.csv; logFC = Luminal-Basal or Lung-Breast.)\n")
for (task in unique(summary$task)) {
  for (domain in c("cptac", "ccle")) {
    da_file <- file.path(REPO, sprintf(
      "reports/benchmark_master/benchmark_results/raw/%s/representation_da/%s/da_limma_result.csv",
      task, domain
    ))
    if (file.exists(da_file)) {
      da <- fread(da_file)
      sig <- da[!is.na(adj.P.Val) & adj.P.Val < 0.05]
      cat(sprintf("  %s / %s / raw:\n", task, toupper(domain)))
      cat(sprintf("    Total genes tested: %d\n", nrow(da)))
      cat(sprintf("    Significant (FDR<0.05): %d\n", nrow(sig)))
      cat(sprintf("    Up (logFC>0): %d\n", sum(sig$logFC > 0, na.rm = TRUE)))
      cat(sprintf("    Down (logFC<0): %d\n", sum(sig$logFC < 0, na.rm = TRUE)))
      cat(sprintf(
        "    Median |logFC| (sig genes): %.3f\n",
        median(abs(sig$logFC), na.rm = TRUE)
      ))
      if ("t" %in% names(sig) && any(!is.na(sig$t) & sig$t != 0)) {
        se <- abs(sig$logFC / sig$t)
        cat(sprintf("    Median SE approx (sig genes): %.4f\n", median(se, na.rm = TRUE)))
      }
      cat("\n")
    }
  }
}

cat("\n=== SLIDE 12: BREAST SUBTYPE CROSS-DOMAIN ===\n")
sub <- summary[task == "breast_subtype"]
cols_to_show <- c(
  "method",
  "fc_correlation_intersection", "fc_same_dir_intersection",
  "n_genes_intersection",
  "permutation_z_fc_corr", "permutation_p_fc_corr",
  "concordance_ceiling_fc_corr", "calibrated_fc_corr_intersection",
  "marker_sanity_cptac", "marker_sanity_ccle",
  "biology_destruction_retention", "biology_destruction_fc_shrinkage"
)
cols_present <- intersect(cols_to_show, names(sub))
print(as.data.frame(sub[, ..cols_present]))

cat("\n=== SLIDES 13-14: BREAST VS LUNG ===\n")
bvl <- summary[task == "breast_vs_lung"]
cat("\nGeometry metrics (subset):\n")
geom_cols <- grep("struct_|classification|silhouette|knn", names(bvl), value = TRUE)
geom_cols <- c("method", intersect(geom_cols, names(bvl)))
if (length(geom_cols) > 1L) print(as.data.frame(bvl[, ..geom_cols]))

cat("\nDA metrics:\n")
da_cols <- c(
  "method",
  "fc_correlation_intersection", "fc_correlation_union",
  "fc_same_dir_intersection", "fc_same_dir_union",
  "permutation_z_fc_corr", "permutation_p_fc_corr",
  "biology_destruction_retention", "biology_destruction_fc_shrinkage"
)
da_present <- intersect(da_cols, names(bvl))
print(as.data.frame(bvl[, ..da_present]))

cat("\n=== SLIDE 16: DISCONNECT SCORES ===\n")
disc_file <- file.path(REPO, "reports/benchmark_master/benchmark_results/disconnect_scores.csv")
if (file.exists(disc_file)) {
  print(fread(disc_file))
} else {
  cat("(disconnect_scores.csv not found)\n")
}

fwrite(summary, file.path(PRES_OUT, "tables/comparison_summary_full.csv"))
cat("\nWrote ", file.path(PRES_OUT, "tables/comparison_summary_full.csv"), "\n", sep = "")
