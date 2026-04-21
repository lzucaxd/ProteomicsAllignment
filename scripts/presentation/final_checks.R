#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1]))) else normalizePath(file.path(getwd(), "scripts", "presentation"))
source(file.path(sd, "presentation_paths.R"))

cat("\n")
cat("========================================================================\n")
cat("        PRESENTATION READINESS CHECKLIST\n")
cat("========================================================================\n\n")

checks <- list()

f <- file.path(REPO, "reports/benchmark_master/benchmark_results/comparison_summary.csv")
if (file.exists(f)) {
  s <- fread(f)
  checks[["comparison_summary"]] <- sprintf("PASS — %d rows (expect 8 = 4 methods x 2 tasks)", nrow(s))
} else {
  checks[["comparison_summary"]] <- "FAIL — comparison_summary.csv not found"
}

cell_da <- file.path(REPO, "reports/benchmark_master/benchmark_results/celligner/breast_subtype/representation_da/cptac/da_limma_result.csv")
raw_da <- file.path(REPO, "reports/benchmark_master/benchmark_results/raw/breast_subtype/representation_da/cptac/da_limma_result.csv")
if (file.exists(cell_da) && file.exists(raw_da)) {
  c_da <- fread(cell_da)
  r_da <- fread(raw_da)
  m <- merge(c_da, r_da, by = "gene", suffixes = c("_cell", "_raw"))
  cor_val <- cor(m$logFC_cell, m$logFC_raw, use = "complete.obs")
  if (!is.na(cor_val) && cor_val > 0.999) {
    checks[["celligner"]] <- sprintf("WARNING — logFC r=%.4f vs raw (possible scaffold)", cor_val)
  } else {
    checks[["celligner"]] <- sprintf("PASS — logFC r=%.4f vs raw", cor_val)
  }
} else {
  checks[["celligner"]] <- "FAIL — DA files missing"
}

null_file <- file.path(REPO, "reports/benchmark_master/benchmark_results/raw/breast_subtype/calibration/observed_vs_null_summary.csv")
checks[["permutation_nulls"]] <- if (file.exists(null_file)) {
  "PASS — null summary exists"
} else {
  "WARNING — observed_vs_null_summary missing"
}

panel_file <- file.path(PRES_OUT, "tables/marker_panel_subtype.csv")
if (file.exists(panel_file)) {
  panel <- fread(panel_file)
  n_found <- sum(!is.na(panel$cptac_logFC))
  checks[["markers"]] <- sprintf("PASS — %d/%d markers with CPTAC logFC", n_found, nrow(panel))
} else {
  checks[["markers"]] <- "NOT RUN — marker_panel_subtype.csv missing"
}

fig_dir <- file.path(REPO, "reports/benchmark_master/meeting/figures")
if (dir.exists(fig_dir)) {
  figs <- list.files(fig_dir, pattern = "\\.(png|pdf)$", ignore.case = TRUE)
  checks[["meeting_figures"]] <- sprintf("PASS — %d meeting figures", length(figs))
} else {
  checks[["meeting_figures"]] <- "FAIL — meeting/figures missing"
}

prof_count <- if (dir.exists(file.path(PRES_OUT, "figures"))) {
  length(list.files(file.path(PRES_OUT, "figures"), pattern = "^profile_.*\\.pdf$"))
} else 0L
checks[["profile_plots"]] <- sprintf("%d profile_*.pdf in presentation_materials/figures", prof_count)

fcse_file <- file.path(PRES_OUT, "tables/fc_se_summary.csv")
checks[["fc_se_summary"]] <- if (file.exists(fcse_file)) "PASS" else "NOT RUN — fc_se_summary.csv missing"

for (name in names(checks)) {
  status <- checks[[name]]
  icon <- if (grepl("^PASS", status)) "[ok]" else if (grepl("^WARNING", status)) "[!]" else if (grepl("^FAIL", status)) "[x]" else "[ ]"
  cat(sprintf("  %s  %-22s %s\n", icon, name, status))
}

cat("\n=== KEY NUMBERS (raw) ===\n\n")
if (file.exists(f)) {
  s <- fread(f)
  raw_sub <- s[method == "raw" & task == "breast_subtype"]
  raw_bvl <- s[method == "raw" & task == "breast_vs_lung"]

  if (nrow(raw_sub) == 1L) {
    cat("Subtype (raw):\n")
    cat(sprintf("  FC corr intersection: %s\n", raw_sub$fc_correlation_intersection))
    cat(sprintf("  Perm z: %s\n", raw_sub$permutation_z_fc_corr))
    cat(sprintf("  Ceiling: %s\n", raw_sub$concordance_ceiling_fc_corr))
    cat(sprintf("  Marker sanity CPTAC: %s\n", raw_sub$marker_sanity_cptac))
    cat(sprintf("  Marker sanity CCLE: %s\n", raw_sub$marker_sanity_ccle))
    cat(sprintf("  n_ccle_samples: %s\n", raw_sub$n_ccle_samples))
  }
  if (nrow(raw_bvl) == 1L) {
    cat("\nBreast vs lung (raw):\n")
    cat(sprintf("  FC corr intersection: %s\n", raw_bvl$fc_correlation_intersection))
    cat(sprintf("  Perm p: %s\n", raw_bvl$permutation_p_fc_corr))
    cat(sprintf("  Domain R2 PC1: %s\n", raw_bvl$struct_domain_r2_pc1))
  }
}

cat("\n")
