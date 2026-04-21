#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) dirname(normalizePath(sub("^--file=", "", ff[1]))) else normalizePath(file.path(getwd(), "scripts", "presentation"))
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

out_md <- file.path(PRES_OUT, "backup_slides", "backup_content.md")
sink(out_md)

cat("# Backup Slides Content\n\n")
cat("Generated:", format(Sys.time(), usetz = TRUE), "\n\n")
cat("Repository:", REPO, "\n\n")

cat("## Backup 1: Full Marker Panel (Subtype)\n\n")
pf <- file.path(PRES_OUT, "tables/marker_panel_subtype.csv")
if (file.exists(pf)) {
  panel <- fread(pf)
  cat("```\n")
  print(as.data.frame(panel))
  cat("```\n\n")
} else {
  cat("(Run extract_marker_panel.R first.)\n\n")
}

cat("## Backup 2: Fold Change and SE Summary\n\n")
ffcsv <- file.path(PRES_OUT, "tables/fc_se_summary.csv")
if (file.exists(ffcsv)) {
  fc_se <- fread(ffcsv)
  cat("```\n")
  print(as.data.frame(fc_se))
  cat("```\n\n")
} else {
  cat("(Run fc_se_summary.R first.)\n\n")
}

cat("## Backup 3: Bridge Normalization — Per-Protein Statistics\n\n")
for (fn in c(
  "tables/bridge_analysis_raw_breast_subtype.csv",
  "tables/bridge_analysis_bridge_shift_breast_subtype.csv"
)) {
  bf <- file.path(PRES_OUT, fn)
  if (file.exists(bf)) {
    cat("### ", basename(bf), "\n\n```\n", sep = "")
    print(as.data.frame(fread(bf)))
    cat("```\n\n")
  }
}

cat("**Bridge shift (concept):** per-gene location alignment using bridge summaries; within-domain contrasts preserved for pure shift.\n\n")

cat("## Backup 4: Residual Dependence\n\n")
for (method in c("raw", "bridge_shift", "celligner")) {
  for (task in c("breast_subtype")) {
    for (domain in c("cptac", "ccle")) {
      f <- file.path(REPO, sprintf(
        "reports/benchmark_master/benchmark_results/%s/%s/calibration/residual_dependence_%s.csv",
        method, task, domain
      ))
      if (file.exists(f)) {
        rd <- fread(f)
        cat(sprintf("### %s / %s / %s\n```\n", method, task, domain))
        print(as.data.frame(rd))
        cat("```\n\n")
      }
    }
  }
}

cat("## Backup 5: Biology Destruction Sensitivity (grid head)\n\n")
for (method in c("raw", "bridge_shift", "bridge_scale", "celligner")) {
  f <- file.path(REPO, sprintf(
    "reports/benchmark_master/benchmark_results/%s/breast_subtype/calibration/destruction_grid_%s.csv",
    method, method
  ))
  if (file.exists(f)) {
    grid <- fread(f)
    cat(sprintf("### %s\n```\n", method))
    print(head(as.data.frame(grid), 20L))
    cat("```\n\n")
  }
}

cat("## Backup 6: Stratified FC\n\n")
strat_files <- Sys.glob(file.path(REPO, "reports/benchmark_master/diagnostics", "fc_stratified*.csv"))
for (f in strat_files) {
  strat <- fread(f)
  cat(sprintf("### %s\n```\n", basename(f)))
  print(as.data.frame(strat))
  cat("```\n\n")
}

cat("## Backup 7: Permutation Null Summaries\n\n")
for (method in c("raw", "bridge_shift", "bridge_scale", "celligner")) {
  for (task in c("breast_subtype", "breast_vs_lung")) {
    f <- file.path(REPO, sprintf(
      "reports/benchmark_master/benchmark_results/%s/%s/calibration/observed_vs_null_summary.csv",
      method, task
    ))
    if (file.exists(f)) {
      cat(sprintf("### %s / %s\n```\n", method, task))
      print(as.data.frame(fread(f)))
      cat("```\n\n")
    }
  }
}

cat("## Backup 8: Gene Coverage Audit\n\n")
for (task in c("breast_subtype", "breast_vs_lung")) {
  f <- file.path(REPO, sprintf("reports/benchmark_master/diagnostics/gene_coverage_audit_%s.csv", task))
  if (file.exists(f)) {
    audit <- fread(f)
    cat(sprintf("### %s\n```\n", task))
    print(table(audit$category))
    cat("```\n\n")
  }
}

cat("## Backup 9: CCLE Annotation (processed)\n\n")
ann_file <- file.path(REPO, "data/processed/ccle_breast_subtype_annotation_processed.csv")
if (file.exists(ann_file)) {
  ann <- fread(ann_file)
  cols <- intersect(c("cell_line", "BvL_group", "subtype_detail", "plexes", "n_plexes"), names(ann))
  cat("```\n")
  print(as.data.frame(ann[, ..cols]))
  cat("```\n\n")
  if ("BvL_group" %in% names(ann)) {
    cat(sprintf("Luminal: %d lines\n", sum(ann$BvL_group == "Luminal", na.rm = TRUE)))
    cat(sprintf("Basal: %d lines\n", sum(ann$BvL_group == "Basal", na.rm = TRUE)))
  }
} else {
  cat("(Annotation file not found.)\n\n")
}

sink()
message("Wrote ", out_md)
