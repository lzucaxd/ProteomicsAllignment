#!/usr/bin/env Rscript
# Extract all presentation numbers into presentation_materials/ALL_NUMBERS.md
# Reads existing outputs only; does not run new analyses.

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) {
  dirname(normalizePath(sub("^--file=", "", ff[1L])))
} else {
  normalizePath(file.path(getwd(), "scripts", "presentation"), mustWork = FALSE)
}
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

OUT_MD <- file.path(PRES_OUT, "ALL_NUMBERS.md")

methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
tasks <- c("breast_subtype", "breast_vs_lung")
domains <- c("cptac", "ccle")

BENCH <- file.path(REPO, "reports", "benchmark_master", "benchmark_results")

safe_fread <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(fread(path), error = function(e) NULL)
}

da_path <- function(method, task, domain) {
  p1 <- file.path(BENCH, method, task, "representation_da", domain, "da_limma_result.csv")
  if (file.exists(p1)) return(p1)
  p2 <- file.path(BENCH, method, task, "representation_da", sprintf("da_%s.csv", domain))
  if (file.exists(p2)) return(p2)
  NA_character_
}

calib_path <- function(method, task, fname) {
  file.path(BENCH, method, task, "calibration", fname)
}

fmt3 <- function(x) ifelse(is.finite(as.numeric(x)), sprintf("%.3f", as.numeric(x)), "NA")
fmt2 <- function(x) ifelse(is.finite(as.numeric(x)), sprintf("%.2f", as.numeric(x)), "NA")
fmt1pct <- function(x) ifelse(is.finite(as.numeric(x)), sprintf("%.1f%%", 100 * as.numeric(x)), "NA")
fmt3s <- function(x) ifelse(is.finite(as.numeric(x)), sprintf("%+.3f", as.numeric(x)), "NA")
fmt4 <- function(x) ifelse(is.finite(as.numeric(x)), sprintf("%.4f", as.numeric(x)), "NA")

counts_by <- function(meta, domain_val, cond_vals) {
  dm <- meta[toupper(domain) == toupper(domain_val)]
  out <- setNames(integer(0), character(0))
  for (cv in cond_vals) out[[cv]] <- sum(dm$condition == cv, na.rm = TRUE)
  out
}

get_marker_row <- function(da, gene) {
  if (is.null(da) || nrow(da) == 0) return(NULL)
  gcol <- intersect(c("gene", "Gene"), names(da))
  if (length(gcol) == 0) return(NULL)
  gcol <- gcol[1]
  # data.table scoping: avoid collision with a column named "gene"
  gene_q <- gene
  row <- da[get(gcol) == gene_q]
  if (nrow(row) == 0) return(NULL)
  row[1]
}

marker_val <- function(row) {
  if (is.null(row)) return("not quantified")
  fccol <- intersect(c("logFC", "log2FC"), names(row))
  pcol <- intersect(c("adj.P.Val", "padj", "FDR"), names(row))
  if (length(fccol) == 0 || length(pcol) == 0) return("not quantified")
  fc <- as.numeric(row[[fccol[1]]])
  pv <- as.numeric(row[[pcol[1]]])
  sig <- is.finite(pv) && pv < 0.05
  sprintf("%s (%s)", fmt3s(fc), ifelse(sig, "sig", "no"))
}

write_slide_counts <- function() {
  meta_sub <- safe_fread(file.path(REPO, "data", "processed", "union", "sample_meta_breast_subtype.csv"))
  meta_bvl <- safe_fread(file.path(REPO, "data", "processed", "union", "sample_meta_breast_vs_lung.csv"))

  cat("SLIDE 3 — SAMPLE COUNTS\n\n")
  if (!is.null(meta_sub)) {
    cat("Breast Subtype:\n")
    cpt <- counts_by(meta_sub, "CPTAC", c("Luminal", "Basal"))
    ccl <- counts_by(meta_sub, "CCLE", c("Luminal", "Basal"))
    cat(sprintf("  CPTAC: %d Luminal + %d Basal = %d total\n", cpt[["Luminal"]], cpt[["Basal"]], sum(cpt)))
    cat(sprintf("  CCLE:  %d Luminal + %d Basal = %d total\n", ccl[["Luminal"]], ccl[["Basal"]], sum(ccl)))
    cat(sprintf("  Combined: %d\n\n", nrow(meta_sub)))
  }

  if (!is.null(meta_bvl)) {
    cat("Breast vs Lung:\n")
    cpt <- counts_by(meta_bvl, "CPTAC", c("Breast", "Lung"))
    ccl <- counts_by(meta_bvl, "CCLE", c("Breast", "Lung"))
    cat(sprintf("  CPTAC: %d Breast + %d Lung = %d total\n", cpt[["Breast"]], cpt[["Lung"]], sum(cpt)))
    cat(sprintf("  CCLE:  %d Breast + %d Lung = %d total\n", ccl[["Breast"]], ccl[["Lung"]], sum(ccl)))
    cat(sprintf("  Combined: %d\n\n", nrow(meta_bvl)))
  }
}

write_gene_coverage <- function() {
  cat("SLIDE 5 — GENE COVERAGE\n\n")
  for (task in c("breast_subtype", "breast_vs_lung")) {
    audit <- safe_fread(file.path(REPO, "reports", "benchmark_master", "diagnostics",
                                 sprintf("gene_coverage_audit_%s.csv", task)))
    if (is.null(audit) || !"category" %in% names(audit)) next
    both <- sum(audit$category == "both_domains", na.rm = TRUE)
    cat(sprintf("%s: %d intersection genes\n", ifelse(task == "breast_subtype", "Subtype", "BvL"), both))
  }
  cat("\n")
}

write_raw_da_counts <- function() {
  cat("SLIDE 9 — DA GENE COUNTS (raw method)\n\n")
  cat("| Metric | CPTAC subtype | CCLE subtype | CPTAC BvL | CCLE BvL |\n")
  cat("|---|---:|---:|---:|---:|\n")

  cols <- list(
    `CPTAC subtype` = list(task = "breast_subtype", domain = "cptac"),
    `CCLE subtype`  = list(task = "breast_subtype", domain = "ccle"),
    `CPTAC BvL`     = list(task = "breast_vs_lung", domain = "cptac"),
    `CCLE BvL`      = list(task = "breast_vs_lung", domain = "ccle")
  )

  metrics <- list()
  for (nm in names(cols)) {
    p <- da_path("raw", cols[[nm]]$task, cols[[nm]]$domain)
    da <- safe_fread(p)
    if (is.null(da)) {
      metrics[[nm]] <- list(n = NA, nsig = NA, up = NA, down = NA, med_abs = NA, med_se = NA)
      next
    }
    pcol <- intersect(c("adj.P.Val", "padj", "FDR"), names(da))
    fccol <- intersect(c("logFC", "log2FC"), names(da))
    tcol <- intersect(c("t", "t_stat"), names(da))
    if (length(pcol) == 0 || length(fccol) == 0) {
      metrics[[nm]] <- list(n = nrow(da), nsig = NA, up = NA, down = NA, med_abs = NA, med_se = NA)
      next
    }
    pcol <- pcol[1]; fccol <- fccol[1]
    sig <- da[is.finite(get(pcol)) & get(pcol) < 0.05]
    med_abs <- median(abs(sig[[fccol]]), na.rm = TRUE)
    med_se <- NA_real_
    if (length(tcol) > 0) {
      tcol <- tcol[1]
      med_se <- median(abs(sig[[fccol]] / sig[[tcol]]), na.rm = TRUE)
    }
    metrics[[nm]] <- list(
      n = nrow(da),
      nsig = nrow(sig),
      up = sum(sig[[fccol]] > 0, na.rm = TRUE),
      down = sum(sig[[fccol]] < 0, na.rm = TRUE),
      med_abs = med_abs,
      med_se = med_se
    )
  }

  rowv <- function(field, fmt_fun = identity) {
    c(
      fmt_fun(metrics[["CPTAC subtype"]][[field]]),
      fmt_fun(metrics[["CCLE subtype"]][[field]]),
      fmt_fun(metrics[["CPTAC BvL"]][[field]]),
      fmt_fun(metrics[["CCLE BvL"]][[field]])
    )
  }
  cat(sprintf("| Genes tested | %s |\n", paste(rowv("n", function(x) ifelse(is.na(x), "NA", as.character(x))), collapse = " | ")))
  cat(sprintf("| Significant | %s |\n", paste(rowv("nsig", function(x) ifelse(is.na(x), "NA", as.character(x))), collapse = " | ")))
  cat(sprintf("| Up (logFC>0) | %s |\n", paste(rowv("up", function(x) ifelse(is.na(x), "NA", as.character(x))), collapse = " | ")))
  cat(sprintf("| Down (logFC<0) | %s |\n", paste(rowv("down", function(x) ifelse(is.na(x), "NA", as.character(x))), collapse = " | ")))
  cat(sprintf("| Median \\|logFC\\| sig | %s |\n", paste(rowv("med_abs", fmt3), collapse = " | ")))
  cat(sprintf("| Median SE sig | %s |\n\n", paste(rowv("med_se", fmt4), collapse = " | ")))
}

write_subtype_cross_domain <- function() {
  cs <- safe_fread(file.path(BENCH, "comparison_summary.csv"))
  if (is.null(cs)) return(invisible(NULL))
  sub <- cs[task == "breast_subtype"]
  if (nrow(sub) == 0) return(invisible(NULL))

  # map method names to slide display names
  disp <- function(m) {
    if (m == "raw") return("Raw")
    if (m == "bridge_shift") return("Dom.Shift")
    if (m == "bridge_scale") return("Dom.Scale")
    if (m == "celligner") return("Celligner")
    m
  }

  pick <- function(m, col) {
    r <- sub[method == m]
    if (nrow(r) == 0 || !col %in% names(r)) return(NA_character_)
    as.character(r[[col]][1])
  }

  cat("SLIDE 10 — SUBTYPE CROSS-DOMAIN\n\n")
  cat("| Metric | Raw | Dom.Shift | Dom.Scale | Celligner |\n")
  cat("|---|---:|---:|---:|---:|\n")
  cat(sprintf("| FC corr (int) | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "fc_correlation_intersection")),
              fmt3(pick("bridge_shift", "fc_correlation_intersection")),
              fmt3(pick("bridge_scale", "fc_correlation_intersection")),
              fmt3(pick("celligner", "fc_correlation_intersection"))))
  cat(sprintf("| Same-dir %% | %s | %s | %s | %s |\n",
              fmt1pct(pick("raw", "fc_same_dir_intersection")),
              fmt1pct(pick("bridge_shift", "fc_same_dir_intersection")),
              fmt1pct(pick("bridge_scale", "fc_same_dir_intersection")),
              fmt1pct(pick("celligner", "fc_same_dir_intersection"))))
  cat(sprintf("| n genes | %s | %s | %s | %s |\n",
              pick("raw", "n_genes_intersection"),
              pick("bridge_shift", "n_genes_intersection"),
              pick("bridge_scale", "n_genes_intersection"),
              pick("celligner", "n_genes_intersection")))
  cat(sprintf("| Perm z | %s | %s | %s | %s |\n",
              fmt2(pick("raw", "permutation_z_fc_corr")),
              fmt2(pick("bridge_shift", "permutation_z_fc_corr")),
              fmt2(pick("bridge_scale", "permutation_z_fc_corr")),
              fmt2(pick("celligner", "permutation_z_fc_corr"))))
  cat(sprintf("| Perm p | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "permutation_p_fc_corr")),
              fmt3(pick("bridge_shift", "permutation_p_fc_corr")),
              fmt3(pick("bridge_scale", "permutation_p_fc_corr")),
              fmt3(pick("celligner", "permutation_p_fc_corr"))))
  cat(sprintf("| Ceiling | %s | %s | %s | %s |\n",
              fmt2(pick("raw", "concordance_ceiling_fc_corr")),
              fmt2(pick("bridge_shift", "concordance_ceiling_fc_corr")),
              fmt2(pick("bridge_scale", "concordance_ceiling_fc_corr")),
              fmt2(pick("celligner", "concordance_ceiling_fc_corr"))))
  cat(sprintf("| Calibrated | %s | %s | %s | %s |\n",
              fmt1pct(pick("raw", "calibrated_fc_corr_intersection")),
              fmt1pct(pick("bridge_shift", "calibrated_fc_corr_intersection")),
              fmt1pct(pick("bridge_scale", "calibrated_fc_corr_intersection")),
              fmt1pct(pick("celligner", "calibrated_fc_corr_intersection"))))
  cat(sprintf("| Marker CPTAC | %s | %s | %s | %s |\n",
              fmt1pct(pick("raw", "marker_sanity_cptac")),
              fmt1pct(pick("bridge_shift", "marker_sanity_cptac")),
              fmt1pct(pick("bridge_scale", "marker_sanity_cptac")),
              fmt1pct(pick("celligner", "marker_sanity_cptac"))))
  cat(sprintf("| Marker CCLE | %s | %s | %s | %s |\n",
              fmt1pct(pick("raw", "marker_sanity_ccle")),
              fmt1pct(pick("bridge_shift", "marker_sanity_ccle")),
              fmt1pct(pick("bridge_scale", "marker_sanity_ccle")),
              fmt1pct(pick("celligner", "marker_sanity_ccle"))))

  # destruction columns may be NA for raw
  ret <- function(m) pick(m, "biology_destruction_retention")
  shr <- function(m) pick(m, "biology_destruction_fc_shrinkage")
  cat(sprintf("| Gene retention | — | %s | %s | %s |\n",
              fmt1pct(ret("bridge_shift")), fmt1pct(ret("bridge_scale")), fmt1pct(ret("celligner"))))
  cat(sprintf("| FC shrinkage | — | %s | %s | %s |\n\n",
              fmt1pct(shr("bridge_shift")), fmt1pct(shr("bridge_scale")), fmt1pct(shr("celligner"))))
}

write_marker_table <- function(task, markers, expected_labels, expected_sign) {
  cat(sprintf("SLIDE %s — %s MARKERS\n\n",
              ifelse(task == "breast_subtype", "11", "13"),
              ifelse(task == "breast_subtype", "SUBTYPE", "BvL")))

  # header
  cols <- c("Gene", "Expected",
            "Raw-CPTAC", "Raw-CCLE",
            "Shift-CPTAC", "Shift-CCLE",
            "Scale-CPTAC", "Scale-CCLE",
            "Cell-CPTAC", "Cell-CCLE")
  cat(paste(cols, collapse = " | "), "\n")
  cat(paste(rep("---", length(cols)), collapse = " | "), "\n")

  for (gi in seq_along(markers)) {
    g <- markers[[gi]]
    exp <- expected_labels[[gi]]
    rowvals <- c(g, exp)
    for (m in c("raw", "bridge_shift", "bridge_scale", "celligner")) {
      for (d in c("cptac", "ccle")) {
        p <- da_path(m, task, d)
        da <- safe_fread(p)
        row <- get_marker_row(da, g)
        rowvals <- c(rowvals, marker_val(row))
      }
    }
    cat(paste(rowvals, collapse = " | "), "\n")
  }
  cat("\n")
}

write_bvl_geometry <- function() {
  cs <- safe_fread(file.path(BENCH, "comparison_summary.csv"))
  if (is.null(cs)) return(invisible(NULL))
  bvl <- cs[task == "breast_vs_lung"]
  if (nrow(bvl) == 0) return(invisible(NULL))

  pick <- function(m, col) {
    r <- bvl[method == m]
    if (nrow(r) == 0 || !col %in% names(r)) return(NA_character_)
    as.character(r[[col]][1])
  }

  cat("SLIDE 12 — BvL GEOMETRY\n\n")
  cat("| Metric | Raw | Dom.Shift | Dom.Scale | Celligner |\n")
  cat("|---|---:|---:|---:|---:|\n")
  cat(sprintf("| Domain R² PC1 | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "struct_domain_r2_pc1")),
              fmt3(pick("bridge_shift", "struct_domain_r2_pc1")),
              fmt3(pick("bridge_scale", "struct_domain_r2_pc1")),
              fmt3(pick("celligner", "struct_domain_r2_pc1"))))
  cat(sprintf("| Domain silhouette | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "struct_silhouette_domain")),
              fmt3(pick("bridge_shift", "struct_silhouette_domain")),
              fmt3(pick("bridge_scale", "struct_silhouette_domain")),
              fmt3(pick("celligner", "struct_silhouette_domain"))))
  cat(sprintf("| kNN purity (domain) | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "struct_knn_purity_domain")),
              fmt3(pick("bridge_shift", "struct_knn_purity_domain")),
              fmt3(pick("bridge_scale", "struct_knn_purity_domain")),
              fmt3(pick("celligner", "struct_knn_purity_domain"))))
  cat(sprintf("| Condition R² PC1 | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "struct_condition_r2_pc1")),
              fmt3(pick("bridge_shift", "struct_condition_r2_pc1")),
              fmt3(pick("bridge_scale", "struct_condition_r2_pc1")),
              fmt3(pick("celligner", "struct_condition_r2_pc1"))))
  cat(sprintf("| Condition silhouette | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "struct_silhouette_condition")),
              fmt3(pick("bridge_shift", "struct_silhouette_condition")),
              fmt3(pick("bridge_scale", "struct_silhouette_condition")),
              fmt3(pick("celligner", "struct_silhouette_condition"))))
  cat(sprintf("| kNN purity (condition) | %s | %s | %s | %s |\n\n",
              fmt3(pick("raw", "struct_knn_purity_condition")),
              fmt3(pick("bridge_shift", "struct_knn_purity_condition")),
              fmt3(pick("bridge_scale", "struct_knn_purity_condition")),
              fmt3(pick("celligner", "struct_knn_purity_condition"))))
}

write_bvl_da <- function() {
  cs <- safe_fread(file.path(BENCH, "comparison_summary.csv"))
  if (is.null(cs)) return(invisible(NULL))
  bvl <- cs[task == "breast_vs_lung"]
  if (nrow(bvl) == 0) return(invisible(NULL))

  pick <- function(m, col) {
    r <- bvl[method == m]
    if (nrow(r) == 0 || !col %in% names(r)) return(NA_character_)
    as.character(r[[col]][1])
  }

  cat("SLIDE 14 — BvL DA\n\n")
  cat("| Metric | Raw | Dom.Shift | Dom.Scale | Celligner |\n")
  cat("|---|---:|---:|---:|---:|\n")
  cat(sprintf("| FC corr (int) | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "fc_correlation_intersection")),
              fmt3(pick("bridge_shift", "fc_correlation_intersection")),
              fmt3(pick("bridge_scale", "fc_correlation_intersection")),
              fmt3(pick("celligner", "fc_correlation_intersection"))))
  cat(sprintf("| Same-dir %% | %s | %s | %s | %s |\n",
              fmt1pct(pick("raw", "fc_same_dir_intersection")),
              fmt1pct(pick("bridge_shift", "fc_same_dir_intersection")),
              fmt1pct(pick("bridge_scale", "fc_same_dir_intersection")),
              fmt1pct(pick("celligner", "fc_same_dir_intersection"))))
  cat(sprintf("| Perm z | %s | %s | %s | %s |\n",
              fmt2(pick("raw", "permutation_z_fc_corr")),
              fmt2(pick("bridge_shift", "permutation_z_fc_corr")),
              fmt2(pick("bridge_scale", "permutation_z_fc_corr")),
              fmt2(pick("celligner", "permutation_z_fc_corr"))))
  cat(sprintf("| Perm p | %s | %s | %s | %s |\n",
              fmt3(pick("raw", "permutation_p_fc_corr")),
              fmt3(pick("bridge_shift", "permutation_p_fc_corr")),
              fmt3(pick("bridge_scale", "permutation_p_fc_corr")),
              fmt3(pick("celligner", "permutation_p_fc_corr"))))
  cat(sprintf("| Gene retention | — | %s | %s | %s |\n",
              fmt1pct(pick("bridge_shift", "biology_destruction_retention")),
              fmt1pct(pick("bridge_scale", "biology_destruction_retention")),
              fmt1pct(pick("celligner", "biology_destruction_retention"))))
  cat(sprintf("| FC shrinkage | — | %s | %s | %s |\n\n",
              fmt1pct(pick("bridge_shift", "biology_destruction_fc_shrinkage")),
              fmt1pct(pick("bridge_scale", "biology_destruction_fc_shrinkage")),
              fmt1pct(pick("celligner", "biology_destruction_fc_shrinkage"))))
}

write_the_disconnect_two_numbers <- function() {
  cs <- safe_fread(file.path(BENCH, "comparison_summary.csv"))
  if (is.null(cs)) return(invisible(NULL))
  bvl <- cs[task == "breast_vs_lung"]
  if (nrow(bvl) == 0) return(invisible(NULL))

  pick <- function(m, col) {
    r <- bvl[method == m]
    if (nrow(r) == 0 || !col %in% names(r)) return(NA_real_)
    as.numeric(r[[col]][1])
  }

  cat("SLIDE 15 — THE DISCONNECT\n\n")
  cat(sprintf("Raw domain R² PC1 (BvL):        %s\n", fmt3(pick("raw", "struct_domain_r2_pc1"))))
  cat(sprintf("Celligner domain R² PC1 (BvL):  %s\n\n", fmt3(pick("celligner", "struct_domain_r2_pc1"))))
  cat(sprintf("Raw FC corr intersection (BvL):       %s\n", fmt3(pick("raw", "fc_correlation_intersection"))))
  cat(sprintf("Celligner FC corr intersection (BvL): %s\n\n", fmt3(pick("celligner", "fc_correlation_intersection"))))
}

write_destruction_table <- function() {
  cs <- safe_fread(file.path(BENCH, "comparison_summary.csv"))
  if (is.null(cs)) return(invisible(NULL))

  cat("SLIDE 16 — DESTRUCTION\n\n")
  cat("| Method | Subtype retention | Subtype shrinkage | BvL retention | BvL shrinkage |\n")
  cat("|---|---:|---:|---:|---:|\n")
  for (m in c("raw", "bridge_shift", "bridge_scale", "celligner")) {
    sub <- cs[task == "breast_subtype" & method == m]
    bvl <- cs[task == "breast_vs_lung" & method == m]
    ret_s <- if (nrow(sub) && "biology_destruction_retention" %in% names(sub) && is.finite(sub$biology_destruction_retention[1])) fmt1pct(sub$biology_destruction_retention[1]) else "—"
    shr_s <- if (nrow(sub) && "biology_destruction_fc_shrinkage" %in% names(sub) && is.finite(sub$biology_destruction_fc_shrinkage[1])) fmt1pct(sub$biology_destruction_fc_shrinkage[1]) else "—"
    ret_b <- if (nrow(bvl) && "biology_destruction_retention" %in% names(bvl) && is.finite(bvl$biology_destruction_retention[1])) fmt1pct(bvl$biology_destruction_retention[1]) else "—"
    shr_b <- if (nrow(bvl) && "biology_destruction_fc_shrinkage" %in% names(bvl) && is.finite(bvl$biology_destruction_fc_shrinkage[1])) fmt1pct(bvl$biology_destruction_fc_shrinkage[1]) else "—"
    disp <- if (m == "raw") "Raw" else if (m == "bridge_shift") "Dom.Shift" else if (m == "bridge_scale") "Dom.Scale" else "Celligner"
    cat(sprintf("| %s | %s | %s | %s | %s |\n", disp, ret_s, shr_s, ret_b, shr_b))
  }
  cat("\n")
}

write_disconnect_scatter_table <- function() {
  disc <- safe_fread(file.path(BENCH, "disconnect_scores.csv"))
  if (is.null(disc)) return(invisible(NULL))
  cat("SLIDE 17 — DISCONNECT SCATTERPLOT\n\n")
  cat("| Method | Task | Geom | DA | Cost |\n")
  cat("|---|---|---:|---:|---:|\n")
  for (i in seq_len(nrow(disc))) {
    m <- as.character(disc$method[i])
    t <- as.character(disc$task[i])
    geom <- disc$geom_improvement[i] %||% disc$geom[i]
    da <- disc$da_improvement[i] %||% disc$da[i]
    cost <- disc$biology_cost[i] %||% disc$cost[i]
    disp <- if (m == "bridge_shift") "Dom.Shift" else if (m == "bridge_scale") "Dom.Scale" else if (m == "celligner") "Celligner" else m
    cat(sprintf("| %s | %s | %s | %s | %s |\n", disp, t, fmt2(geom), fmt2(da), fmt2(cost)))
  }
  cat("\n")
}

write_stratified_fc <- function() {
  cat("SLIDE 18 — STRATIFIED FC\n\n")
  for (task in c("breast_subtype", "breast_vs_lung")) {
    f <- file.path(REPO, "reports", "benchmark_master", "diagnostics", sprintf("fc_stratified_%s.csv", task))
    dt <- safe_fread(f)
    if (is.null(dt)) next
    cat(sprintf("%s:\n", ifelse(task == "breast_subtype", "Subtype", "BvL")))
    # Use raw method rows only for slide table
    if ("method" %in% names(dt)) dt <- dt[method == "raw"]
    sc <- "stratum"; nc <- "n_genes"; fc <- "fc_correlation"
    # print as markdown table
    cat("| Stratum | n_genes | FC_corr |\n|---|---:|---:|\n")
    for (i in seq_len(nrow(dt))) {
      cat(sprintf("| %s | %s | %s |\n", as.character(dt[[sc]][i]), as.character(dt[[nc]][i]), fmt3(dt[[fc]][i])))
    }
    cat("\n")
  }
}

write_fc_se_across_methods <- function() {
  cat("BACKUP B5 — FC/SE ACROSS METHODS (CPTAC subtype)\n\n")
  cat("| Method | n_sig | Median|logFC| | MedianSE |\n")
  cat("|---|---:|---:|---:|\n")
  for (m in methods) {
    p <- da_path(m, "breast_subtype", "cptac")
    da <- safe_fread(p)
    if (is.null(da)) next
    pcol <- intersect(c("adj.P.Val", "padj", "FDR"), names(da))[1]
    fccol <- intersect(c("logFC", "log2FC"), names(da))[1]
    tcol <- intersect(c("t", "t_stat"), names(da))[1]
    if (is.na(pcol) || is.na(fccol) || is.na(tcol)) next
    sig <- da[is.finite(get(pcol)) & get(pcol) < 0.05]
    n_sig <- nrow(sig)
    med_abs <- median(abs(sig[[fccol]]), na.rm = TRUE)
    med_se <- median(abs(sig[[fccol]] / sig[[tcol]]), na.rm = TRUE)
    disp <- if (m == "raw") "Raw" else if (m == "bridge_shift") "Dom.Shift" else if (m == "bridge_scale") "Dom.Scale" else "Celligner"
    cat(sprintf("| %s | %d | %s | %s |\n", disp, n_sig, fmt4(med_abs), fmt4(med_se)))
  }
  cat("\n")
}

read_celligner_all_header <- function() {
  p <- file.path(REPO, "reports", "benchmark_master", "celligner_all", "celligner_aligned_matrix.csv")
  if (!file.exists(p)) return(list(n_genes = NA_integer_, present = function(x) FALSE))
  hdr <- tryCatch(readLines(p, n = 1, warn = FALSE), error = function(e) character(0))
  if (length(hdr) == 0) return(list(n_genes = NA_integer_, present = function(x) FALSE))
  cols <- strsplit(hdr[[1]], ",", fixed = TRUE)[[1]]
  genes <- cols[-1] # first is sample id col (blank header)
  gset <- unique(genes)
  list(
    n_genes = length(gset),
    present = function(x) x %in% gset
  )
}

write_gene_set_comparison <- function() {
  cat("GENE SET COMPARISON\n\n")

  # Use DA headers for raw + celligner in v2 outputs (these are the actual benchmark gene sets)
  raw_c <- safe_fread(da_path("raw", "breast_subtype", "cptac"))
  cel_c <- safe_fread(da_path("celligner", "breast_subtype", "cptac"))
  gcol <- function(dt) intersect(c("gene", "Gene"), names(dt))[1]

  raw_genes <- if (!is.null(raw_c) && !is.na(gcol(raw_c))) unique(as.character(raw_c[[gcol(raw_c)]])) else character(0)
  cel_genes <- if (!is.null(cel_c) && !is.na(gcol(cel_c))) unique(as.character(cel_c[[gcol(cel_c)]])) else character(0)

  cat(sprintf("Celligner genes: %s\n", ifelse(length(cel_genes), as.character(length(cel_genes)), "NA")))
  cat(sprintf("Raw genes: %s\n", ifelse(length(raw_genes), as.character(length(raw_genes)), "NA")))
  cat(sprintf("Shared: %s\n", as.character(length(intersect(raw_genes, cel_genes)))))
  cat(sprintf("Raw only: %s\n", as.character(length(setdiff(raw_genes, cel_genes)))))
  cat(sprintf("Celligner only: %s\n\n", as.character(length(setdiff(cel_genes, raw_genes)))))

  subtype_markers <- c("ESR1", "PGR", "GATA3", "FOXA1", "EGFR", "KRT5", "KRT17", "FOXC1")
  bvl_markers <- c("NKX2-1", "NAPSA", "SFTPB", "GATA3", "FOXA1", "ESR1")

  cat("Subtype markers in Celligner:\n")
  for (g in subtype_markers) cat(sprintf("  %s: %s\n", g, ifelse(g %in% cel_genes, "yes", "no")))
  cat("\nBvL markers in Celligner:\n")
  # use bvl celligner CPTAC DA genes for membership
  cel_b <- safe_fread(da_path("celligner", "breast_vs_lung", "cptac"))
  cel_b_genes <- if (!is.null(cel_b) && !is.na(gcol(cel_b))) unique(as.character(cel_b[[gcol(cel_b)]])) else cel_genes
  for (g in bvl_markers) cat(sprintf("  %s: %s\n", g, ifelse(g %in% cel_b_genes, "yes", "no")))
  cat("\n")

  # Also check the all-data celligner header (legacy context)
  hdr <- read_celligner_all_header()
  cat(sprintf("Celligner-all (legacy) genes: %s\n", as.character(hdr$n_genes)))
  cat("Subtype markers in Celligner-all:\n")
  for (g in subtype_markers) cat(sprintf("  %s: %s\n", g, ifelse(hdr$present(g), "yes", "no")))
  cat("\n")
}

write_celligner_status <- function() {
  cat("CELLIGNER STATUS\n\n")
  c_path <- da_path("celligner", "breast_subtype", "cptac")
  r_path <- da_path("raw", "breast_subtype", "cptac")
  c_da <- safe_fread(c_path)
  r_da <- safe_fread(r_path)
  if (is.null(c_da) || is.null(r_da)) {
    cat("logFC correlation with raw (CPTAC subtype): NA\nStatus: NA\n\n")
    return(invisible(NULL))
  }
  gcol <- intersect(c("gene", "Gene"), names(c_da))[1]
  if (is.na(gcol)) {
    cat("logFC correlation with raw (CPTAC subtype): NA\nStatus: NA\n\n")
    return(invisible(NULL))
  }
  m <- merge(
    c_da[, .(gene = get(gcol), fc_cell = logFC)],
    r_da[, .(gene = get(gcol), fc_raw = logFC)],
    by = "gene"
  )
  cor_val <- suppressWarnings(cor(m$fc_cell, m$fc_raw, use = "complete.obs"))
  cat(sprintf("logFC correlation with raw (CPTAC subtype): %s\n", sprintf("%.4f", cor_val)))
  cat(sprintf("Status: %s\n\n", ifelse(is.finite(cor_val) && cor_val > 0.999, "SCAFFOLD", "REAL")))
}

write_ccle_annotation <- function() {
  cat("CCLE ANNOTATION (backup B10)\n\n")
  p1 <- file.path(REPO, "data", "processed", "ccle_breast_subtype_annotation_processed.csv")
  p2 <- file.path(REPO, "data", "ccle", "ccle_breast_subtype_annotations_v2.csv")
  ann <- safe_fread(p1)
  src <- p1
  if (is.null(ann)) { ann <- safe_fread(p2); src <- p2 }
  if (is.null(ann)) {
    cat("Annotation file missing.\n\n")
    return(invisible(NULL))
  }
  cat(sprintf("Source: %s\n", src))

  grp_col <- intersect(c("BvL_group", "group", "subtype", "Subtype", "label"), names(ann))[1]
  if (!is.na(grp_col)) {
    tab <- sort(table(ann[[grp_col]]), decreasing = TRUE)
    cat("Counts by subtype/group:\n")
    for (nm in names(tab)) cat(sprintf("  %s: %d\n", nm, tab[[nm]]))
  }

  plex_col <- intersect(c("plex", "Plex", "mixture", "Mixture", "TMTplex", "tmt_plex"), names(ann))
  if (length(plex_col) > 0) {
    pc <- plex_col[1]
    dup <- ann[, .N, by = pc][N > 1][order(-N)]
    cat(sprintf("\nLines sharing %s (N>1):\n", pc))
    if (nrow(dup) == 0) cat("  none\n")
    else print(dup)
  }
  cat("\n")
}

write_talk_cheat_sheet <- function() {
  cs <- safe_fread(file.path(BENCH, "comparison_summary.csv"))
  if (is.null(cs)) return(invisible(NULL))
  sub_raw <- cs[method == "raw" & task == "breast_subtype"]
  sub_cell <- cs[method == "celligner" & task == "breast_subtype"]
  bvl_raw <- cs[method == "raw" & task == "breast_vs_lung"]
  bvl_cell <- cs[method == "celligner" & task == "breast_vs_lung"]

  cat("TALK CHEAT SHEET\n\n")
  cat(sprintf("1. Subtype FC corr (raw):           %s\n", fmt3(sub_raw$fc_correlation_intersection[1])))
  cat(sprintf("2. Subtype perm z (raw):            %s\n", fmt2(sub_raw$permutation_z_fc_corr[1])))
  cat(sprintf("3. Subtype ceiling:                 %s\n", fmt2(sub_raw$concordance_ceiling_fc_corr[1])))
  cat(sprintf("4. Subtype calibrated:              %s\n", fmt1pct(sub_raw$calibrated_fc_corr_intersection[1])))
  cat(sprintf("5. Celligner subtype FC corr:       %s\n", fmt3(sub_cell$fc_correlation_intersection[1])))
  cat(sprintf("6. Celligner subtype retention:     %s\n", fmt1pct(sub_cell$biology_destruction_retention[1])))
  cat(sprintf("7. BvL domain R² raw:               %s\n", fmt3(bvl_raw$struct_domain_r2_pc1[1])))
  cat(sprintf("8. BvL domain R² Celligner:         %s\n", fmt3(bvl_cell$struct_domain_r2_pc1[1])))
  cat(sprintf("9. BvL FC corr raw:                 %s\n", fmt3(bvl_raw$fc_correlation_intersection[1])))
  cat(sprintf("10. BvL FC corr Celligner:          %s\n", fmt3(bvl_cell$fc_correlation_intersection[1])))
  cat("\n")
}

con <- file(OUT_MD, open = "wt")
on.exit(close(con), add = TRUE)
sink(con)
on.exit(sink(NULL), add = TRUE)

write_slide_counts()
cat("\n")
write_gene_coverage()
write_raw_da_counts()
write_subtype_cross_domain()

# Slide 11 subtype markers
sub_markers <- c("ESR1", "PGR", "GATA3", "FOXA1", "EGFR", "KRT5", "KRT17", "FOXC1")
sub_exp <- c("Lum ↑ (+)", "Lum ↑ (+)", "Lum ↑ (+)", "Lum ↑ (+)", "Bas ↑ (−)", "Bas ↑ (−)", "Bas ↑ (−)", "Bas ↑ (−)")
write_marker_table("breast_subtype", sub_markers, sub_exp, c(1,1,1,1,-1,-1,-1,-1))

write_bvl_geometry()

# Slide 13 BvL markers
cat("NOTE: DA files use logFC = Lung − Breast. Expected signs below reflect that.\n\n")
bvl_markers <- c("NKX2-1", "NAPSA", "SFTPB", "GATA3", "FOXA1", "ESR1")
bvl_exp <- c("Lung ↑ (+)", "Lung ↑ (+)", "Lung ↑ (+)", "Breast ↑ (−)", "Breast ↑ (−)", "Breast ↑ (−)")
write_marker_table("breast_vs_lung", bvl_markers, bvl_exp, c(1,1,1,-1,-1,-1))

write_bvl_da()
write_the_disconnect_two_numbers()
write_destruction_table()
write_disconnect_scatter_table()
write_stratified_fc()
write_fc_se_across_methods()
write_gene_set_comparison()
write_celligner_status()
write_ccle_annotation()
write_talk_cheat_sheet()

sink(NULL)
cat(sprintf("Wrote %s\n", OUT_MD))

