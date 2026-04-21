#!/usr/bin/env Rscript
# Extract every number needed for the 20-slide presentation.
# Output: presentation_materials/SLIDE_DATA.md
#
# Reads from:
# - reports/benchmark_master/benchmark_results/ (v2 tables + DA)
# - reports/benchmark_master/ (rerun_full_benchmark outputs)

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

BENCH <- file.path(REPO, "reports", "benchmark_master", "benchmark_results")
DIAG  <- file.path(REPO, "reports", "benchmark_master", "diagnostics")
UNION <- file.path(REPO, "data", "processed", "union")

methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
tasks <- c("breast_subtype", "breast_vs_lung")

sink_file <- file.path(PRES_OUT, "SLIDE_DATA.md")

safe_read <- function(path) {
  if (!file.exists(path)) {
    cat(sprintf("MISSING: %s\n", path))
    return(NULL)
  }
  tryCatch(fread(path), error = function(e) {
    cat(sprintf("ERROR reading %s: %s\n", path, e$message))
    NULL
  })
}

safe_read_lines <- function(path) {
  if (!file.exists(path)) {
    cat(sprintf("MISSING: %s\n", path))
    return(NULL)
  }
  tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
}

md_h2 <- function(txt) cat("\n## ", txt, "\n\n", sep = "")
md_h3 <- function(txt) cat("\n### ", txt, "\n\n", sep = "")
md_rule <- function() cat("\n---\n\n")

comparison_path <- file.path(BENCH, "comparison_summary.csv")
comparison_tiered_path <- file.path(BENCH, "comparison_summary_tiered.csv")
disconnect_path <- file.path(BENCH, "disconnect_scores.csv")

cs <- safe_read(comparison_path)
if (is.null(cs)) cs <- safe_read(comparison_tiered_path)
disc <- safe_read(disconnect_path)

da_path_v2 <- function(method, task, domain) {
  p1 <- file.path(BENCH, method, task, "representation_da", domain, "da_limma_result.csv")
  if (file.exists(p1)) return(p1)
  p2 <- file.path(BENCH, method, task, "representation_da", sprintf("da_%s.csv", domain))
  if (file.exists(p2)) return(p2)
  NA_character_
}

write_md <- function() {
  cat("# SLIDE DATA — extracted ", format(Sys.time()), "\n\n", sep = "")
  cat("Copy-paste numbers from here into your slides.\n\n")
  cat("- Repo: `", REPO, "`\n", sep = "")
  cat("- v2 benchmark: `reports/benchmark_master/benchmark_results/`\n")
  cat("- Union metadata: `data/processed/union/`\n\n")

  md_rule()
  md_h2("SLIDE 3: DATA — sample counts")
  for (task in tasks) {
    meta_path <- file.path(UNION, sprintf("sample_meta_%s.csv", task))
    meta <- safe_read(meta_path)
    if (is.null(meta)) next
    md_h3(task)
    cat(sprintf("- Total samples: **%d**\n", nrow(meta)))
    if ("domain" %in% names(meta)) {
      dom_tab <- sort(table(meta$domain), decreasing = TRUE)
      cat("- Domain counts:\n")
      for (d in names(dom_tab)) cat(sprintf("  - %s: **%d**\n", d, dom_tab[[d]]))
    }
    if (all(c("domain", "condition") %in% names(meta))) {
      cat("\n- Condition × Domain:\n\n")
      print(as.data.table(as.matrix(table(meta$condition, meta$domain))))
    }
    cat("\n")
  }

  md_rule()
  md_h2("SLIDE 3: Gene coverage")
  for (task in tasks) {
    audit_path <- file.path(DIAG, sprintf("gene_coverage_audit_%s.csv", task))
    audit <- safe_read(audit_path)
    if (!is.null(audit)) {
      md_h3(task)
      if ("category" %in% names(audit)) {
        cat_tab <- sort(table(audit$category), decreasing = TRUE)
        for (ct in names(cat_tab)) cat(sprintf("- %s: **%d**\n", ct, cat_tab[[ct]]))
        cat(sprintf("- Total genes (rows): **%d**\n\n", nrow(audit)))
      } else {
        cat(sprintf("- Rows: **%d**\n\n", nrow(audit)))
      }
    }
    int_path <- file.path(REPO, "data", "processed", sprintf("intersection_genes_%s.txt", task))
    if (file.exists(int_path)) {
      n <- length(safe_read_lines(int_path) %||% character(0))
      cat(sprintf("- Intersection genes list (%s): **%d**\n\n", task, n))
    }
  }

  md_rule()
  md_h2("SLIDE 6: DA gene counts per domain (raw, v2 DA tables)")
  for (task in tasks) {
    for (domain in c("cptac", "ccle")) {
      p <- da_path_v2("raw", task, domain)
      if (is.na(p)) next
      da <- safe_read(p)
      if (is.null(da)) next
      pcol <- intersect(c("adj.P.Val", "padj", "FDR", "qval"), names(da))
      fccol <- intersect(c("logFC", "log2FC", "log2FoldChange"), names(da))
      if (length(pcol) == 0L || length(fccol) == 0L) next
      pcol <- pcol[1L]; fccol <- fccol[1L]
      sig <- da[is.finite(get(pcol)) & get(pcol) < 0.05]

      md_h3(sprintf("%s / %s", task, toupper(domain)))
      cat(sprintf("- Genes tested: **%d**\n", nrow(da)))
      cat(sprintf("- Significant (FDR<0.05): **%d**\n", nrow(sig)))
      cat(sprintf("- Up among sig (logFC>0): **%d**\n", sum(sig[[fccol]] > 0, na.rm = TRUE)))
      cat(sprintf("- Down among sig (logFC<0): **%d**\n", sum(sig[[fccol]] < 0, na.rm = TRUE)))
      cat(sprintf("- Median |logFC| all: **%.4f**\n", median(abs(da[[fccol]]), na.rm = TRUE)))
      cat(sprintf("- Median |logFC| sig: **%.4f**\n", median(abs(sig[[fccol]]), na.rm = TRUE)))
      cat("\n")
    }
  }

  md_rule()
  md_h2("SLIDE 10: Breast subtype — cross-domain summary (v2 comparison_summary)")
  if (!is.null(cs) && all(c("method", "task") %in% names(cs))) {
    sub <- cs[task == "breast_subtype"]
    if (nrow(sub) > 0L) {
      key_cols <- c(
        "method",
        "fc_correlation_intersection", "fc_same_dir_intersection", "n_genes_intersection",
        "permutation_z_fc_corr", "permutation_p_fc_corr",
        "concordance_ceiling_fc_corr", "calibrated_fc_corr_intersection",
        "marker_sanity_cptac", "marker_sanity_ccle",
        "biology_destruction_retention", "biology_destruction_fc_shrinkage",
        "disconnect_score"
      )
      present <- intersect(key_cols, names(sub))
      for (i in seq_len(nrow(sub))) {
        md_h3(as.character(sub$method[i]))
        for (col in present[present != "method"]) {
          cat(sprintf("- %s: **%s**\n", col, as.character(sub[[col]][i])))
        }
        cat("\n")
      }
    }
  }

  md_rule()
  md_h2("SLIDE 14: Disconnect scores (v2)")
  if (!is.null(disc)) {
    print(disc)
    cat("\n")
  }
}

con <- file(sink_file, open = "wt")
on.exit(close(con), add = TRUE)
sink(con)
on.exit(sink(NULL), add = TRUE)

write_md()

# Copy 4 main figures (if they exist)
main_dir <- file.path(PRES_OUT, "figures", "main")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
copy_if_exists <- function(src) {
  if (!file.exists(src)) return(FALSE)
  file.copy(src, file.path(main_dir, basename(src)), overwrite = TRUE)
}
copy_if_exists(file.path(PRES_OUT, "figures", "meeting", "01_contrast_validation_table.png"))
copy_if_exists(file.path(PRES_OUT, "figures", "meeting", "02_concordance_ceiling_context.png"))
copy_if_exists(file.path(PRES_OUT, "figures", "meeting", "03_comparison_table.png"))
copy_if_exists(file.path(PRES_OUT, "figures", "marker_agreement", "marker_agreement_raw.png"))

sink(NULL)
cat(sprintf("Wrote %s\n", sink_file))
cat(sprintf("Copied main figures to %s\n", main_dir))

