#!/usr/bin/env Rscript
# One figure per method × task: marker logFC CPTAC vs CCLE, colored by expected
# direction under the task contrast.
#
# Reads:  reports/benchmark_master/representation_level_da/{method}/{task}/{cptac,ccle}/marker_summary.csv
# Writes: presentation_materials/figures/marker_agreement/marker_agreement_{task}_{method}.pdf (+ .png)

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(a, b) if (is.null(a)) b else a

.args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
if ("--repo-root" %in% .args) {
  i <- which(.args == "--repo-root")[[1L]]
  if (i < length(.args)) repo_root <- .args[[i + 1L]]
}
repo_root <- normalizePath(repo_root, mustWork = TRUE)

base <- file.path(repo_root, "reports", "benchmark_master", "representation_level_da")
out_root <- file.path(repo_root, "presentation_materials", "figures", "marker_agreement")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

task_spec <- list(
  breast_subtype = list(
    label = "Luminal − Basal",
    genes = c("FOXA1", "GATA3", "ESR1", "PGR", "ERBB2", "CDH1", "EGFR", "KRT5", "KRT14", "KRT17"),
    expected = c(FOXA1 = 1L, GATA3 = 1L, ESR1 = 1L, PGR = 1L, ERBB2 = 1L, CDH1 = 1L,
                 EGFR = -1L, KRT5 = -1L, KRT14 = -1L, KRT17 = -1L)
  ),
  breast_vs_lung = list(
    label = "Lung − Breast",
    # Requested minimal marker panel for slides
    genes = c("NKX2-1", "SFTPB", "NAPSA", "GATA3", "FOXA1", "ESR1"),
    # +1 means higher in Lung; -1 higher in Breast
    expected = setNames(
      c(1L, 1L, 1L, -1L, -1L, -1L),
      c("NKX2-1", "SFTPB", "NAPSA", "GATA3", "FOXA1", "ESR1")
    )
  )
)

col_ok <- "#2E7D32"
col_bad <- "#C62828"
col_miss <- "#BDBDBD"

agrees <- function(logfc, exp_sign) {
  s <- sign(logfc)
  is.finite(logfc) && s != 0 && as.integer(s) == exp_sign
}

plot_one_method_task <- function(method, task) {
  spec <- task_spec[[task]]
  genes <- spec$genes
  expected <- spec$expected
  out_dir <- file.path(out_root, task)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Display names for slides (keep on-disk method folder names unchanged)
  method_display <- switch(
    method,
    bridge_shift = "domain_shift",
    bridge_scale = "domain_scale",
    method
  )

  f_c <- file.path(base, method, task, "cptac", "marker_summary.csv")
  f_e <- file.path(base, method, task, "ccle", "marker_summary.csv")
  if (!file.exists(f_c) || !file.exists(f_e)) {
    message("Skip ", task, " / ", method, " — missing marker_summary.csv")
    return(invisible(NULL))
  }
  dc <- fread(f_c, select = c("gene", "logFC"))
  de <- fread(f_e, select = c("gene", "logFC"))
  mi <- match(genes, dc$gene)
  fc <- ifelse(is.na(mi), NA_real_, as.numeric(dc$logFC[mi]))
  me <- match(genes, de$gene)
  fe <- ifelse(is.na(me), NA_real_, as.numeric(de$logFC[me]))

  col_c <- vapply(seq_along(genes), function(i) {
    v <- fc[i]
    if (is.na(v)) return(col_miss)
    if (agrees(v, expected[[genes[i]]])) col_ok else col_bad
  }, character(1L))
  col_e <- vapply(seq_along(genes), function(i) {
    v <- fe[i]
    if (is.na(v)) return(col_miss)
    if (agrees(v, expected[[genes[i]]])) col_ok else col_bad
  }, character(1L))

  n_ok_c <- sum(mapply(function(g, v) !is.na(v) && agrees(v, expected[[g]]), genes, fc, SIMPLIFY = TRUE))
  n_pr_c <- sum(!is.na(fc))
  n_ok_e <- sum(mapply(function(g, v) !is.na(v) && agrees(v, expected[[g]]), genes, fe, SIMPLIFY = TRUE))
  n_pr_e <- sum(!is.na(fe))

  pdf_path <- file.path(out_dir, paste0("marker_agreement_", method_display, ".pdf"))
  pdf(pdf_path, width = 11, height = 5.5)

  par(mar = c(7, 4, 4, 2), mgp = c(2.2, 0.6, 0))
  n <- length(genes)
  x <- seq_len(n)
  w <- 0.32
  off <- 0.2
  ylim <- range(c(0, fc, fe), na.rm = TRUE)
  if (!is.finite(diff(ylim))) ylim <- c(-1, 1)
  ylim <- ylim + c(-1, 1) * 0.08 * max(diff(ylim), 0.5)

  plot(NA, xlim = c(0.5, n + 0.5), ylim = ylim, xaxt = "n", xlab = "",
       ylab = paste0("log2 FC (", spec$label, ")"),
       main = paste0("Marker direction vs biology — ", method_display, " (", task, ")"))
  abline(h = 0, col = "black", lwd = 1)
  for (i in seq_len(n)) {
    xc <- i - off
    xe <- i + off
    if (!is.na(fc[i])) {
      rect(xc - w / 2, 0, xc + w / 2, fc[i], col = col_c[i], border = "white", lwd = 0.5)
    } else {
      text(xc, 0, "—", col = "#616161", cex = 1.1)
    }
    if (!is.na(fe[i])) {
      rect(xe - w / 2, 0, xe + w / 2, fe[i], col = col_e[i], border = "white", lwd = 0.5)
    } else {
      text(xe, 0, "—", col = "#616161", cex = 1.1)
    }
  }
  axis(1, at = x, labels = FALSE)
  text(x, par("usr")[[3L]] - 0.08 * diff(par("usr")[3:4]), genes, srt = 35, adj = c(1, 1), xpd = TRUE, cex = 0.85)

  legend("topright",
         legend = c("Sign matches biology", "Sign opposite", "Absent from DA table"),
         fill = c(col_ok, col_bad, col_miss), border = NA, bty = "n", cex = 0.85)
  mtext(paste0("CPTAC: ", n_ok_c, "/", n_pr_c, " correct (present)  |  CCLE: ", n_ok_e, "/", n_pr_e, " correct (present)"),
        side = 1, line = 5.2, cex = 0.78, col = "gray20")
  mtext("Left bar = CPTAC, right bar = CCLE per gene", side = 3, line = 0.3, cex = 0.78, col = "gray30")

  dev.off()
  message("Wrote ", pdf_path)

  png_path <- file.path(out_dir, paste0("marker_agreement_", method_display, ".png"))
  tryCatch({
    png(png_path, width = 1100, height = 550, res = 120)
    par(mar = c(7, 4, 4, 2), mgp = c(2.2, 0.6, 0))
    plot(NA, xlim = c(0.5, n + 0.5), ylim = ylim, xaxt = "n", xlab = "",
         ylab = paste0("log2 FC (", spec$label, ")"),
         main = paste0("Marker direction vs biology — ", method_display, " (", task, ")"))
    abline(h = 0, col = "black", lwd = 1)
    for (i in seq_len(n)) {
      xc <- i - off
      xe <- i + off
      if (!is.na(fc[i])) {
        rect(xc - w / 2, 0, xc + w / 2, fc[i], col = col_c[i], border = "white", lwd = 0.5)
      } else {
        text(xc, 0, "—", col = "#616161", cex = 1.1)
      }
      if (!is.na(fe[i])) {
        rect(xe - w / 2, 0, xe + w / 2, fe[i], col = col_e[i], border = "white", lwd = 0.5)
      } else {
        text(xe, 0, "—", col = "#616161", cex = 1.1)
      }
    }
    axis(1, at = x, labels = FALSE)
    text(x, par("usr")[[3L]] - 0.08 * diff(par("usr")[3:4]), genes, srt = 35, adj = c(1, 1), xpd = TRUE, cex = 0.85)
    legend("topright",
           legend = c("Sign matches biology", "Sign opposite", "Absent from DA table"),
           fill = c(col_ok, col_bad, col_miss), border = NA, bty = "n", cex = 0.85)
    mtext(paste0("CPTAC: ", n_ok_c, "/", n_pr_c, " correct (present)  |  CCLE: ", n_ok_e, "/", n_pr_e, " correct (present)"),
          side = 1, line = 5.2, cex = 0.78, col = "gray20")
    mtext("Left bar = CPTAC, right bar = CCLE per gene", side = 3, line = 0.3, cex = 0.78, col = "gray30")
    dev.off()
    message("Wrote ", png_path)
  }, error = function(e) message("PNG skipped: ", conditionMessage(e)))
}

for (task in names(task_spec)) {
  for (m in c("raw", "celligner", "bridge_shift", "bridge_scale")) {
    plot_one_method_task(m, task)
  }
}

message("Done. Output: ", out_root)
