#!/usr/bin/env Rscript
# =============================================================================
# Method-Agnostic Sample Profile Plots for Benchmark Marker Proteins
# =============================================================================
#
# INPUT CONTRACT
# ==============
# Each representation is a named list element containing:
#   $matrix   — numeric matrix, genes (rows) × samples (cols), rownames = gene IDs
#   $name     — character, display name for this method
#
# The plotting system additionally receives:
#   task_meta — data.table with columns: sample_id, condition, domain
#   markers   — character vector of gene identifiers to plot
#   task_name — character label for the task
#   outdir    — path for output files
#
# OUTPUT: One combined summary figure per task (markers × methods grid),
#         plus a compact boxplot summary. Replaces per-marker files.
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ─── Build deterministic sample order ────────────────────────────────────
build_sample_order <- function(task_meta, block_order = NULL) {
  stopifnot(all(c("sample_id", "condition", "domain") %in% names(task_meta)))
  m <- copy(task_meta)
  if (!is.null(block_order)) {
    m[, block := paste(domain, condition, sep = " ")]
    m[, block := factor(block, levels = block_order)]
  } else {
    m[, block := factor(paste(domain, condition, sep = " "))]
  }
  m <- m[order(block, sample_id)]
  m[, x_pos := seq_len(.N)]
  m
}

# ─── Block geometry helpers ──────────────────────────────────────────────
compute_blocks <- function(blocks, block_levels) {
  block_ends <- integer(0)
  block_mids <- numeric(0)
  active_levels <- character(0)
  for (bl in block_levels) {
    idx <- which(blocks == bl)
    if (length(idx) == 0) next
    block_ends <- c(block_ends, max(idx))
    block_mids <- c(block_mids, mean(range(idx)))
    active_levels <- c(active_levels, bl)
  }
  names(block_mids) <- active_levels
  list(ends = block_ends, mids = block_mids, active = active_levels)
}

block_palette <- c(
  "#2196F3", "#E91E63", "#64B5F6", "#F48FB1",
  "#4CAF50", "#FF9800", "#81C784", "#FFB74D",
  "#9C27B0", "#00BCD4", "#795548", "#607D8B"
)

# ─── Extract values for one marker in one method ────────────────────────
get_marker_vals <- function(marker, mat, sids) {
  if (!marker %in% rownames(mat)) return(NULL)
  common <- intersect(sids, colnames(mat))
  v <- rep(NA_real_, length(sids))
  names(v) <- sids
  v[common] <- mat[marker, common]
  v
}

# ─── Draw one cell of the grid (one marker × one method) ────────────────
draw_cell <- function(v, ylim, blocks, block_levels, block_info, block_cols,
                       sample_cols, n_samples, marker_label, method_label,
                       show_x_labels, show_y_axis, is_absent) {
  if (is_absent) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
         xlim = c(0.5, n_samples + 0.5), ylim = ylim)
    text(n_samples / 2, mean(ylim), "n/a", cex = 0.8, col = "gray60")
    box(col = "gray85")
    return(invisible(NULL))
  }

  plot(seq_len(n_samples), v, type = "n",
       xlim = c(0.5, n_samples + 0.5), ylim = ylim,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n")

  if (show_y_axis) axis(2, cex.axis = 0.55, las = 1, tck = -0.03, mgp = c(2, 0.3, 0))

  # Block backgrounds
  for (bl in block_levels) {
    idx <- which(blocks == bl)
    if (length(idx) == 0) next
    rect(min(idx) - 0.5, ylim[1], max(idx) + 0.5, ylim[2],
         col = adjustcolor(block_cols[bl], alpha.f = 0.08), border = NA)
  }

  # Block separators
  for (be in block_info$ends[-length(block_info$ends)]) {
    abline(v = be + 0.5, col = "gray50", lwd = 0.8, lty = 3)
  }

  # Group medians
  for (bl in block_levels) {
    idx <- which(blocks == bl)
    if (length(idx) == 0) next
    gmed <- median(v[idx], na.rm = TRUE)
    if (is.finite(gmed)) {
      segments(min(idx) - 0.3, gmed, max(idx) + 0.3, gmed,
               col = adjustcolor(block_cols[bl], alpha.f = 0.7),
               lwd = 1.8)
    }
  }

  points(seq_len(n_samples), v, pch = 16, cex = 0.3,
         col = adjustcolor(sample_cols, alpha.f = 0.6))
}

# ═══════════════════════════════════════════════════════════════════════════
# COMBINED GRID: markers (rows) × methods (columns) — sample profiles
# ═══════════════════════════════════════════════════════════════════════════
make_sample_profile_grid <- function(representations, task_meta, markers,
                                      task_name, outdir,
                                      block_order = NULL,
                                      fmt = c("png", "pdf")) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ordered_meta <- build_sample_order(task_meta, block_order)

  sids <- ordered_meta$sample_id
  blocks <- ordered_meta$block
  block_levels <- levels(blocks)
  n_samples <- length(sids)
  method_names <- names(representations)
  n_methods <- length(method_names)

  block_info <- compute_blocks(blocks, block_levels)
  block_cols <- setNames(block_palette[seq_along(block_levels)], block_levels)
  sample_cols <- block_cols[as.character(blocks)]

  # Pre-extract all values and compute per-marker y-ranges
  all_vals <- list()
  avail <- matrix(FALSE, nrow = length(markers), ncol = n_methods,
                   dimnames = list(markers, method_names))
  ylims <- list()

  for (g in markers) {
    all_vals[[g]] <- list()
    finite_vals <- numeric(0)
    for (nm in method_names) {
      v <- get_marker_vals(g, representations[[nm]]$matrix, sids)
      all_vals[[g]][[nm]] <- v
      if (!is.null(v)) {
        avail[g, nm] <- TRUE
        finite_vals <- c(finite_vals, v[is.finite(v)])
      }
    }
    if (length(finite_vals) >= 2) {
      rng <- range(finite_vals)
      pad <- diff(rng) * 0.08
      ylims[[g]] <- rng + c(-pad, pad)
    } else {
      ylims[[g]] <- c(0, 1)
    }
  }

  plottable <- markers[rowSums(avail) > 0]
  if (length(plottable) == 0) {
    cat("  No plottable markers.\n")
    return(invisible(NULL))
  }
  n_markers <- length(plottable)

  # Dimensions
  cell_w <- max(2.5, n_samples * 0.022 + 0.8)
  cell_h <- 1.3
  left_margin <- 0.7
  top_margin <- 0.8
  bottom_margin <- 0.6
  fig_w <- left_margin + cell_w * n_methods + 0.3
  fig_h <- top_margin + cell_h * n_markers + bottom_margin

  do_grid <- function(device_fn, path, ...) {
    device_fn(path, width = fig_w, height = fig_h, ...)
    layout_mat <- matrix(seq_len(n_markers * n_methods), nrow = n_markers,
                          ncol = n_methods, byrow = TRUE)
    layout(layout_mat,
           widths = rep(cell_w, n_methods),
           heights = rep(cell_h, n_markers))

    for (ri in seq_along(plottable)) {
      g <- plottable[ri]
      for (ci in seq_along(method_names)) {
        nm <- method_names[ci]
        show_y <- (ci == 1)
        show_x <- (ri == n_markers)

        lm <- c(0.3, ifelse(show_y, 2.8, 0.5), ifelse(ri == 1, 1.8, 0.3), 0.3)
        par(mar = lm, mgp = c(1.5, 0.3, 0))

        v <- all_vals[[g]][[nm]]
        draw_cell(v, ylims[[g]], blocks, block_levels, block_info, block_cols,
                   sample_cols, n_samples, g, representations[[nm]]$name,
                   show_x, show_y, is.null(v))

        # Row label (marker) on left edge of first column
        if (ci == 1) {
          mtext(g, side = 2, line = 1.8, cex = 0.55, font = 2, las = 1)
        }

        # Column label (method) on top of first row
        if (ri == 1) {
          mtext(representations[[nm]]$name, side = 3, line = 0.4, cex = 0.55, font = 2)
        }

        # Block labels at bottom of last row
        if (show_x) {
          for (bl in block_info$active) {
            mtext(bl, side = 1, line = 0.2, at = block_info$mids[bl],
                  cex = 0.35, col = block_cols[bl], font = 2)
          }
        }
      }
    }

    # Overall title
    mtext(paste0("Sample Profiles — ", gsub("_", " ", task_name)),
          side = 3, outer = TRUE, line = -1.2, cex = 0.9, font = 2)

    dev.off()
  }

  paths <- character(0)
  fname <- paste0("profile_grid_", gsub(" ", "_", task_name))
  if ("png" %in% fmt) {
    p <- file.path(outdir, paste0(fname, ".png"))
    do_grid(png, p, res = 160, units = "in")
    paths <- c(paths, p)
    cat("  Grid (png):", basename(p), "\n")
  }
  if ("pdf" %in% fmt) {
    p <- file.path(outdir, paste0(fname, ".pdf"))
    do_grid(pdf, p)
    paths <- c(paths, p)
    cat("  Grid (pdf):", basename(p), "\n")
  }

  list(paths = paths, plottable = plottable,
       skipped = setdiff(markers, plottable), avail = avail)
}

# ═══════════════════════════════════════════════════════════════════════════
# COMPACT BOXPLOT SUMMARY: markers (rows) × methods (columns)
# ═══════════════════════════════════════════════════════════════════════════
make_boxplot_summary <- function(representations, task_meta, markers,
                                  task_name, outdir,
                                  block_order = NULL,
                                  fmt = c("png", "pdf")) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ordered_meta <- build_sample_order(task_meta, block_order)

  sids <- ordered_meta$sample_id
  blocks <- ordered_meta$block
  block_levels <- levels(blocks)
  method_names <- names(representations)
  n_methods <- length(method_names)

  block_cols <- setNames(block_palette[seq_along(block_levels)], block_levels)

  # Check availability
  plottable <- character(0)
  for (g in markers) {
    for (nm in method_names) {
      if (g %in% rownames(representations[[nm]]$matrix)) { plottable <- c(plottable, g); break }
    }
  }
  plottable <- unique(plottable)
  n_markers <- length(plottable)
  if (n_markers == 0) return(invisible(NULL))

  n_blocks <- length(block_levels)
  cell_w <- max(1.8, n_blocks * 0.45)
  cell_h <- 1.4
  fig_w <- 0.8 + cell_w * n_methods + 0.3
  fig_h <- 0.9 + cell_h * n_markers + 0.5

  do_bp <- function(device_fn, path, ...) {
    device_fn(path, width = fig_w, height = fig_h, ...)
    layout_mat <- matrix(seq_len(n_markers * n_methods), nrow = n_markers,
                          ncol = n_methods, byrow = TRUE)
    layout(layout_mat,
           widths = rep(cell_w, n_methods),
           heights = rep(cell_h, n_markers))

    for (ri in seq_along(plottable)) {
      g <- plottable[ri]

      # Compute shared y-range for this marker across all methods
      all_finite <- numeric(0)
      for (nm in method_names) {
        v <- get_marker_vals(g, representations[[nm]]$matrix, sids)
        if (!is.null(v)) all_finite <- c(all_finite, v[is.finite(v)])
      }
      if (length(all_finite) < 2) { ylim <- c(0, 1) } else {
        rng <- range(all_finite); pad <- diff(rng) * 0.1
        ylim <- rng + c(-pad, pad)
      }

      for (ci in seq_along(method_names)) {
        nm <- method_names[ci]
        show_y <- (ci == 1)

        lm <- c(0.3, ifelse(show_y, 2.8, 0.5), ifelse(ri == 1, 1.8, 0.3), 0.3)
        par(mar = lm, mgp = c(1.5, 0.3, 0))

        v <- get_marker_vals(g, representations[[nm]]$matrix, sids)
        if (is.null(v)) {
          plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
               xlim = c(0.5, n_blocks + 0.5), ylim = ylim)
          text(n_blocks / 2, mean(ylim), "n/a", cex = 0.7, col = "gray60")
          box(col = "gray85")
        } else {
          gvals <- lapply(block_levels, function(bl) {
            idx <- which(blocks == bl)
            v[idx]
          })
          names(gvals) <- block_levels

          boxplot(gvals, ylim = ylim, outline = FALSE, xaxt = "n",
                  yaxt = ifelse(show_y, "s", "n"), cex.axis = 0.5,
                  border = block_cols[block_levels],
                  col = adjustcolor(block_cols[block_levels], alpha.f = 0.15),
                  lwd = 0.8, medlwd = 1.5, whisklty = 1)
          stripchart(gvals, vertical = TRUE, method = "jitter",
                     jitter = 0.15, pch = 16, cex = 0.25,
                     col = adjustcolor(block_cols[block_levels], alpha.f = 0.4),
                     add = TRUE)
        }

        if (ci == 1) mtext(g, side = 2, line = 1.8, cex = 0.55, font = 2, las = 1)
        if (ri == 1) mtext(representations[[nm]]$name, side = 3, line = 0.4, cex = 0.55, font = 2)
        if (ri == n_markers) {
          short_labels <- gsub("CPTAC ", "C:", gsub("CCLE ", "E:", block_levels))
          axis(1, at = seq_along(block_levels), labels = short_labels,
               cex.axis = 0.35, tck = -0.02, mgp = c(1, 0.1, 0), las = 2)
        }
      }
    }

    mtext(paste0("Marker Boxplots — ", gsub("_", " ", task_name)),
          side = 3, outer = TRUE, line = -1.2, cex = 0.9, font = 2)
    dev.off()
  }

  paths <- character(0)
  fname <- paste0("boxplot_summary_", gsub(" ", "_", task_name))
  if ("png" %in% fmt) {
    p <- file.path(outdir, paste0(fname, ".png"))
    do_bp(png, p, res = 160, units = "in")
    paths <- c(paths, p)
    cat("  Boxplot (png):", basename(p), "\n")
  }
  if ("pdf" %in% fmt) {
    p <- file.path(outdir, paste0(fname, ".pdf"))
    do_bp(pdf, p)
    paths <- c(paths, p)
    cat("  Boxplot (pdf):", basename(p), "\n")
  }

  list(paths = paths, plottable = plottable)
}

# ═══════════════════════════════════════════════════════════════════════════
# TOP-LEVEL: generate all outputs for a task
# ═══════════════════════════════════════════════════════════════════════════
make_sample_profile_plots <- function(representations, task_meta, markers,
                                        task_name, outdir,
                                        block_order = NULL,
                                        fmt = c("png", "pdf")) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  cat("  Generating sample profile grid...\n")
  grid_res <- make_sample_profile_grid(representations, task_meta, markers,
                                         task_name, outdir, block_order, fmt)

  cat("  Generating boxplot summary...\n")
  bp_res <- make_boxplot_summary(representations, task_meta, markers,
                                   task_name, outdir, block_order, fmt)

  # Warnings log
  mnames <- names(representations)
  avail <- grid_res$avail
  lines <- c(
    paste0("# Profile Plot Warnings — ", task_name),
    paste0("\nGenerated: ", Sys.time()),
    paste0("Methods: ", paste(mnames, collapse = ", ")),
    paste0("Markers requested: ", length(markers)),
    paste0("Markers plotted: ", length(grid_res$plottable)),
    paste0("Markers skipped: ", length(grid_res$skipped)),
    ""
  )
  if (length(grid_res$skipped) > 0) {
    lines <- c(lines, "## Skipped Markers", "")
    for (g in grid_res$skipped) {
      absent <- mnames[!avail[g, ]]
      lines <- c(lines, paste0("- **", g, "**: absent in all methods"))
    }
    lines <- c(lines, "")
  }

  hdr <- paste0("| Marker | ", paste(mnames, collapse = " | "), " |")
  sepr <- paste0("|--------|", paste(rep("---", length(mnames)), collapse = "|"), "|")
  lines <- c(lines, "## Marker x Method Availability", "", hdr, sepr)
  for (g in markers) {
    cells <- sapply(mnames, function(nm) {
      if (g %in% rownames(avail) && avail[g, nm]) "present" else "**ABSENT**"
    })
    lines <- c(lines, paste0("| ", g, " | ", paste(cells, collapse = " | "), " |"))
  }
  writeLines(lines, file.path(outdir, "profile_plot_warnings.md"))

  invisible(list(grid = grid_res, boxplot = bp_res))
}
