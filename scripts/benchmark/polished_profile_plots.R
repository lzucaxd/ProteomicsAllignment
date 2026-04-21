#!/usr/bin/env Rscript
# =============================================================================
# Polished Marker Sample Profile Plots — Meeting / Paper Quality
# =============================================================================
#
# Design: markers (rows) × methods (columns), max 4 markers per figure page.
# Each cell: faint individual sample points + bold block-median trace.
# Blocks get minimum visual width so even small groups (e.g. CCLE) are readable.
# Y-axis shared across all methods for each marker.
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# ─── Colour definitions ──────────────────────────────────────────────────
BLOCK_COLOURS <- list(
  "CPTAC Luminal" = "#1565C0", "CPTAC Basal"   = "#C62828",
  "CCLE Luminal"  = "#42A5F5", "CCLE Basal"    = "#EF5350",
  "CPTAC Breast"  = "#1565C0", "CPTAC Lung"    = "#2E7D32",
  "CCLE Breast"   = "#42A5F5", "CCLE Lung"     = "#66BB6A"
)
BLOCK_FILLS <- list(
  "CPTAC Luminal" = "#E3F2FD", "CPTAC Basal"   = "#FCE4EC",
  "CCLE Luminal"  = "#F1F8E9", "CCLE Basal"    = "#FFF8E1",
  "CPTAC Breast"  = "#E3F2FD", "CPTAC Lung"    = "#E8F5E9",
  "CCLE Breast"   = "#FFF8E1", "CCLE Lung"     = "#F3E5F5"
)

# Short labels for block annotations (two-line format)
BLOCK_LABELS_SHORT <- list(
  "CPTAC Luminal" = "CPTAC\nLuminal",  "CPTAC Basal"   = "CPTAC\nBasal",
  "CCLE Luminal"  = "CCLE\nLuminal",   "CCLE Basal"    = "CCLE\nBasal",
  "CPTAC Breast"  = "CPTAC\nBreast",   "CPTAC Lung"    = "CPTAC\nLung",
  "CCLE Breast"   = "CCLE\nBreast",    "CCLE Lung"     = "CCLE\nLung"
)

# ─── Build sample order with expanded x-positions ────────────────────────
# Gives each block a minimum visual width so small groups aren't invisible.
build_fixed_sample_order <- function(task_meta, block_order, min_block_width = 12) {
  m <- copy(task_meta)
  m[, block := paste(domain, condition, sep = " ")]
  m[, block := factor(block, levels = block_order)]
  m <- m[order(block, sample_id)]

  # Compute block widths: actual sample count, but at least min_block_width
  block_counts <- m[, .N, by = block]
  setkey(block_counts, block)
  block_counts[, visual_width := pmax(N, min_block_width)]

  # Assign x-positions with gaps between blocks
  gap <- 3
  offset <- 0
  m[, x_pos := NA_real_]
  block_info <- list()
  for (bl in block_order) {
    idx <- which(m$block == bl)
    if (length(idx) == 0) {
      block_info[[bl]] <- list(start = NA, end = NA, mid = NA, n = 0)
      next
    }
    n_bl <- length(idx)
    vis_w <- block_counts[block == bl, visual_width]
    spacing <- vis_w / n_bl
    positions <- offset + (seq_len(n_bl) - 0.5) * spacing
    m$x_pos[idx] <- positions
    bl_start <- offset
    bl_end <- offset + vis_w
    block_info[[bl]] <- list(start = bl_start, end = bl_end,
                              mid = (bl_start + bl_end) / 2, n = n_bl)
    offset <- bl_end + gap
  }
  attr(m, "block_info") <- block_info
  attr(m, "x_max") <- offset - gap
  m
}

# ─── Extract values for one marker from one matrix ───────────────────────
.get_vals <- function(marker, mat, sids) {
  if (!marker %in% rownames(mat)) return(NULL)
  common <- intersect(sids, colnames(mat))
  if (length(common) == 0) return(NULL)
  v <- rep(NA_real_, length(sids))
  names(v) <- sids
  v[common] <- mat[marker, common]
  v
}

# ─── Draw one polished cell ──────────────────────────────────────────────
plot_marker_profile_panel <- function(v, ylim, ordered_meta, block_order,
                                       marker, method_name, is_absent,
                                       show_marker_label, show_method_label,
                                       show_block_labels) {
  block_info <- attr(ordered_meta, "block_info")
  x_max <- attr(ordered_meta, "x_max")
  xpos <- ordered_meta$x_pos
  blocks <- as.character(ordered_meta$block)
  n <- nrow(ordered_meta)

  if (is_absent) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
         xlim = c(-1, x_max + 1), ylim = ylim)
    rect(-1, ylim[1] - 100, x_max + 1, ylim[2] + 100, col = "gray96", border = NA)
    text(x_max / 2, mean(ylim),
         paste0(marker, "\nnot available\nin ", method_name),
         cex = 0.9, col = "gray50", font = 3)
    box(col = "gray85", lwd = 0.6)
    if (show_marker_label)
      mtext(marker, side = 2, line = 0.8, cex = 1.0, font = 2, las = 1)
    if (show_method_label)
      mtext(method_name, side = 3, line = 0.6, cex = 0.95, font = 2, col = "gray20")
    return(invisible(NULL))
  }

  plot(xpos, v, type = "n",
       xlim = c(-1, x_max + 1), ylim = ylim,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n")

  # Y-axis: clean ticks
  pretty_y <- pretty(ylim, n = 5)
  pretty_y <- pretty_y[pretty_y >= ylim[1] & pretty_y <= ylim[2]]
  axis(2, at = pretty_y, labels = formatC(pretty_y, format = "f", digits = 1),
       las = 1, cex.axis = 0.7, tck = -0.02, mgp = c(2, 0.35, 0), col.axis = "gray30")

  # Block shading
  for (bl in block_order) {
    bi <- block_info[[bl]]
    if (is.na(bi$start)) next
    fill <- BLOCK_FILLS[[bl]]
    if (is.null(fill)) fill <- "gray97"
    rect(bi$start, ylim[1] - 100, bi$end, ylim[2] + 100,
         col = fill, border = NA)
  }

  # Block separators (dashed vertical lines between blocks)
  for (i in seq_along(block_order)[-1]) {
    bl <- block_order[i]
    bi <- block_info[[bl]]
    if (is.na(bi$start)) next
    prev_bl <- block_order[i - 1]
    prev_bi <- block_info[[prev_bl]]
    if (is.na(prev_bi$end)) next
    sep_x <- (prev_bi$end + bi$start) / 2
    abline(v = sep_x, col = "gray45", lwd = 1.4, lty = 2)
  }

  # Subtle horizontal grid
  for (yy in pretty_y) {
    abline(h = yy, col = "gray85", lwd = 0.3)
  }

  # Faint individual sample points
  for (bl in block_order) {
    idx <- which(blocks == bl)
    if (length(idx) == 0) next
    col_base <- BLOCK_COLOURS[[bl]]
    if (is.null(col_base)) col_base <- "gray50"
    ok <- !is.na(v[idx])
    points(xpos[idx][ok], v[idx][ok], pch = 16, cex = 0.55,
           col = adjustcolor(col_base, alpha.f = 0.20))
  }

  # Bold block-median trace
  block_medians <- numeric(0)
  block_mids <- numeric(0)
  for (bl in block_order) {
    bi <- block_info[[bl]]
    if (is.na(bi$start)) next
    idx <- which(blocks == bl)
    vals <- v[idx]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) next
    gmed <- median(vals)
    col_base <- BLOCK_COLOURS[[bl]]
    if (is.null(col_base)) col_base <- "gray30"

    segments(bi$start + 0.5, gmed, bi$end - 0.5, gmed,
             col = col_base, lwd = 3.5, lend = 1)

    # Confidence band: IQR
    q25 <- quantile(vals, 0.25)
    q75 <- quantile(vals, 0.75)
    rect(bi$start + 0.5, q25, bi$end - 0.5, q75,
         col = adjustcolor(col_base, alpha.f = 0.08), border = NA)

    block_medians <- c(block_medians, gmed)
    block_mids <- c(block_mids, bi$mid)
  }

  # Connect medians across blocks
  if (length(block_mids) > 1) {
    lines(block_mids, block_medians, col = "gray35", lwd = 1.0, lty = 3)
  }

  box(lwd = 0.5, col = "gray70")

  if (show_marker_label)
    mtext(marker, side = 2, line = 0.8, cex = 1.0, font = 2, las = 1)
  if (show_method_label)
    mtext(method_name, side = 3, line = 0.6, cex = 0.95, font = 2, col = "gray20")

  if (show_block_labels) {
    for (bl in block_order) {
      bi <- block_info[[bl]]
      if (is.na(bi$mid)) next
      lbl <- BLOCK_LABELS_SHORT[[bl]]
      if (is.null(lbl)) lbl <- bl
      col_base <- BLOCK_COLOURS[[bl]]
      if (is.null(col_base)) col_base <- "gray30"
      mtext(lbl, side = 1, line = 1.0, at = bi$mid,
            cex = 0.55, col = col_base, font = 2)
    }
  }
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN: Polished multi-marker figure (3-4 markers per page)
# ═══════════════════════════════════════════════════════════════════════════
make_polished_sample_profile_plots <- function(representations, task_meta,
                                                 markers, task_name, outdir,
                                                 block_order,
                                                 markers_per_page = 4,
                                                 fmt = c("png", "pdf")) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ordered_meta <- build_fixed_sample_order(task_meta, block_order)

  sids <- ordered_meta$sample_id
  x_max <- attr(ordered_meta, "x_max")
  method_names <- names(representations)
  n_methods <- length(method_names)

  # Pre-extract values and compute y-ranges per marker (shared across methods)
  all_vals <- list()
  avail <- list()
  ylims <- list()

  for (g in markers) {
    all_vals[[g]] <- list()
    avail[[g]] <- list()
    finite_vals <- numeric(0)
    for (nm in method_names) {
      v <- .get_vals(g, representations[[nm]]$matrix, sids)
      all_vals[[g]][[nm]] <- v
      avail[[g]][[nm]] <- !is.null(v)
      if (!is.null(v)) finite_vals <- c(finite_vals, v[is.finite(v)])
    }
    if (length(finite_vals) >= 2) {
      rng <- range(finite_vals)
      pad <- diff(rng) * 0.12
      ylims[[g]] <- rng + c(-pad, pad)
    } else {
      ylims[[g]] <- c(-1, 1)
    }
  }

  plottable <- markers[sapply(markers, function(g) any(unlist(avail[[g]])))]
  skipped <- setdiff(markers, plottable)

  if (length(plottable) == 0) {
    cat("  No plottable markers for", task_name, "\n")
    return(invisible(list(plottable = character(0), skipped = markers, paths = character(0))))
  }

  pages <- split(plottable, ceiling(seq_along(plottable) / markers_per_page))
  paths <- character(0)

  for (pi in seq_along(pages)) {
    page_markers <- pages[[pi]]
    n_mk <- length(page_markers)

    # Panel sizing: cap width so large-sample tasks don't produce absurd images
    cell_w <- min(5.5, max(4.0, x_max * 0.035 + 1.0))
    cell_h <- 2.8
    left_margin_in <- 1.6
    fig_w <- left_margin_in + cell_w * n_methods + 0.3
    fig_h <- 0.9 + cell_h * n_mk + 1.0

    do_plot <- function(device_fn, path, ...) {
      device_fn(path, width = fig_w, height = fig_h, ...)

      layout_mat <- matrix(seq_len(n_mk * n_methods),
                            nrow = n_mk, ncol = n_methods, byrow = TRUE)
      layout(layout_mat,
             widths = c(cell_w + left_margin_in - 0.8,
                        rep(cell_w, n_methods - 1)),
             heights = rep(cell_h, n_mk))

      for (ri in seq_along(page_markers)) {
        g <- page_markers[ri]
        for (ci in seq_along(method_names)) {
          nm <- method_names[ci]
          show_marker <- (ci == 1)
          show_method <- (ri == 1)
          show_blocks <- (ri == n_mk)

          left_m <- if (show_marker) 5.5 else 2.5
          top_m  <- if (show_method) 2.2 else 0.6
          bot_m  <- if (show_blocks) 3.0 else 0.6

          par(mar = c(bot_m, left_m, top_m, 0.4), family = "sans")

          v <- all_vals[[g]][[nm]]
          plot_marker_profile_panel(
            v = v, ylim = ylims[[g]], ordered_meta = ordered_meta,
            block_order = block_order, marker = g,
            method_name = representations[[nm]]$name,
            is_absent = is.null(v),
            show_marker_label = show_marker,
            show_method_label = show_method,
            show_block_labels = show_blocks
          )
        }
      }

      dev.off()
    }

    page_suffix <- if (length(pages) > 1) paste0("_p", pi) else ""
    fname <- paste0("profiles_", gsub(" ", "_", task_name), page_suffix)

    if ("png" %in% fmt) {
      p <- file.path(outdir, paste0(fname, ".png"))
      do_plot(png, p, res = 200, units = "in")
      paths <- c(paths, p)
      cat("    ", basename(p), "(", paste(page_markers, collapse = ", "), ")\n")
    }
    if ("pdf" %in% fmt) {
      p <- file.path(outdir, paste0(fname, ".pdf"))
      do_plot(pdf, p)
      paths <- c(paths, p)
    }
  }

  # Warnings / availability log
  warn_lines <- c(
    paste0("# Polished Profile Plot Warnings — ", task_name),
    paste0("\nGenerated: ", Sys.time()),
    paste0("Methods: ", paste(sapply(representations, `[[`, "name"), collapse = ", ")),
    paste0("Markers requested: ", length(markers)),
    paste0("Markers plotted: ", length(plottable)),
    paste0("Markers skipped: ", length(skipped)),
    if (length(skipped) > 0) paste0("  Skipped: ", paste(skipped, collapse = ", ")) else NULL,
    "",
    "## Per-marker Availability",
    ""
  )
  hdr <- paste0("| Marker | ", paste(sapply(representations, `[[`, "name"), collapse = " | "), " |")
  sep_line <- paste0("|", paste(rep("------|", n_methods + 1), collapse = ""))
  warn_lines <- c(warn_lines, hdr, sep_line)
  for (g in markers) {
    cells <- sapply(method_names, function(nm) {
      if (isTRUE(avail[[g]][[nm]])) "present" else "**absent**"
    })
    warn_lines <- c(warn_lines,
      paste0("| ", g, " | ", paste(cells, collapse = " | "), " |"))
  }

  warn_lines <- c(warn_lines, "",
    "## Y-axis Scaling Notes",
    "",
    "For each marker, the y-axis range is shared across all methods to enable",
    "direct visual comparison. If a method produces values on a drastically",
    "different scale (e.g. Celligner z-scores vs Raw log2 abundance), the",
    "shared range accommodates both. Check axis tick labels to compare magnitudes."
  )
  writeLines(warn_lines, file.path(outdir, "polished_profile_plot_warnings.md"))

  invisible(list(plottable = plottable, skipped = skipped, paths = paths))
}
