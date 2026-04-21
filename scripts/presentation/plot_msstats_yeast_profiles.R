#!/usr/bin/env Rscript
# NOTE: Some machines source renv/activate.R from an R profile even when renv
# isn't present in this repo. If you hit that error, run:
#   Rscript --vanilla scripts/presentation/plot_msstats_yeast_profiles.R
#
# MSstatsTMT-style "yeast" profile plots:
# - x: evenly spaced samples within ordered blocks
# - y: log2 abundance (already on transformed matrices)
# - thin line connects consecutive samples WITHIN each block (sorted by abundance)
#
# Outputs:
#   presentation_materials/figures/msstats_profiles/profiles_<task>.{pdf,png}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
sd <- if (length(ff)) {
  dirname(normalizePath(sub("^--file=", "", ff[1L])))
} else {
  normalizePath(file.path(getwd(), "scripts", "presentation"), mustWork = FALSE)
}
source(file.path(sd, "presentation_paths.R"))
pres_ensure_dirs()

OUT_DIR <- file.path(PRES_OUT, "figures", "msstats_profiles")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

METHODS_DIR <- file.path(REPO, "data", "processed", "methods")
UNION_DIR <- file.path(REPO, "data", "processed", "union")

methods <- c("raw", "bridge_shift", "bridge_scale", "celligner")
method_labels <- c(
  raw = "Raw",
  bridge_shift = "Domain shift-only",
  bridge_scale = "Domain shift+scale",
  celligner = "Celligner"
)

markers_subtype <- c("PGR", "GATA3", "FOXA1", "EGFR", "KRT5", "KRT17", "FOXC1")
markers_bvl <- c("NKX2-1", "NAPSA", "SFTPB", "GATA3", "FOXA1", "ESR1")

block_fills <- c(
  "CPTAC\nLuminal" = "#DBEAF9",
  "CCLE\nLuminal"  = "#D9F0D3",
  "CPTAC\nBasal"   = "#FADBD8",
  "CCLE\nBasal"    = "#FEF5C8",
  "CPTAC\nBreast"  = "#DBEAF9",
  "CCLE\nBreast"   = "#D9F0D3",
  "CPTAC\nLung"    = "#FADBD8",
  "CCLE\nLung"     = "#FEF5C8"
)

read_gene_matrix <- function(path) {
  dt <- fread(path)
  if (ncol(dt) < 2L) stop("Matrix has <2 columns: ", path)
  gene_col <- names(dt)[1L]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- as.character(dt[[gene_col]])
  colnames(mat) <- names(dt)[-1L]

  # If genes are columns (wide genes), transpose to genes-as-rows.
  if (!any(rownames(mat) %in% c(markers_subtype, markers_bvl))) {
    maybe_genes <- colnames(mat)
    if (length(maybe_genes) > 0L && any(maybe_genes %in% c(markers_subtype, markers_bvl))) {
      mat <- t(mat)
    }
  }
  mat
}

make_yeast_profile <- function(gene, method, task) {
  mat_path <- file.path(METHODS_DIR, method, paste0("transformed_", task, ".csv"))
  meta_path <- file.path(UNION_DIR, paste0("sample_meta_", task, ".csv"))
  if (!file.exists(mat_path) || !file.exists(meta_path)) return(NULL)

  mat <- read_gene_matrix(mat_path)
  meta <- fread(meta_path)
  if (!all(c("sample_id", "domain", "condition") %in% names(meta))) {
    stop("sample_meta must contain sample_id, domain, condition: ", meta_path)
  }
  meta[, domain := toupper(as.character(domain))]

  if (!gene %in% rownames(mat)) {
    return(
      ggplot() +
        annotate(
          "text", x = 0.5, y = 0.5,
          label = paste0(gene, "\nnot available\nin ", method_labels[[method]]),
          size = 4, color = "grey50", hjust = 0.5
        ) +
        theme_void() +
        theme(plot.background = element_rect(fill = "grey96", color = "grey80"))
    )
  }

  sample_cols <- intersect(colnames(mat), as.character(meta$sample_id))
  if (length(sample_cols) < 4L) return(NULL)

  vals <- as.numeric(mat[gene, sample_cols, drop = TRUE])
  df <- data.table(
    sample_id = sample_cols,
    abundance = vals
  )
  df <- merge(df, meta[, .(sample_id, domain, condition)], by = "sample_id", all.x = TRUE)
  df <- df[is.finite(abundance)]
  df[, block := paste0(domain, "\n", condition)]

  if (identical(task, "breast_subtype")) {
    block_order <- c("CPTAC\nLuminal", "CCLE\nLuminal", "CPTAC\nBasal", "CCLE\nBasal")
  } else if (identical(task, "breast_vs_lung")) {
    block_order <- c("CPTAC\nBreast", "CCLE\nBreast", "CPTAC\nLung", "CCLE\nLung")
  } else {
    stop("Unknown task: ", task)
  }

  df[, block := factor(block, levels = block_order)]
  df <- df[!is.na(block)]
  if (nrow(df) < 4L) return(NULL)

  gap <- 2
  block_counts <- df[, .N, keyby = block]
  setnames(block_counts, "N", "n_block")

  x_start <- 1
  block_ranges <- vector("list", length(block_order))
  names(block_ranges) <- block_order

  for (b in block_order) {
    n <- as.integer(block_counts[block == b, n_block][1L])
    if (is.na(n)) n <- 0L
    if (n > 0L) {
      block_ranges[[b]] <- list(
        start = x_start,
        end = x_start + n - 1L,
        mid = x_start + (n - 1L) / 2
      )
      x_start <- x_start + n + gap
    } else {
      block_ranges[[b]] <- list(start = x_start, end = x_start, mid = x_start)
      x_start <- x_start + gap
    }
  }

  df <- df[order(block, abundance)]
  df[, x_pos := {
    b <- as.character(block[1L])
    r <- block_ranges[[b]]
    n <- .N
    seq(r$start, r$end, length.out = n)
  }, by = block]

  rect_df <- data.table(
    block = factor(block_order, levels = block_order),
    xmin = vapply(block_ranges, function(r) r$start - 0.5, numeric(1)),
    xmax = vapply(block_ranges, function(r) r$end + 0.5, numeric(1))
  )
  rect_df[, fill := block_fills[as.character(block)]]

  y_min <- min(df$abundance, na.rm = TRUE)
  y_max <- max(df$abundance, na.rm = TRUE)
  y_pad <- (y_max - y_min) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  label_df <- data.table(
    block = factor(block_order, levels = block_order),
    x_mid = vapply(block_ranges, function(r) r$mid, numeric(1))
  )

  sep_between_conditions <- (block_ranges[[block_order[2L]]]$end +
    block_ranges[[block_order[3L]]]$start) / 2
  sep_within_lum <- (block_ranges[[block_order[1L]]]$end +
    block_ranges[[block_order[2L]]]$start) / 2
  sep_within_bas <- (block_ranges[[block_order[3L]]]$end +
    block_ranges[[block_order[4L]]]$start) / 2

  ggplot() +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = y_min - y_pad, ymax = y_max + y_pad),
      fill = rect_df$fill, alpha = 0.6, inherit.aes = FALSE
    ) +
    geom_vline(xintercept = sep_between_conditions, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(
      xintercept = c(sep_within_lum, sep_within_bas),
      linetype = "dashed", color = "grey70", linewidth = 0.3
    ) +
    geom_line(
      data = df,
      aes(x = x_pos, y = abundance, group = block),
      color = "grey50", linewidth = 0.3, alpha = 0.6
    ) +
    geom_point(
      data = df,
      aes(x = x_pos, y = abundance),
      color = "grey30", size = 0.8, alpha = 0.7
    ) +
    scale_x_continuous(breaks = label_df$x_mid, labels = as.character(label_df$block)) +
    coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad)) +
    labs(
      title = gene,
      y = "Log2 abundance",
      x = NULL,
      subtitle = paste0(method_labels[[method]], " — ", task)
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle = element_text(size = 8, hjust = 0, color = "grey35"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(size = 7, lineheight = 0.85),
      plot.margin = margin(5, 5, 5, 5)
    )
}

write_task_grid <- function(task, markers) {
  all_plots <- list()
  for (gene in markers) {
    for (method in methods) {
      key <- paste(gene, method, sep = "___")
      p <- tryCatch(
        make_yeast_profile(gene, method, task),
        error = function(e) {
          ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = paste0(gene, "\nerror:\n", e$message), size = 3) +
            theme_void()
        }
      )
      if (is.null(p)) {
        p <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste0(gene, "\nnot available"), size = 3, color = "grey50") +
          theme_void()
      }
      all_plots[[key]] <- p
    }
  }

  n_genes <- length(markers)
  n_methods <- length(methods)

  plot_list_ordered <- vector("list", n_genes * n_methods)
  k <- 0L
  for (gene in markers) {
    for (method in methods) {
      k <- k + 1L
      plot_list_ordered[[k]] <- all_plots[[paste(gene, method, sep = "___")]]
    }
  }

  grid <- plot_grid(plotlist = plot_list_ordered, nrow = n_genes, ncol = n_methods, align = "hv")

  headers <- lapply(methods, function(m) {
    ggdraw() + draw_label(method_labels[[m]], size = 12, fontface = "bold")
  })
  header_row <- plot_grid(plotlist = headers, nrow = 1)

  row_labels <- lapply(markers, function(g) {
    ggdraw() + draw_label(g, size = 11, fontface = "bold", angle = 90)
  })
  label_col <- plot_grid(plotlist = row_labels, ncol = 1)

  main <- plot_grid(header_row, grid, ncol = 1, rel_heights = c(0.04, 0.96))
  final <- plot_grid(label_col, main, nrow = 1, rel_widths = c(0.04, 0.96))

  out_pdf <- file.path(OUT_DIR, sprintf("profiles_%s.pdf", task))
  out_png <- file.path(OUT_DIR, sprintf("profiles_%s.png", task))

  ggsave(out_pdf, final, width = 4.5 * n_methods, height = 2.8 * n_genes)
  ggsave(out_png, final, width = 4.5 * n_methods, height = 2.8 * n_genes, dpi = 150)

  message("Wrote ", out_pdf, " (", n_genes, " genes x ", n_methods, " methods)")
}

write_task_grid("breast_subtype", markers_subtype)
write_task_grid("breast_vs_lung", markers_bvl)

message("Done. Outputs in ", OUT_DIR)
