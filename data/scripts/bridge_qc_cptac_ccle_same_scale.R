#!/usr/bin/env Rscript
# Bridge QC: within each dataset, compare bridge (Norm) channel to all other TMT channels.
# Uses msstats_input.tsv — all PSM rows for per-(Run,Channel) medians; bridge channel
# identified from Condition == "Norm" rows (same channel as MSstatsTMT reference).
# Run-level summaries (bridge-only) remain for cross-run stability.
#
# Heavy step: median log2 by (Run, Channel) over full msstats_input (~40M rows) via DuckDB.
#
# Run from repo root:
#   Rscript --no-init-file data/scripts/bridge_qc_cptac_ccle_same_scale.R
#   Rscript --no-init-file data/scripts/bridge_qc_cptac_ccle_same_scale.R --scale raw
#   Rscript --no-init-file data/scripts/bridge_qc_cptac_ccle_same_scale.R --scale hybrid
#
# --scale log2 (default): both cohorts — median log2(intensity + 1) per run/channel.
# --scale raw: both cohorts — median linear Intensity (CPTAC reporter intensity; CCLE rq_*_sn as Intensity).
#   Pooled y-axis on combined plots; units still differ by cohort.
# --scale hybrid: CPTAC log2(intensity+1); CCLE linear Intensity. Combined plots use facet_wrap(..., scales = "free_y").

suppressPackageStartupMessages({
  if (!requireNamespace("duckdb", quietly = TRUE))
    install.packages("duckdb", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(DBI)
  library(duckdb)
})

repo <- getwd()
if (!file.exists(file.path(repo, "data", "results", "PDC000120", "msstats_input.tsv"))) {
  repo <- normalizePath(file.path(getwd(), ".."))
}
if (!file.exists(file.path(repo, "data", "results", "PDC000120", "msstats_input.tsv"))) {
  repo <- normalizePath(file.path(getwd(), "..", ".."))
}

path_cptac <- file.path(repo, "data", "results", "PDC000120", "msstats_input.tsv")
path_ccle <- file.path(repo, "data", "results", "CCLE_corrected", "msstats_input.tsv")
if (!file.exists(path_ccle)) {
  path_ccle <- file.path(repo, "data", "results", "CCLE", "msstats_input.tsv")
}

if (!file.exists(path_cptac)) stop("Missing CPTAC msstats_input.tsv: ", path_cptac)
if (!file.exists(path_ccle)) stop("Missing CCLE msstats_input.tsv (tried CCLE_corrected then CCLE): ", path_ccle)

out_dir <- file.path(repo, "reports", "bridge_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

args_all <- commandArgs(trailingOnly = TRUE)
scale_mode <- if ("--scale" %in% args_all) {
  i <- match("--scale", args_all)
  v <- if (length(args_all) >= i + 1L) tolower(trimws(args_all[i + 1L])) else "log2"
  if (!v %in% c("log2", "raw", "hybrid")) stop("--scale must be log2, raw, or hybrid, got: ", v)
  v
} else {
  "log2"
}
message("Scale mode: ", scale_mode)

log2_int <- function(x) {
  log2(pmax(as.numeric(x), 0, na.rm = TRUE) + 1)
}

# Fast: which TMT channel carries Norm/bridge (typically 131)
detect_bridge_channel <- function(msstats_path) {
  awk_cmd <- sprintf(
    "awk -F'\\t' 'NR>1 && $8==\"Norm\" {print $7}' %s | sort -u",
    shQuote(msstats_path, type = "sh")
  )
  ch <- tryCatch(
    trimws(system(awk_cmd, intern = TRUE)),
    error = function(e) character(0)
  )
  ch <- unique(ch[nzchar(ch)])
  if (length(ch) >= 1L) return(as.character(ch[[1]]))
  "131"
}

load_norm_long <- function(msstats_path) {
  awk_cmd <- sprintf(
    "awk -F'\\t' 'NR==1 || $8==\"Norm\"' %s",
    shQuote(msstats_path, type = "sh")
  )
  message("Loading Norm rows (bridge only): ", basename(msstats_path), " ...")
  bn <- fread(
    cmd = awk_cmd,
    sep = "\t",
    quote = "",
    showProgress = TRUE,
    select = c("Run", "Channel", "Intensity")
  )
  if (nrow(bn) == 0L) stop("No Norm rows in ", msstats_path)
  bn[, Channel := as.character(Channel)]
  bn[, Intensity := as.numeric(Intensity)]
  bn[, log2_intensity := log2_int(Intensity)]
  bn
}

summarize_bridge_run <- function(bn, label, scale_mode) {
  use_raw <- identical(scale_mode, "raw") ||
    (identical(scale_mode, "hybrid") && label != "CPTAC")
  if (use_raw) {
    bn[, .(
      median_log2_intensity = as.numeric(median(pmax(Intensity, 0, na.rm = TRUE), na.rm = TRUE)),
      n_psm = .N
    ), by = Run][, dataset := label]
  } else {
    bn[, .(
      median_log2_intensity = as.numeric(median(log2_intensity, na.rm = TRUE)),
      n_psm = .N
    ), by = Run][, dataset := label]
  }
}

#' Median log2(Intensity+1) per (Run, Channel) over ALL PSM rows in msstats_input (DuckDB).
aggregate_all_channels_duckdb <- function(msstats_path, label) {
  message("DuckDB: per-(Run,Channel) medians over full file: ", basename(msstats_path), " ...")
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  path_lit <- dbQuoteString(con, normalizePath(msstats_path, winslash = "/", mustWork = TRUE))
  sql <- paste0(
    "SELECT \"Run\" AS run, CAST(\"Channel\" AS VARCHAR) AS channel, ",
    "median(log2(1 + greatest(COALESCE(TRY_CAST(\"Intensity\" AS DOUBLE), 0), 0))) AS median_log2_intensity, ",
    "COUNT(*) AS n_psm ",
    "FROM read_csv_auto(", path_lit, ", delim='\t', header=true, quote='', nullstr='') ",
    "GROUP BY \"Run\", \"Channel\""
  )
  d <- as.data.table(dbGetQuery(con, sql))
  setnames(d, c("Run", "Channel", "median_log2_intensity", "n_psm"))
  d[, dataset := label]
  d
}

#' Median linear Intensity per (Run, Channel) — DuckDB over full file.
aggregate_all_channels_duckdb_raw <- function(msstats_path, label) {
  message("DuckDB (raw linear): per-(Run,Channel) medians: ", basename(msstats_path), " ...")
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  path_lit <- dbQuoteString(con, normalizePath(msstats_path, winslash = "/", mustWork = TRUE))
  sql <- paste0(
    "SELECT \"Run\" AS run, CAST(\"Channel\" AS VARCHAR) AS channel, ",
    "median(greatest(COALESCE(TRY_CAST(\"Intensity\" AS DOUBLE), 0), 0)) AS median_log2_intensity, ",
    "COUNT(*) AS n_psm ",
    "FROM read_csv_auto(", path_lit, ", delim='\t', header=true, quote='', nullstr='') ",
    "GROUP BY \"Run\", \"Channel\""
  )
  d <- as.data.table(dbGetQuery(con, sql))
  setnames(d, c("Run", "Channel", "median_log2_intensity", "n_psm"))
  d[, dataset := label]
  d
}

tmt_channel_order <- function(ch) {
  tmt <- c(
    "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C",
    "131", "132N", "132C", "133N", "133C", "134N"
  )
  u <- unique(as.character(ch))
  c(intersect(tmt, u), sort(setdiff(u, tmt)))
}

robust_ylim <- function(x, lo = 0.01, hi = 0.99, pad = 0.04) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(0, 1))
  q <- as.numeric(stats::quantile(x, probs = c(lo, hi), na.rm = TRUE, names = FALSE))
  span <- diff(q)
  if (!is.finite(span) || span < 1e-8) span <- max(0.05, stats::mad(x, na.rm = TRUE) * 3)
  c(q[1] - pad * span, q[2] + pad * span)
}

plot_run_box <- function(run_dt, title, ylim, subtitle = NULL, scale_mode = "log2") {
  ylab <- if (identical(scale_mode, "raw")) {
    "Median intensity (linear)\nper run (Norm PSMs only)"
  } else {
    "Median log2(intensity + 1)\nper run (Norm PSMs only)"
  }
  ggplot(run_dt, aes(x = "Bridge (Norm) per run", y = median_log2_intensity)) +
    geom_boxplot(width = 0.35, fill = "grey88", colour = "grey35", outlier.size = 1.2) +
    geom_jitter(width = 0.08, height = 0, size = 0.9, alpha = 0.35, colour = "grey25") +
    coord_cartesian(ylim = ylim) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = ylab
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
}

#' Boxplots: one box per TMT channel; bridge channel highlighted vs sample channels.
plot_channel_bridge_vs_others <- function(ch_dt, bridge_ch, title, ylim, subtitle = NULL, scale_mode = "log2") {
  ylab <- if (identical(scale_mode, "raw")) {
    "Median intensity (linear)\nper run and channel (all PSMs)"
  } else {
    "Median log2(intensity + 1)\nper run and channel (all PSMs)"
  }
  ch_dt <- copy(ch_dt)
  ch_dt[, Channel := as.character(Channel)]
  br <- as.character(bridge_ch)
  br_lab <- paste0("Bridge (ch ", br, ")")
  ch_dt[, role := ifelse(Channel == br, br_lab, "Sample channels")]
  ch_dt[, Channel := factor(Channel, levels = tmt_channel_order(ch_dt$Channel))]
  fills <- c("steelblue3", "grey92")
  names(fills) <- c(br_lab, "Sample channels")
  ggplot(ch_dt, aes(x = Channel, y = median_log2_intensity, fill = role)) +
    geom_boxplot(width = 0.65, alpha = 0.92, colour = "grey35", outlier.size = 0.65) +
    geom_jitter(width = 0.1, size = 0.35, alpha = 0.2, colour = "grey25") +
    coord_cartesian(ylim = ylim) +
    scale_fill_manual(name = NULL, values = fills) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "TMT channel (all channels in plex)",
      y = ylab
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}

plot_combined_runs <- function(run_cptac, run_ccle, ylim, scale_mode = "log2") {
  sub <- if (identical(scale_mode, "raw")) {
    "Each point = one MS run; median linear intensity over Norm PSM rows only"
  } else {
    "Each point = one MS run; median log2 over Norm PSM rows only"
  }
  ylab <- if (identical(scale_mode, "raw")) {
    "Median intensity (linear) per run"
  } else {
    "Median log2(intensity + 1) per run"
  }
  d <- rbindlist(list(
    run_cptac[, .(dataset, median_log2_intensity)],
    run_ccle[, .(dataset, median_log2_intensity)]
  ))
  ggplot(d, aes(x = dataset, y = median_log2_intensity)) +
    geom_boxplot(width = 0.45, fill = "grey90", colour = "grey30", outlier.size = 1) +
    geom_jitter(width = 0.07, size = 0.7, alpha = 0.3, colour = "grey25") +
    coord_cartesian(ylim = ylim) +
    labs(
      title = "Bridge (Norm) — run-level median spread (same y scale)",
      subtitle = sub,
      x = NULL,
      y = ylab
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

#' Combined bridge run spread: CPTAC log2 vs CCLE linear; **free y** per cohort.
plot_combined_runs_hybrid <- function(run_cptac, run_ccle) {
  d <- rbindlist(list(
    run_cptac[, .(dataset, median_log2_intensity)],
    run_ccle[, .(dataset, median_log2_intensity)]
  ))
  ggplot(d, aes(x = "Norm / bridge", y = median_log2_intensity)) +
    geom_boxplot(width = 0.35, fill = "grey90", colour = "grey30", outlier.size = 1) +
    geom_jitter(width = 0.06, size = 0.7, alpha = 0.35, colour = "grey25") +
    facet_wrap(~dataset, scales = "free_y", ncol = 2L) +
    labs(
      title = "Bridge (Norm) — run-level medians (CPTAC log2, CCLE linear)",
      subtitle = "Left: median log2(Intensity+1); right: median linear Intensity — separate y scales",
      x = NULL,
      y = "Median (scale differs by panel)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

#' All channels, CPTAC log2 vs CCLE linear; **free y** per cohort.
plot_combined_channels_hybrid <- function(ch_cptac, ch_ccle, bridge_cpt, bridge_cc, cc_dataset_label) {
  d <- rbindlist(list(
    ch_cptac[, .(dataset, Channel, median_log2_intensity)],
    ch_ccle[, .(dataset, Channel, median_log2_intensity)]
  ))
  d[, Channel := as.character(Channel)]
  d[, role := ifelse(
    (dataset == "CPTAC" & Channel == bridge_cpt) |
      (dataset == cc_dataset_label & Channel == bridge_cc),
    "Bridge", "Sample"
  )]
  d[, Channel := factor(Channel, levels = tmt_channel_order(d$Channel))]
  d[, grp := paste(dataset, role, sep = ".")]
  fill_cols <- c(
    "steelblue3", "grey85", "darkseagreen3", "wheat"
  )
  names(fill_cols) <- c(
    "CPTAC.Bridge", "CPTAC.Sample",
    paste0(cc_dataset_label, ".Bridge"),
    paste0(cc_dataset_label, ".Sample")
  )
  ggplot(d, aes(x = Channel, y = median_log2_intensity, fill = grp)) +
    geom_boxplot(position = position_dodge(width = 0.88), width = 0.82, outlier.size = 0.45, alpha = 0.9, colour = "grey30") +
    facet_wrap(~dataset, scales = "free_y", ncol = 2L) +
    scale_fill_manual(name = NULL, values = fill_cols) +
    labs(
      title = "All TMT channels — CPTAC log2 vs CCLE linear (free y per cohort)",
      subtitle = "CPTAC: median log2(Intensity+1); CCLE: median linear Intensity — compare shape within each panel",
      x = "Channel",
      y = "Median (scale differs by panel)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}

plot_combined_channels <- function(ch_cptac, ch_ccle, bridge_cpt, bridge_cc, ylim, cc_dataset_label, scale_mode = "log2") {
  sub <- if (identical(scale_mode, "raw")) {
    "Per (run, channel): median linear PSM intensity; bridge channel filled darker per cohort | CPTAC vs CCLE units differ"
  } else {
    "Per (run, channel): median log2 PSM intensity; bridge channel filled darker per cohort"
  }
  ylab <- if (identical(scale_mode, "raw")) {
    "Median intensity (linear)"
  } else {
    "Median log2(intensity + 1)"
  }
  d <- rbindlist(list(
    ch_cptac[, .(dataset, Channel, median_log2_intensity)],
    ch_ccle[, .(dataset, Channel, median_log2_intensity)]
  ))
  d[, Channel := as.character(Channel)]
  d[, role := ifelse(
    (dataset == "CPTAC" & Channel == bridge_cpt) |
      (dataset == cc_dataset_label & Channel == bridge_cc),
    "Bridge", "Sample"
  )]
  d[, Channel := factor(Channel, levels = tmt_channel_order(d$Channel))]
  d[, grp := paste(dataset, role, sep = ".")]
  fill_cols <- c(
    "steelblue3", "grey85", "darkseagreen3", "wheat"
  )
  names(fill_cols) <- c(
    "CPTAC.Bridge", "CPTAC.Sample",
    paste0(cc_dataset_label, ".Bridge"),
    paste0(cc_dataset_label, ".Sample")
  )
  ggplot(d, aes(x = Channel, y = median_log2_intensity, fill = grp)) +
    geom_boxplot(position = position_dodge(width = 0.88), width = 0.82, outlier.size = 0.45, alpha = 0.9, colour = "grey30") +
    coord_cartesian(ylim = ylim) +
    scale_fill_manual(name = NULL, values = fill_cols) +
    labs(
      title = "All TMT channels — within-plex profile (same y scale)",
      subtitle = sub,
      x = "Channel",
      y = ylab
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}

# ---- Bridge channel labels ----
bridge_cpt <- detect_bridge_channel(path_cptac)
bridge_cc <- detect_bridge_channel(path_ccle)
message("Bridge channel (Norm): CPTAC=", bridge_cpt, " CCLE=", bridge_cc)

# ---- Norm-only: run stability ----
bn_cpt <- load_norm_long(path_cptac)
bn_ccle <- load_norm_long(path_ccle)

lab_cc <- if (grepl("CCLE_corrected", path_ccle)) "CCLE_corrected" else "CCLE"
run_cpt <- summarize_bridge_run(bn_cpt, "CPTAC", scale_mode)
run_cc <- summarize_bridge_run(bn_ccle, lab_cc, scale_mode)

suffix <- if (identical(scale_mode, "raw")) {
  "_raw"
} else if (identical(scale_mode, "hybrid")) {
  "_hybrid"
} else {
  ""
}

fwrite(run_cpt, file.path(out_dir, paste0("cptac_bridge_run_summary", suffix, ".tsv")), sep = "\t")
fwrite(run_cc, file.path(out_dir, paste0("ccle_bridge_run_summary", suffix, ".tsv")), sep = "\t")

# ---- All channels: DuckDB ----
if (identical(scale_mode, "raw")) {
  ch_cpt <- aggregate_all_channels_duckdb_raw(path_cptac, "CPTAC")
  ch_cc <- aggregate_all_channels_duckdb_raw(path_ccle, lab_cc)
} else if (identical(scale_mode, "hybrid")) {
  ch_cpt <- aggregate_all_channels_duckdb(path_cptac, "CPTAC")
  ch_cc <- aggregate_all_channels_duckdb_raw(path_ccle, lab_cc)
} else {
  ch_cpt <- aggregate_all_channels_duckdb(path_cptac, "CPTAC")
  ch_cc <- aggregate_all_channels_duckdb(path_ccle, lab_cc)
}

fwrite(ch_cpt, file.path(out_dir, paste0("cptac_run_channel_medians_all_channels", suffix, ".tsv")), sep = "\t")
fwrite(ch_cc, file.path(out_dir, paste0("ccle_run_channel_medians_all_channels", suffix, ".tsv")), sep = "\t")

# Legacy filenames (now = all-channel table)
fwrite(ch_cpt, file.path(out_dir, paste0("cptac_bridge_channel_summary", suffix, ".tsv")), sep = "\t")
fwrite(ch_cc, file.path(out_dir, paste0("ccle_bridge_channel_summary", suffix, ".tsv")), sep = "\t")

# ---- Shared y limits ----
if (identical(scale_mode, "hybrid")) {
  ylim_run_cpt <- robust_ylim(run_cpt$median_log2_intensity)
  ylim_run_cc <- robust_ylim(run_cc$median_log2_intensity)
  ylim_ch_cpt <- robust_ylim(ch_cpt$median_log2_intensity)
  ylim_ch_cc <- robust_ylim(ch_cc$median_log2_intensity)
  ylim_run <- range(c(ylim_run_cpt, ylim_run_cc), na.rm = TRUE)
  ylim_ch <- range(c(ylim_ch_cpt, ylim_ch_cc), na.rm = TRUE)
} else {
  ylim_run <- robust_ylim(c(run_cpt$median_log2_intensity, run_cc$median_log2_intensity))
  ylim_ch <- robust_ylim(c(ch_cpt$median_log2_intensity, ch_cc$median_log2_intensity))
}

ccle_plot_label <- if (grepl("CCLE_corrected", path_ccle)) "CCLE (corrected)" else "CCLE"

note <- c(
  "Shared y-axis limits for bridge QC plots",
  "========================================",
  "",
  paste("Scale mode:", scale_mode),
  if (identical(scale_mode, "raw")) {
    c(
      "",
      "RAW LINEAR: median(Intensity) per run or per (run, channel).",
      "  CPTAC Intensity = TMT reporter intensity from PSM table.",
      "  CCLE Intensity = reporter-ion S/N (rq_*_sn) ingested as MSstats Intensity.",
      "  Pooled y-limits are for visualization only; absolute units are NOT comparable across cohorts."
    )
  } else if (identical(scale_mode, "hybrid")) {
    c(
      "",
      "HYBRID: CPTAC median log2(Intensity+1); CCLE median linear Intensity.",
      "  Combined channel/run figures: facet_wrap with free_y (separate y per cohort).",
      "  Per-cohort PNGs: CPTAC uses log2 y-axis; CCLE uses linear y-axis."
    )
  } else {
    c(
      "",
      "LOG2: median(log2(Intensity+1)) per run or per (run, channel)."
    )
  },
  "",
  "msstats_input files:",
  paste("  CPTAC:", path_cptac),
  paste("  CCLE: ", path_ccle),
  "",
  "Bridge channel (from Condition == Norm rows):",
  paste("  CPTAC:", bridge_cpt),
  paste("  CCLE: ", bridge_cc),
  "",
  "Run-level plots:",
  "  Source: Norm/bridge PSM rows only (Condition == Norm).",
  if (identical(scale_mode, "raw")) {
    "  y = median linear Intensity pooled within each run."
  } else if (identical(scale_mode, "hybrid")) {
    "  y = CPTAC: median log2(Intensity+1); CCLE: median linear Intensity per run."
  } else {
    "  y = median log2(Intensity+1) pooled within each run."
  },
  "",
  "Channel-level plots:",
  "  Source: ALL PSM rows in msstats_input (every TMT channel).",
  if (identical(scale_mode, "raw")) {
    "  Per (Run, Channel): median linear Intensity."
  } else if (identical(scale_mode, "hybrid")) {
    "  Per (Run, Channel): CPTAC log2 medians; CCLE linear medians."
  } else {
    "  Per (Run, Channel): median log2(Intensity+1)."
  },
  "  Bridge channel is highlighted vs sample channels (same plex, within dataset).",
  "",
  "Y-axis rule:",
  if (identical(scale_mode, "hybrid")) {
    c(
      "  Hybrid: separate robust_ylim per cohort for single-cohort PNGs; combined uses free_y facets.",
      sprintf("  CPTAC run ylim: [%.6f, %.6f]", ylim_run_cpt[1], ylim_run_cpt[2]),
      sprintf("  CCLE run ylim:  [%.6f, %.6f]", ylim_run_cc[1], ylim_run_cc[2]),
      sprintf("  CPTAC channel ylim: [%.6f, %.6f]", ylim_ch_cpt[1], ylim_ch_cpt[2]),
      sprintf("  CCLE channel ylim:  [%.6f, %.6f]", ylim_ch_cc[1], ylim_ch_cc[2])
    )
  } else {
    c(
      "  1st–99th percentile + 4% padding (see robust_ylim in script), cohorts pooled.",
      sprintf("Run plots ylim: [%.6f, %.6f]", ylim_run[1], ylim_run[2]),
      sprintf("Channel plots ylim: [%.6f, %.6f]", ylim_ch[1], ylim_ch[2])
    )
  },
  "",
  "Aggregation: DuckDB median() over full TSV (GROUP BY Run, Channel).",
  ""
)
writeLines(note, file.path(out_dir, paste0("shared_y_axis_note", suffix, ".txt")))

# ---- Plots ----
if (identical(scale_mode, "hybrid")) {
  p_cpt_run <- plot_run_box(
    run_cpt,
    "CPTAC — bridge (Norm) run medians",
    ylim_run_cpt,
    subtitle = "Hybrid mode: log2 scale (see combined figure for CCLE linear)",
    scale_mode = "log2"
  )
  ggsave(file.path(out_dir, paste0("cptac_bridge_run_boxplot_same_scale", suffix, ".png")), p_cpt_run, width = 6.5, height = 4.9, dpi = 200)

  p_cc_run <- plot_run_box(
    run_cc,
    paste(ccle_plot_label, "— bridge (Norm) run medians"),
    ylim_run_cc,
    subtitle = "Hybrid mode: linear S/N scale (see combined figure for CPTAC log2)",
    scale_mode = "raw"
  )
  ggsave(file.path(out_dir, paste0("ccle_bridge_run_boxplot_same_scale", suffix, ".png")), p_cc_run, width = 6.5, height = 4.9, dpi = 200)

  p_cpt_ch <- plot_channel_bridge_vs_others(
    ch_cpt,
    bridge_cpt,
    "CPTAC — bridge vs all TMT channels (within study)",
    ylim_ch_cpt,
    subtitle = "log2 medians per run × channel (hybrid: CPTAC only)",
    scale_mode = "log2"
  )
  ggsave(file.path(out_dir, paste0("cptac_bridge_channel_boxplot_same_scale", suffix, ".png")), p_cpt_ch, width = 10, height = 5.8, dpi = 200)

  p_cc_ch <- plot_channel_bridge_vs_others(
    ch_cc,
    bridge_cc,
    paste(ccle_plot_label, "— bridge vs all TMT channels (within study)"),
    ylim_ch_cc,
    subtitle = "linear medians per run × channel (hybrid: CCLE only)",
    scale_mode = "raw"
  )
  ggsave(file.path(out_dir, paste0("ccle_bridge_channel_boxplot_same_scale", suffix, ".png")), p_cc_ch, width = 10, height = 5.8, dpi = 200)

  ggsave(
    file.path(out_dir, paste0("bridge_run_boxplot_combined", suffix, ".png")),
    plot_combined_runs_hybrid(run_cpt, run_cc),
    width = 8.5,
    height = 4.9,
    dpi = 200
  )

  ggsave(
    file.path(out_dir, paste0("bridge_channel_boxplot_combined", suffix, ".png")),
    plot_combined_channels_hybrid(ch_cpt, ch_cc, bridge_cpt, bridge_cc, lab_cc),
    width = 12,
    height = 6,
    dpi = 200
  )
} else {
  p_cpt_run <- plot_run_box(
    run_cpt,
    "CPTAC — bridge (Norm) run medians",
    ylim_run,
    subtitle = "Same y scale as CCLE run plot",
    scale_mode = scale_mode
  )
  ggsave(file.path(out_dir, paste0("cptac_bridge_run_boxplot_same_scale", suffix, ".png")), p_cpt_run, width = 6.5, height = 4.9, dpi = 200)

  p_cc_run <- plot_run_box(
    run_cc,
    paste(ccle_plot_label, "— bridge (Norm) run medians"),
    ylim_run,
    subtitle = "Same y scale as CPTAC run plot",
    scale_mode = scale_mode
  )
  ggsave(file.path(out_dir, paste0("ccle_bridge_run_boxplot_same_scale", suffix, ".png")), p_cc_run, width = 6.5, height = 4.9, dpi = 200)

  p_cpt_ch <- plot_channel_bridge_vs_others(
    ch_cpt,
    bridge_cpt,
    "CPTAC — bridge vs all TMT channels (within study)",
    ylim_ch,
    subtitle = "Medians per run x channel; bridge channel (Norm) vs sample channels | same y as CCLE",
    scale_mode = scale_mode
  )
  ggsave(file.path(out_dir, paste0("cptac_bridge_channel_boxplot_same_scale", suffix, ".png")), p_cpt_ch, width = 10, height = 5.8, dpi = 200)

  p_cc_ch <- plot_channel_bridge_vs_others(
    ch_cc,
    bridge_cc,
    paste(ccle_plot_label, "— bridge vs all TMT channels (within study)"),
    ylim_ch,
    subtitle = "Same y scale as CPTAC channel plot",
    scale_mode = scale_mode
  )
  ggsave(file.path(out_dir, paste0("ccle_bridge_channel_boxplot_same_scale", suffix, ".png")), p_cc_ch, width = 10, height = 5.8, dpi = 200)

  ggsave(
    file.path(out_dir, paste0("bridge_run_boxplot_combined", suffix, ".png")),
    plot_combined_runs(run_cpt, run_cc, ylim_run, scale_mode = scale_mode),
    width = 6,
    height = 4.9,
    dpi = 200
  )

  ggsave(
    file.path(out_dir, paste0("bridge_channel_boxplot_combined", suffix, ".png")),
    plot_combined_channels(ch_cpt, ch_cc, bridge_cpt, bridge_cc, ylim_ch, lab_cc, scale_mode = scale_mode),
    width = 12,
    height = 6,
    dpi = 200
  )
}

interp <- c(
  "Bridge vs all channels — interpretation",
  "=====================================",
  "",
  "What this shows (within each dataset)",
  "- msstats_input.tsv lists PSM-level reporter intensities for every TMT channel in the plex.",
  "- The bridge/reference is the channel where Condition == \"Norm\" (pooled reference in CPTAC; sample-sheet Norm in CCLE).",
  "- For each LC-MS run and each channel, we take the median of log2(Intensity+1) across all PSMs in that run and channel.",
  "- The channel boxplots compare the bridge channel to all other (sample) channels in the same study: typical TMT imbalance and how the reference sits relative to samples.",
  "",
  "Why not protein_summary.tsv",
  "- After summarization, the reference channel is usually not reported as a separate protein column; use msstats_input for channel-level PSM data.",
  "",
  "Run boxplots (still Norm-only)",
  "- One value per run: median log2 over Norm PSM rows only — cross-run stability of the bridge.",
  "",
  "Cross-study (CPTAC vs CCLE) same y-axis",
  "- Matched limits help compare spread/outlier behavior, not absolute calibration.",
  "- Do not interpret vertical offsets between cohorts as biological truth.",
  "",
  "Files: reports/bridge_qc/*.tsv, *.png, shared_y_axis_note.txt",
  ""
)
writeLines(interp, file.path(out_dir, paste0("bridge_qc_same_scale_interpretation", suffix, ".txt")))

message("Done. Outputs in ", out_dir)
