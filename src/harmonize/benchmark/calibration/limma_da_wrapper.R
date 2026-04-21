#!/usr/bin/env Rscript
# =============================================================================
# limma DA Wrapper — CLI-callable per-domain limma for the Python benchmark
# =============================================================================
#
# Contrast direction (explicit via --contrast-a and --contrast-b):
#   makeContrasts( paste0(contrast_b, " - ", contrast_a), ... )
#   with model.matrix(~ 0 + group_factor) and
#   group_factor <- factor(..., levels = c(contrast_a, contrast_b)).
#
# Breast subtype (defaults from run_all_limma_da.R):
#   contrast_a = Basal, contrast_b = Luminal, name = Luminal_vs_Basal
#   => logFC = mean(Luminal) - mean(Basal). Positive logFC => higher in Luminal.
#
# Breast vs lung:
#   contrast_a = Breast, contrast_b = Lung => logFC = Lung - Breast.
#
# Outputs (same rows/columns in both paths):
#   <outdir>/<domain>/da_limma_result.csv   — canonical
#   <outdir>/da_<domain>.csv                — flat copy (do not use stale copies)
#
# Usage:
#   Rscript limma_da_wrapper.R \
#     --matrix input_matrix.csv \
#     --meta   input_meta.csv \
#     --contrast-a Basal \
#     --contrast-b Luminal \
#     --outdir /path/to/output
#
# Runs limma separately for each domain (CPTAC, CCLE) found in metadata.
# Contrast direction is always contrast_b - contrast_a (explicit, not
# dependent on metadata ordering).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ── Parse CLI args ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  args[idx + 1]
}

matrix_path   <- parse_arg("--matrix")
meta_path     <- parse_arg("--meta")
contrast_a    <- parse_arg("--contrast-a")
contrast_b    <- parse_arg("--contrast-b")
outdir        <- parse_arg("--outdir")
contrast_name <- parse_arg("--contrast-name")

if (is.null(matrix_path) || is.null(meta_path) || is.null(outdir))
  stop("Required: --matrix, --meta, --outdir")
if (is.null(contrast_a) || is.null(contrast_b))
  stop("Required: --contrast-a and --contrast-b (explicit contrast levels)")

if (is.null(contrast_name))
  contrast_name <- paste0(contrast_b, "_vs_", contrast_a)

# ── Load data ────────────────────────────────────────────────────────────────
mat  <- as.matrix(fread(matrix_path), rownames = 1)
meta <- fread(meta_path)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("limma DA wrapper\n")
cat("  Matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")
cat("  Contrast:", contrast_b, "-", contrast_a, "\n")

# ── Per-domain limma ─────────────────────────────────────────────────────────
run_domain_limma <- function(mat, meta, domain, contrast_a, contrast_b,
                              contrast_name, outdir) {
  dom_sids <- meta[toupper(domain_col) == toupper(domain), sample_id]
  dom_cols <- intersect(dom_sids, colnames(mat))

  if (length(dom_cols) < 4) {
    cat("  ", domain, ": only", length(dom_cols), "samples — skipping\n")
    return(NULL)
  }

  dom_mat    <- mat[, dom_cols, drop = FALSE]
  dom_groups <- meta[match(dom_cols, sample_id), condition]

  uniq <- unique(dom_groups)
  if (!contrast_a %in% uniq || !contrast_b %in% uniq) {
    cat("  ", domain, ": condition levels", paste(uniq, collapse = ", "),
        "do not contain both", contrast_a, "and", contrast_b, "— skipping\n")
    return(NULL)
  }

  # Keep only samples from the two contrast levels
  keep_idx <- dom_groups %in% c(contrast_a, contrast_b)
  dom_mat    <- dom_mat[, keep_idx, drop = FALSE]
  dom_groups <- dom_groups[keep_idx]

  # Filter genes with >50% NA in either group
  for (g in c(contrast_a, contrast_b)) {
    g_cols <- which(dom_groups == g)
    if (length(g_cols) > 0) {
      na_frac <- rowMeans(is.na(dom_mat[, g_cols, drop = FALSE]))
      dom_mat[na_frac > 0.5, g_cols] <- NA
    }
  }
  keep_genes <- rowSums(!is.na(dom_mat)) >= ncol(dom_mat) * 0.5
  dom_mat <- dom_mat[keep_genes, , drop = FALSE]

  if (nrow(dom_mat) < 10) {
    cat("  ", domain, ": only", nrow(dom_mat), "genes after filtering — skipping\n")
    return(NULL)
  }

  cat("  ", domain, ":", ncol(dom_mat), "samples,", nrow(dom_mat), "genes\n")

  # Explicit factor with controlled level order so contrast_b - contrast_a works
  group_factor <- factor(dom_groups, levels = c(contrast_a, contrast_b))
  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- levels(group_factor)

  fit <- lmFit(dom_mat, design)
  contrast_str <- paste0(contrast_b, " - ", contrast_a)
  cm <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)

  tt <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"))
  tt$gene <- rownames(topTable(fit2, number = Inf, sort.by = "none"))
  tt$contrast <- contrast_name
  tt$domain <- domain
  tt$inference_type <- "representation_level_limma"
  setcolorder(tt, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
                     "B", "contrast", "domain", "inference_type"))

  dom_dir <- file.path(outdir, tolower(domain))
  dir.create(dom_dir, recursive = TRUE, showWarnings = FALSE)
  sub_path <- file.path(dom_dir, "da_limma_result.csv")
  fwrite(tt, sub_path)
  cat("    Saved:", sub_path, "\n")
  # Flat alias at representation_da/da_<domain>.csv (identical to subfolder file)
  flat <- file.path(outdir, paste0("da_", tolower(domain), ".csv"))
  fwrite(tt, flat)
  cat("    Saved:", flat, "\n")

  tt
}

# Rename domain column for consistent access
if ("domain" %in% names(meta)) {
  setnames(meta, "domain", "domain_col")
} else {
  stop("Metadata must have a 'domain' column")
}

domains_present <- unique(toupper(meta$domain_col))
cat("  Domains found:", paste(domains_present, collapse = ", "), "\n")

for (dom in c("CPTAC", "CCLE")) {
  if (toupper(dom) %in% toupper(domains_present)) {
    run_domain_limma(mat, meta, dom, contrast_a, contrast_b, contrast_name, outdir)
  }
}

cat("limma DA wrapper complete.\n")
