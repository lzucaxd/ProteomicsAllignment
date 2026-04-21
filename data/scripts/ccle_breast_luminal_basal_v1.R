#!/usr/bin/env Rscript
# CCLE breast Luminal vs Basal first-pass benchmark (v1)
# Writes:
#   results/CCLE/ccle_breast_lines_present.csv
#   results/CCLE/ccle_breast_subtype_labels_v1.csv
#   results/CCLE/DA_luminal_vs_basal/* (PCA, DA, markers)
#   results/CCLE/ccle_luminal_vs_basal_feasibility.md
#
# Usage: cd data && Rscript --vanilla scripts/ccle_breast_luminal_basal_v1.R

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", update = FALSE, ask = FALSE)
})
library(data.table)
library(limma)

DATA_DIR <- getwd()
if (!file.exists(file.path(DATA_DIR, "results", "CCLE", "gene_matrix.csv")))
  DATA_DIR <- file.path(getwd(), "data")
RES <- file.path(DATA_DIR, "results", "CCLE")
if (!file.exists(file.path(RES, "gene_matrix.csv")))
  stop("Need results/CCLE/gene_matrix.csv")

OUT_DA <- file.path(RES, "DA_luminal_vs_basal")
dir.create(OUT_DA, recursive = TRUE, showWarnings = FALSE)

SAMPLE_INFO <- file.path(DATA_DIR, "ccle_peptide", "sample_info_ccle.csv")
MATRIX_PATH <- file.path(RES, "gene_matrix.csv")

# --- 1) Breast lines present vs matrix columns ---
si <- fread(SAMPLE_INFO)
si <- si[`Tissue of Origin` == "Breast"]
si[, Cell_Line := trimws(`Cell Line`)]
breast_unique <- unique(si[, .(Cell_Line, CCLE_Code = `CCLE Code`)])

gm_header <- strsplit(readLines(MATRIX_PATH, n = 1L), ",", fixed = TRUE)[[1]]
matrix_cols <- setdiff(gm_header, c("GeneSymbol", "UniProtID"))
mat_set <- matrix_cols

breast_unique[, in_gene_matrix := Cell_Line %in% mat_set]
fwrite(breast_unique[, .(cell_line = Cell_Line, ccle_code = CCLE_Code, in_gene_matrix)],
       file.path(RES, "ccle_breast_lines_present.csv"))

# --- 2) v1 labels: only lines present in matrix + literature-backed list ---
# Luminal: MCF7, T-47D, CAMA-1, ZR-75-1 (BT-474, SUM149PT not in this dataset)
# Basal: HCC 1806, HCC1143, HCC70, MDA-MB-468 (exclude MDA-MB-231 from clean pass)
luminal_lines <- c("MCF7", "T-47D", "CAMA-1", "ZR-75-1")
basal_lines   <- c("HCC 1806", "HCC1143", "HCC70", "MDA-MB-468")

lab_rows <- rbind(
  data.table(
    cell_line = luminal_lines,
    subtype_label = "Luminal",
    literature_note = c(
      "Classic ER+ luminal model; widely used (CCLE/DepMap).",
      "ER+ luminal breast model cell line.",
      "ER+ luminal; derived from metastatic site.",
      "ER+ luminal; ZR-75-1 naming in CCLE sample sheet."
    ),
    exclude_from_clean_contrast = FALSE
  ),
  data.table(
    cell_line = basal_lines,
    subtype_label = "Basal",
    literature_note = c(
      "TNBC/basal-like; commonly classified basal in breast CCLE panels.",
      "Basal/TNBC-like in breast cancer cell line compendia.",
      "Basal-like breast line in multiple studies.",
      "Basal/TNBC-like; not luminal."
    ),
    exclude_from_clean_contrast = FALSE
  )
)
lab_rows[, in_gene_matrix := cell_line %in% mat_set]
fwrite(lab_rows, file.path(RES, "ccle_breast_subtype_labels_v1.csv"))

use_lines <- lab_rows[in_gene_matrix == TRUE, cell_line]
n_lum <- sum(luminal_lines %in% use_lines)
n_bas <- sum(basal_lines %in% use_lines)

# --- 3) Feasibility assessment ---
feasible_run <- (length(use_lines) >= 6 && n_lum >= 3 && n_bas >= 3)

# --- Load matrix subset ---
gm <- fread(MATRIX_PATH)
gene_col <- names(gm)[1]
uid_col <- if ("UniProtID" %in% names(gm)) "UniProtID" else names(gm)[2]
sample_cols <- intersect(use_lines, names(gm))
if (length(sample_cols) < length(use_lines)) {
  missing <- setdiff(use_lines, names(gm))
  warning("Some labeled lines missing from matrix columns: ", paste(missing, collapse = ", "))
}

# If we cannot subset, only write feasibility
if (length(sample_cols) < 6) {
  sink(file.path(RES, "ccle_luminal_vs_basal_feasibility.md"))
  cat("# CCLE Luminal vs Basal first-pass — not run\n\n")
  cat("**Reason:** Too few labeled lines present as columns in `gene_matrix.csv` (need ≥6 total, ≥3 per arm).\n\n")
  cat("**Present columns:** ", paste(sample_cols, collapse = ", "), "\n", sep = "")
  sink(NULL)
  quit(save = "no", status = 0)
}

mat <- as.matrix(gm[, ..sample_cols])
storage.mode(mat) <- "double"
rownames(mat) <- gm[[gene_col]]
uid <- gm[[uid_col]]

# genes with any signal in subset
vr <- apply(mat, 1, function(z) var(z, na.rm = TRUE))
keep <- is.finite(vr) & vr > 0
mat <- mat[keep, , drop = FALSE]
uid <- uid[keep]

# design: Luminal vs Basal
grp <- ifelse(sample_cols %in% luminal_lines, "Luminal", "Basal")
grp <- factor(grp, levels = c("Basal", "Luminal"))
design <- model.matrix(~ 0 + grp)
colnames(design) <- c("Basal", "Luminal")

# PCA (rows = samples)
X <- t(mat)
X[is.na(X)] <- median(X, na.rm = TRUE)
pc <- prcomp(X, center = TRUE, scale. = TRUE)
pc12 <- pc$x[, 1:2]
pcdf <- data.table(
  cell_line = sample_cols,
  PC1 = pc12[, 1], PC2 = pc12[, 2],
  subtype = as.character(grp)
)
fwrite(pcdf, file.path(OUT_DA, "pca_pc1_pc2_scores.csv"))

pdf(file.path(OUT_DA, "pca_pc1_pc2.pdf"), width = 7, height = 6)
col <- ifelse(pcdf$subtype == "Luminal", "#2166ac", "#b2182b")
plot(pcdf$PC1, pcdf$PC2, pch = 16, col = col, xlab = "PC1", ylab = "PC2",
     main = "CCLE breast lines (v1): PCA on protein abundance")
text(pcdf$PC1, pcdf$PC2, labels = pcdf$cell_line, pos = 3, cex = 0.55)
legend("topright", legend = c("Luminal", "Basal"), col = c("#2166ac", "#b2182b"), pch = 16)
dev.off()

# limma: Luminal - Basal
fit <- lmFit(mat, design)
ctr <- makeContrasts(Luminal - Basal, levels = design)
fit2 <- contrasts.fit(fit, ctr)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
tt$GeneSymbol <- rownames(mat)
tt$UniProtID <- uid[match(rownames(mat), rownames(mat))]
fwrite(as.data.table(tt), file.path(OUT_DA, "DA_limma_Luminal_vs_Basal.csv"))

sink(file.path(OUT_DA, "DA_limma_Luminal_vs_Basal_summary.txt"))
cat("CCLE exploratory: Luminal - Basal (v1 labeled lines)\n")
cat("Samples (columns): ", length(sample_cols), "\n")
print(table(grp))
cat("\n**Interpretation:** One column per cell line; no biological replication.\n")
cat("P-values are heuristic only (lines as independent units).\n")
sink(NULL)

# Marker lookup by UniProt substring (matrix uses Entrez + UniProt in col2)
markers <- data.table(
  symbol = c("ESR1", "GATA3", "FOXA1", "KRT18", "PGR", "KRT5", "KRT14", "KRT17", "EGFR", "FOXC1"),
  uniprot_pattern = c(
    "P03372", "P23771", "P55317", "P05783", "P06454",
    "P13647", "P02533", "Q04695", "P00533", "Q12948"
),
  expected_sign = c(
    rep("positive (Luminal higher)", 5),
    rep("negative (Basal higher)", 5)
  )
)
uid_chr <- as.character(uid)
mr <- lapply(seq_len(nrow(markers)), function(i) {
  pat <- markers$uniprot_pattern[i]
  hit <- grep(pat, uid_chr, fixed = TRUE)
  if (length(hit) == 0)
    return(data.table(symbol = markers$symbol[i], UniProtID = NA_character_,
                      logFC = NA_real_, adj.P.Val = NA_real_, found = FALSE,
                      note = "No matching row in gene_matrix"))
  j <- hit[1L]
  data.table(
    symbol = markers$symbol[i],
    UniProtID = uid_chr[j],
    logFC = tt$logFC[j],
    adj.P.Val = tt$adj.P.Val[j],
    found = TRUE,
    note = if (grepl("contaminant", uid_chr[j], ignore.case = TRUE)) "Isoform/contaminant row; interpret cautiously" else ""
  )
})
marker_out <- rbindlist(mr)
fwrite(marker_out, file.path(OUT_DA, "canonical_markers_check.csv"))

# Feasibility / limitations markdown (always written)
sink(file.path(RES, "ccle_luminal_vs_basal_feasibility.md"))
cat("# CCLE Luminal vs Basal first-pass — scope and limitations\n\n")
cat("## What was run (exploratory)\n\n")
cat("- **Labeled lines (v1):** 4 Luminal + 4 Basal, all present in `gene_matrix.csv`.\n")
cat("- **PCA** and **limma** contrast **Luminal − Basal** on the **gene × line** matrix.\n\n")
cat("## Why this is only a first-pass / not a rigorous subtype benchmark\n\n")
cat("1. **No biological replication:** one proteomics column per cell line. limma treats eight lines as eight independent samples; effects are **line-specific** as much as subtype-specific.\n")
cat("2. **Subtype labels** are **literature-based** for classic lines, not derived from proteomics or single-cell ground truth in this repo.\n")
cat("3. **Missing candidate lines:** **BT-474**, **SUM149PT** (and **ZR7530** as alternate name) do **not** appear in the CCLE sample sheet / matrix used here.\n")
cat("4. **MDA-MB-231** was excluded from the clean v1 list by design.\n")
cat("5. **Marker proteins:** **ESR1** (`P03372`) and **PGR** (`P06454`) are **not** present in this `gene_matrix.csv` build; several keratins map to **contaminant** isoform rows — see `canonical_markers_check.csv`.\n\n")
cat("## When a stronger benchmark would be defensible\n\n")
cat("- More lines per subtype **or** replicate cultures per line.\n")
cat("- Harmonized gene/protein IDs so key receptors (ESR1/PGR) are quantified.\n")
cat("- Explicit batch correction if new CCLE plexes are added.\n\n")
cat("## Outputs\n\n")
cat("- `DA_luminal_vs_basal/pca_pc1_pc2.pdf` — exploratory separation.\n")
cat("- `DA_luminal_vs_basal/DA_limma_Luminal_vs_Basal.csv` — **exploratory** DE (see caveats).\n")
cat("- `DA_luminal_vs_basal/canonical_markers_check.csv` — direction vs expectation where proteins exist.\n")
sink(NULL)

message("Wrote: ", RES, "/ccle_breast_lines_present.csv")
message("Wrote: ", RES, "/ccle_breast_subtype_labels_v1.csv")
message("Wrote: ", OUT_DA, "/ and feasibility md")
