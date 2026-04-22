#!/usr/bin/env Rscript
# =============================================================================
# Tumor vs NAT differential abundance — CPTAC breast (PDC000120) using MSstatsTMT
#
# Input:
#   - results/PDC000120/protein_summary.tsv (from MSstatsTMT proteinSummarization)
#   - biospecimen/PDC_study_biospecimen_03162026_190026.csv (Tumor vs NAT labels)
#
# Output:
#   - results/PDC000120/DA_MSstatsTMT_tumor_vs_NAT.csv
#   - results/PDC000120/volcano_MSstatsTMT.pdf
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  if (!requireNamespace("MSstatsTMT", quietly = TRUE))
    BiocManager::install("MSstatsTMT", update = FALSE, ask = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel", repos = "https://cloud.r-project.org")
})

library(data.table)
library(MSstatsTMT)
library(ggplot2)
library(ggrepel)

# Be conservative with threading (helps avoid OpenMP shared-memory issues on macOS sandboxes)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

if (length(sys.frames()) > 0 && exists("ofile", sys.frame(1))) {
  DATA_DIR <- dirname(dirname(sys.frame(1)$ofile))
} else {
  DATA_DIR <- getwd()
  if (!file.exists(file.path(DATA_DIR, "results", "PDC000120", "protein_summary.tsv")))
    DATA_DIR <- file.path(getwd(), "data")
}
setwd(DATA_DIR)

RESULTS_DIR <- file.path(DATA_DIR, "results", "PDC000120")
PROT_SUM_PATH <- file.path(RESULTS_DIR, "protein_summary.tsv")
BIOSPEC_PATH <- file.path(DATA_DIR, "biospecimen", "PDC_study_biospecimen_03162026_190026.csv")
OUT_CSV <- file.path(RESULTS_DIR, "DA_MSstatsTMT_tumor_vs_NAT.csv")
OUT_PDF <- file.path(RESULTS_DIR, "volcano_MSstatsTMT.pdf")

stopifnot(file.exists(PROT_SUM_PATH))
stopifnot(file.exists(BIOSPEC_PATH))

message("Loading protein summaries")
prot <- fread(PROT_SUM_PATH)
setnames(prot, trimws(gsub("^\uFEFF", "", names(prot))))

req_cols <- c("Protein", "Abundance", "BioReplicate", "Condition", "Mixture", "TechRepMixture", "Run", "Channel")
missing_cols <- setdiff(req_cols, names(prot))
if (length(missing_cols) > 0) stop("protein_summary.tsv missing columns: ", paste(missing_cols, collapse = ", "))

message("Loading biospecimen labels")
bio <- fread(BIOSPEC_PATH)
setnames(bio, trimws(gsub("^\uFEFF", "", names(bio))))
aliquot_col <- names(bio)[grepl("Aliquot.*Submitter|Submitter.*ID", names(bio), ignore.case = TRUE)][1]
type_col <- names(bio)[grepl("Sample.Type", names(bio), ignore.case = TRUE)][1]
if (is.na(aliquot_col)) aliquot_col <- "Aliquot Submitter ID"
if (is.na(type_col)) type_col <- "Sample Type"
aliquot_to_type <- setNames(trimws(bio[[type_col]]), trimws(bio[[aliquot_col]]))

to_group <- function(sample_type) {
  ifelse(sample_type == "Primary Tumor", "Tumor",
         ifelse(sample_type == "Solid Tissue Normal", "NAT", NA_character_))
}

prot[, BioReplicate := trimws(as.character(BioReplicate))]
prot[, Condition := trimws(as.character(Condition))]

# Remap non-Norm channels to Tumor/NAT using biospecimen; keep Norm as "Norm"
prot[, Condition2 := fifelse(tolower(Condition) == "norm", "Norm", to_group(aliquot_to_type[BioReplicate]))]

keep <- !is.na(prot$Condition2)
message("Keeping rows with Condition in {Tumor, NAT, Norm}: ", sum(keep), " / ", nrow(prot))
prot <- prot[keep]
prot[, Condition := Condition2][, Condition2 := NULL]

message("Counts: Tumor ", sum(prot$Condition == "Tumor"), ", NAT ", sum(prot$Condition == "NAT"), ", Norm ", sum(prot$Condition == "Norm"))

# Safety: drop blank/NA Abundance
prot <- prot[!is.na(Abundance) & nzchar(as.character(Abundance))]

message("Running MSstatsTMT groupComparisonTMT (Tumor vs NAT)")
contrast <- matrix(c(1, -1), nrow = 1, dimnames = list("Tumor-NAT", c("Tumor", "NAT")))

gc <- groupComparisonTMT(
  # groupComparisonTMT expects the list output of proteinSummarization()
  # and accesses data$ProteinLevelData internally.
  data = list(ProteinLevelData = prot),
  contrast.matrix = contrast,
  moderated = TRUE,
  adj.method = "BH",
  remove_norm_channel = TRUE,
  remove_empty_channel = TRUE,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = file.path(RESULTS_DIR, "MSstatsTMT_groupComparison.log")
)

res <- as.data.frame(gc$ComparisonResult)

# Normalize column names for downstream compatibility
if (!"log2FC" %in% names(res) && "logFC" %in% names(res)) res$log2FC <- res$logFC
if (!"adj.pvalue" %in% names(res) && "adj.pvalue" %in% names(res)) res$adj.pvalue <- res$adj.pvalue
if (!"adj.pvalue" %in% names(res) && "Adjusted.Pvalue" %in% names(res)) res$adj.pvalue <- res$Adjusted.Pvalue

fwrite(res, OUT_CSV)
message("Wrote: ", OUT_CSV)

message("Plotting volcano")
if (!("Protein" %in% names(res)) && ("ProteinName" %in% names(res))) res$Protein <- res$ProteinName
if (!("Protein" %in% names(res))) res$Protein <- res$ProteinName
if (!("adj.pvalue" %in% names(res)) && ("adj.P.Val" %in% names(res))) res$adj.pvalue <- res$adj.P.Val

plot_dt <- as.data.table(res)
plot_dt[, neglog10FDR := -log10(pmax(adj.pvalue, 1e-300))]
plot_dt[, sig := abs(log2FC) > 1 & adj.pvalue < 0.05]

top_n <- 20L
plot_dt[, rank := rank(adj.pvalue, ties.method = "first")]
plot_dt[, label := ifelse(sig & rank <= top_n, as.character(Protein), NA_character_)]

markers <- c("ESR1", "GATA3", "KRT18", "PTPRC", "COL1A1", "VIM")
plot_dt[Protein %in% markers, label := as.character(Protein)]

p <- ggplot(plot_dt, aes(x = log2FC, y = neglog10FDR)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.1) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55") +
  ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3, box.padding = 0.4, min.segment.length = 0) +
  labs(
    title = "MSstatsTMT Tumor vs NAT (PDC000120)",
    x = "log2FC (Tumor - NAT)",
    y = "-log10(FDR)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(OUT_PDF, p, width = 8.5, height = 6.5)
message("Wrote: ", OUT_PDF)

sig_hits <- plot_dt[sig == TRUE]
message("Significant hits (FDR<0.05, |log2FC|>1): ", nrow(sig_hits))

