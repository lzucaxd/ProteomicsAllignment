#!/usr/bin/env Rscript
# =============================================================================
# Install R dependencies for ProteomicsAllignment (reproducible setup)
# Run from repo root: Rscript install_r_packages.R
# Requires R >= 4.0.
# =============================================================================

repos_cran <- "https://cloud.r-project.org"

message("Checking / installing BiocManager ...")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = repos_cran)

message("Installing CRAN packages: data.table, ggplot2, tidyr ...")
install.packages(c("data.table", "ggplot2", "tidyr"), repos = repos_cran)

message("Installing Bioconductor packages: MSstatsTMT, org.Hs.eg.db ...")
BiocManager::install(c("MSstatsTMT", "org.Hs.eg.db"), update = FALSE, ask = FALSE)

message("Optional (mouse only): org.Mm.eg.db. Uncomment next line if you need mouse gene mapping.")
# BiocManager::install("org.Mm.eg.db", update = FALSE, ask = FALSE)

message("Done. R dependencies are ready.")
