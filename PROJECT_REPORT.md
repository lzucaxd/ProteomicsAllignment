# ProteomicsAllignment — Project report

This document summarizes **what this repository does**: core PDC → MSstatsTMT pipelines, CPTAC BRCA (PDC000120) differential abundance and benchmarking, CCLE utilities, reports, and housekeeping.

**Harmonization benchmark v2 (slides, outputs, paths):** see **[docs/BENCHMARK_V2_AND_PRESENTATION.md](docs/BENCHMARK_V2_AND_PRESENTATION.md)** and `scripts/benchmark/run_overnight_v2.sh`.

---

## 1. Purpose

**Primary goal:** Turn **Proteomic Data Commons (PDC)** TMT PSM-level data into **sample × gene abundance matrices** using **MSstatsTMT**, with reproducible scripts and optional downstream **differential abundance (DA)** for alignment / benchmarking narratives.

**Secondary goal:** Support **cross-study** and **methods** work: standardized matrices, consensus DA across limma / MSstats / MSstatsTMT, subtype-aware CPTAC BRCA analyses, and written reports for presentations.

---

## 2. Repository layout (high level)

| Location | Role |
|----------|------|
| Repo root | `README.md`, `pyproject.toml`, `requirements.txt`, `install_r_packages.R`, this report |
| `docs/` | Design notes; **BENCHMARK_V2_AND_PRESENTATION.md** = slide + output hub |
| `src/harmonize/` | Python package (preprocessing, benchmark utilities) |
| `configs/` | Preprocessing / task YAML |
| `scripts/benchmark/` | **Overnight v2** orchestrator + limma, calibration, figures |
| `models/celligner-master/` | Vendored Celligner + mnnpy |
| `data/` | **Main working directory** for PDC pipelines (run scripts from here) |
| `data/manifests/` | PDC file manifest CSVs per study |
| `data/pdc_psm/` | Downloaded PSM files (large, git-ignored) |
| `data/results/{study_id}/` | Per-study outputs: gene matrix, protein summary, annotation, DA (git-ignored) |
| `data/scripts/` | R/Python DA, CPTAC subtype mapping, CCLE DA, Venn script, utilities |
| `data/biospecimen/` | CPTAC clinical / biospecimen inputs |
| `data/ccle_peptide/` | CCLE peptide → MSstats input utilities |
| `reports/` | Narrative reports, figures, R Markdown, pathway notes, **presentation_subtype_benchmark/** |
| `scripts/` | **`benchmark/`** (overnight v2), **`methods/`** (Celligner driver); see [scripts/README.md](scripts/README.md). Most DA / CPTAC one-offs remain under **`data/scripts/`** |

**Git-ignored by default:** `data/pdc_psm/`, `data/results/`, `*.log`, `.DS_Store` (see `.gitignore`).

---

## 3. Core pipeline: PDC PSM → MSstatsTMT → gene matrix

**Documentation:** [`data/PIPELINE_README.md`](data/PIPELINE_README.md)  
**Orchestration:** `data/run_pipeline_per_manifest.sh`  
**Entry point (typical):**

```bash
cd data
./run_pipeline_per_manifest.sh
```

**Flow:**

1. **Download** PSM files from PDC manifests into `pdc_psm/{study_id}/`.
2. **Parse** PSMs → long format; build **MSstatsTMT** input (Run, Channel, Condition, BioReplicate, Mixture, intensities).
3. **MSstatsTMT:** reference-channel normalization, peptide → **protein** summarization.
4. **Map** proteins → **gene symbols** (`org.Hs.eg.db` / optional mouse DB).
5. **Output** per study under `data/results/{study_id}/`: e.g. `gene_matrix.csv`, `protein_summary.tsv`, `annotation_filled_corrected.csv`, QC plots (optional).

**Important:** The core pipeline **stops at the gene matrix**; it does **not** do ML preprocessing (imputation, batch correction for ML, feature selection)—by design.

**Related scripts:**

| Script | Purpose |
|--------|---------|
| `data/pdc_manifest_downloader.py` | Download PSMs from manifests |
| `data/pdc_psm_to_protein_matrix.R` / `pdc_psm_to_msstatsTMT_protein_matrix.R` | PSM → protein / MSstatsTMT path |
| `data/run_pdc_psm_matrix.sh` / `data/run_msstats_tmt_gene_matrix.sh` | Shell wrappers |
| `data/msstatsTMT_qc_plots.R` | QC plotting |
| `data/check_manifests.py`, `data/check_studies_sample_file.py`, `data/test_pdc_download.py` | Validation helpers |
| `data/cleanup_study_disk.sh` | Optional deletion of huge intermediates (`parsed_psm_long.tsv`, `msstats_input.tsv`, optionally raw PSM) after `gene_matrix.csv` exists |

---

## 4. CPTAC BRCA study PDC000120 — downstream analyses

Study outputs live under **`data/results/PDC000120/`** (when generated). The following scripts assume that directory contains at least `gene_matrix.csv`, `annotation_filled_corrected.csv`, and often `protein_summary.tsv`.

### 4.1 Subtype and sample annotation (Python)

**Script:** `data/scripts/build_PDC000120_subtype_mapping.py`

**Inputs:** CPTAC clinical/biospecimen, bridge/sample mapping, `annotation_filled_corrected.csv`, `gene_matrix.csv`.

**Outputs (examples):** `DA_sample_annotation.csv`, `DA_subtype_tumor_only.csv`, `DA_subtype_counts.csv`, `gene_matrix_subtype_mapping.csv`, `subtype_DA_recommendations.txt`, plus optional intermediate mapping tables when the full script is run.

**Goal:** One analysis-ready row per matrix column: **Tumor vs NAT** from biospecimen, **PAM50** for tumor-only subtype DA.

### 4.2 Tumor vs NAT differential abundance

| Script | What it does | Main outputs |
|--------|----------------|--------------|
| `DA_tumor_vs_NAT_CPTAC_breast.R` | **MSstats** + **limma** on gene matrix; merged comparison | `DA_MSstats_tumor_vs_NAT.csv`, `DA_limma_tumor_vs_NAT.csv`, volcanos → usually under `DA_tumor_vs_NAT/` |
| `DA_tumor_vs_NAT_CPTAC_breast_MSstatsTMT.R` | **MSstatsTMT** on protein-level summary | `DA_MSstatsTMT_tumor_vs_NAT.csv`, volcano |

**Folder:** `data/results/PDC000120/DA_tumor_vs_NAT/` — holds consensus inputs, volcanos (limma, MSstats, MSstatsTMT), scatter plots, upset plot, p-value histograms, enrichment CSVs/PDFs when run.

### 4.3 Consensus DA + pathway enrichment (R)

**Script:** `data/scripts/run_consensus_DA_analysis_PDC000120.R`

**Inputs:** The three Tumor vs NAT result CSVs above.

**Outputs:** `consensus/` — `DA_consensus_table.csv`, high-confidence tumor/NAT lists, method correlations/overlap, **Hallmark / GO BP / Reactome** enrichment tables and plots, summary text (`DA_consensus_summary.txt`).

**Note:** Pathway ORA uses **clusterProfiler** + **msigdbr**; universe choices are documented in [`reports/pathway_enrichment_universe_note.md`](reports/pathway_enrichment_universe_note.md).

### 4.4 Subtype DA (Basal vs Luminal, etc.)

| Script | Method | Notes |
|--------|--------|--------|
| `DA_subtype_MSstats_PDC000120.R` | Gene matrix → long format; **limma** fallback (mixture as batch) | Volcano, summary, `DA_MSstats_<contrast>.csv` |
| `DA_subtype_MSstatsTMT_PDC000120.R` | **MSstatsTMT** `groupComparisonTMT` on `protein_summary.tsv` | Coverage filters; **RefSeq → gene symbol**; `*_gene_symbols.csv`, `*_marker_sanity.csv`, volcano |
| `DA_subtype_CPTAC_breast.R` | Older/generic subtype limma helper | Referenced in recommendations |

**Typical contrast:** **Basal vs Luminal** (LumA + LumB), tumor-only, PAM50-filtered samples.

### 4.5 Other

- `data/scripts/run_consensus_DA_analysis_PDC000120.R` — also generates comparison plots (scatter, upset) referenced in benchmarking.
- **`Rscript --vanilla`** is recommended if `renv/activate.R` is missing locally.

---

## 5. CCLE peptide pipeline (parallel to CPTAC, no PSM download)

**Goal:** Run the **same** MSstatsTMT → protein → gene-matrix stack used for CPTAC, but starting from **CCLE peptide-level TMT tables** already in hand (no PDC manifests, no `.psm` parsing).

**Full detail:** [`data/ccle_peptide/README_CCLE.md`](data/ccle_peptide/README_CCLE.md)

### 5.1 How it differs from the PDC path

| Stage | CPTAC (PDC) | CCLE |
|-------|----------------|------|
| Input | PSM files + `.sample.txt` / manifest | Peptide TSV (`rq_*_sn` reporter columns) + sample spreadsheet |
| Parsing | PSM → long format inside R | **Python** builds MSstatsTMT long input |
| R entry | `pdc_psm_to_msstatsTMT_protein_matrix.R --psm_dir …` | **Same R script** with **`--msstats_input_dir`** (skips PSM parse; loads pre-built `msstats_input.tsv`) |
| Output | `gene_matrix.csv`, `protein_summary.tsv`, … | Same artifacts under e.g. `data/results/CCLE/` |

So: **one R driver** (`data/pdc_psm_to_msstatsTMT_protein_matrix.R`); CCLE only supplies a folder that already contains `msstats_input.tsv` + `annotation_filled.csv`.

### 5.2 Data flow (CCLE)

```
CCLE sample workbook (Sheet2)  ──►  export_sample_sheet2_csv.py  ──►  sample_info_ccle.csv
        +
Peptide TSV (ProteinId, PeptideSequence, rq_*_sn, RunLoadPath, …)
        │
        ▼
ccle_to_msstats_input.py  ──►  results/CCLE/msstats_input.tsv
                            ──►  results/CCLE/annotation_filled.csv
        │
        ▼
Rscript pdc_psm_to_msstatsTMT_protein_matrix.R --msstats_input_dir results/CCLE --outdir results/CCLE
        │
        ▼
protein_summary.tsv, gene_matrix.csv, plots/, qc_summary.txt, …
```

### 5.3 Inputs (typical)

- **`Table_S1_Sample_Information (1).xlsx`** — CCLE sample metadata; **Sheet2** is the channel / plex design (same *role* as CPTAC sample design).
- **`ccle_protein_quant_with_peptides_*.tsv`** — Peptide-level table with **`rq_*_sn`** columns (reporter-ion signal-to-noise); channel columns are **auto-detected** (supports TMT6/10/11/16/TMTpro/18-style labels).

### 5.4 Sample sheet → MSstatsTMT columns

The converter maps CCLE columns to the same annotation schema as CPTAC:

| CCLE / Sheet2 | MSstatsTMT field |
|---------------|------------------|
| Cell Line | **BioReplicate** |
| Protein 10-Plex ID | **Mixture** (plex id) |
| Protein TMT Label | **Channel** (126, 127N, 127C, … canonicalized) |
| Notes containing “Bridge line” | **Condition** = `Norm` (reference); else `Sample` |
| Basename of `RunLoadPath` (from peptide table) | **Run** |
| (fixed) | **TechRepMixture** = 1 |

**Rules:** Each mixture must have **exactly one** bridge (**Norm**) per plex; **(Mixture, Channel)** pairs must be unique; all mixtures should share the **same channel set** (validation enforces this unless overridden by flags—see README).

**Run ↔ Mixture:** Plex IDs come from the sample file, not run order. Runs map to mixtures via run name patterns (e.g. `Prot_01` → mixture) or by matching the set of channels with data to a plex in the sample file (documented in README_CCLE).

### 5.5 Scripts (CCLE-specific)

| Script | Role |
|--------|------|
| `data/ccle_peptide/export_sample_sheet2_csv.py` | Export Sheet2 → `sample_info_ccle.csv` (stdlib). |
| `data/ccle_peptide/ccle_to_msstats_input.py` | TSV + sample CSV → `msstats_input.tsv` + `annotation_filled.csv`; validates channel/Norm consistency before write. |
| `data/pdc_psm_to_msstatsTMT_protein_matrix.R` | **`--msstats_input_dir <dir>`** — load pre-built input, run **proteinSummarization** → **gene mapping** → outputs. |
| `data/scripts/ccle_qc_plots_from_protein_tsv.R` | Light QC **PNG**s from existing `protein_summary.tsv` (no re-summarization). |
| `data/scripts/ccle_regenerate_qc_plots.R` | Re-run summarization (slow) to feed MSstatsTMT’s **QCPlot** / standard PDF QC if needed. |
| `data/scripts/check_ccle_matrix.py` | Sanity-check **`gene_matrix.csv`**: shape, missingness, PCA plot (`pca_plot.png`), annotation checks (needs pandas/sklearn). |

### 5.6 Commands (minimal)

From **`data/`** (paths adjusted to your TSV name):

```bash
# 1) Sample sheet once
python3 ccle_peptide/export_sample_sheet2_csv.py

# 2) Build MSstatsTMT input
python3 ccle_peptide/ccle_to_msstats_input.py \
  --tsv ccle_peptide/ccle_protein_quant_with_peptides_14745.tsv \
  --sample_csv ccle_peptide/sample_info_ccle.csv \
  --outdir results/CCLE

# 3) Same R pipeline as CPTAC, CCLE mode
Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R \
  --msstats_input_dir results/CCLE \
  --outdir results/CCLE
```

**CPTAC is unchanged:** studies that use `--psm_dir` and manifests keep the original PDC flow; only directories passed as **`--msstats_input_dir`** use the pre-built table path.

---

## 6. Reports and documentation (for presentations / methods)

| File | Content |
|------|---------|
| [`reports/Proteomics_Alignment_Benchmarking_Results.md`](reports/Proteomics_Alignment_Benchmarking_Results.md) | Interpreted results: Tumor vs NAT concordance, Basal vs Luminal markers, file references |
| [`reports/pathway_enrichment_universe_note.md`](reports/pathway_enrichment_universe_note.md) | Hypergeometric background / proteomics caveats; what the consensus script actually uses |
| [`reports/README.md`](reports/README.md) | How to render R Markdown report, generate figure PDFs, volcano locations |
| [`reports/CPTAC_BRCA_DA_report.Rmd`](reports/CPTAC_BRCA_DA_report.Rmd) | Optional HTML report (needs Pandoc) |
| [`reports/generate_report_plots.R`](reports/generate_report_plots.R) | Standalone PDF figures under `reports/figures/` |
| [`data/results/PDC000120/README.md`](data/results/PDC000120/README.md) | What to keep vs regenerate in that folder |

**Figures:** `reports/figures/` — method overlap, top hits, enrichment, Basal/Luminal tops, subtype markers; copies of Tumor vs NAT volcanos (`volcano_tumor_vs_NAT_*.pdf`).

---

## 7. Utilities

| File | Purpose |
|------|---------|
| `data/scripts/notify_when_DA_finishes.sh` | macOS notification when MSstatsTMT subtype DA script exits |
| `data/cleanup_study_disk.sh` | Free disk space by removing large intermediates per study |
| `install_r_packages.R` | One-shot CRAN/Bioc installs for the core pipeline |

---

## 8. Dependencies (summary)

- **Python:** `requests`, `pandas` (subtype script), etc. — see `requirements.txt`.
- **R:** `MSstatsTMT`, `MSstats`, `limma`, `data.table`, `ggplot2`, `org.Hs.eg.db`; optional `clusterProfiler`, `msigdbr`, `enrichplot`, `ReactomePA` for consensus enrichment.

---

## 9. What this repo does *not* centralize

- Raw PDC credentials / private data (manifests are user-supplied).
- Long-term archival of **large** `results/` or `pdc_psm/` (ignored in git).
- Single “one command” for every DA variant — scripts are **per analysis** and documented above.

---

## 10. Subtype benchmarking (CPTAC vs CCLE) — MSstatsTMT alignment & slides

This section summarizes **work done to compare Basal vs Luminal signal** across **CPTAC tumors (PDC000120)** and **CCLE breast lines**, with **slide-ready** tables, Venns, and methods notes. Full narrative and file map: **[`reports/BENCHMARKING_AND_SLIDES_REPORT.md`](reports/BENCHMARKING_AND_SLIDES_REPORT.md)**.

**Benchmark v1 (frozen 2026-04-06):** **[`reports/benchmark_v1/`](reports/benchmark_v1/)** — spec, gene-level shared table, raw metrics, diagnostics pack, consistency audit, thesis draft section, Sarah prep.

### 10.1 Unified inference (headline comparison)

For a **fair protein-level comparison** to CPTAC, CCLE subtype DA uses the **same R call** to **`MSstatsTMT::groupComparisonTMT`** as `DA_subtype_MSstatsTMT_PDC000120.R`: same `contrast.matrix` (Luminal − Basal), `moderated = TRUE`, `adj.method = "BH"`, channel cleanup flags. Implementation: **`data/scripts/ccle_DA_luminal_basal_v1.R`** on **`data/results/CCLE_corrected/protein_summary.tsv`**.

### 10.2 CPTAC: mixture-balanced **subset** (clinical benchmark)

**Problem:** In multi-plex TMT data, some plexes contain **only one subtype** (or a singleton minority). Subtype differences can then align with **which plex** you are in, weakening interpretation.

**Response:** **`data/scripts/build_cptac_basal_luminal_mixture_subset.py`** builds a **tumor-level** keep/drop list; the headline subtype run uses annotations restricted to **kept** mixtures. Typical headline counts: **49 Luminal / 26 Basal** tumors. MSstatsTMT outputs: **`data/results/PDC000120/DA_subtype_subset_runs/`** (e.g. `DA_MSstatsTMT_Luminal_vs_Basal.csv`, marker sanity, volcano).

**Canonical marker panel (subset):** 11 **protein rows** (two ESR1 accessions): **11/11** expected direction; **10/11** FDR &lt; 0.05 (KRT14 NS); **4** rows with FDR &lt; 0.05 and |log₂FC| &gt; 1 (both ESR1 isoforms + GATA3 + FOXA1).

### 10.3 CCLE: no mixture subset (design limit)

The eight benchmark lines are **one cell line per TMT plex**; **no** plex contains both a Luminal and a Basal line. A CPTAC-style mixture-balanced subset would remove **all** lines—so it is **not applicable**. **`groupComparisonTMT`** still receives **Mixture / Run / Channel** in `ProteinLevelData`; the limitation is **identifiability** (subtype confounded with plex at line level), not “turning off” MSstatsTMT.

**Sensitivity:** leave-one-line-out diagnostics — **`data/scripts/ccle_subtype_sensitivity_and_influence.R`**, summaries under **`data/results/CCLE_corrected/DA_luminal_vs_basal/`**.

### 10.4 Presentation deliverables (paths)

| Location | Contents |
|----------|----------|
| `reports/presentation_subtype_benchmark/` | Deck context, DA summaries, overlap, QC notes — start with **`README.md`**, **`LLM_CONTEXT_SUMMARY.md`** |
| `reports/presentation_subtype_benchmark/03_ccle_subtype/` | Marker slide CSVs, **Venns** (`venn_figures/*.png`), **headline bullets** (`slide_benchmark_headlines_CPTAC_vs_CCLE.md`), **model / subsetting** (`MSstatsTMT_model_and_subsetting_for_slides.md`, `mixture_subsetting_and_model_note.md`) |
| `data/scripts/plot_subtype_benchmark_venns.py` | Regenerate gene-level FDR Venns (needs local `data/results/` CSVs) |

**Gene vs protein counts (CCLE Venns):** **103** protein rows FDR &lt; 0.05 vs **100** unique **`Gene_symbol`** (MAPT isoforms + one row without symbol)—expected when deduplicating to genes.

### 10.5 Note on two CCLE DE conventions

Some **overlap / Jaccard** tables in `reports/presentation_subtype_benchmark/05_overlap_and_summary/` use **limma on the corrected gene matrix** as the “primary” CCLE row for **gene-list** overlap. The **MSstatsTMT protein-level** CCLE run above is the **method-aligned** arm to CPTAC; use the convention that matches your slide caption.

---

## 11. Quick reference: “what to run for what”

| Goal | Command / action |
|------|------------------|
| Full PDC → gene matrix | `cd data && ./run_pipeline_per_manifest.sh` |
| **CCLE** peptide → MSstats input | `python3 ccle_peptide/ccle_to_msstats_input.py --tsv … --sample_csv … --outdir results/CCLE` (after `export_sample_sheet2_csv.py`; see §5) |
| **CCLE** → protein / gene matrix | `cd data && Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R --msstats_input_dir results/CCLE --outdir results/CCLE` |
| **CCLE** matrix QC (optional) | `python data/scripts/check_ccle_matrix.py --matrix data/results/CCLE/gene_matrix.csv` |
| Subtype mapping (PDC000120) | `python data/scripts/build_PDC000120_subtype_mapping.py` (from `data/` with paths as in script) |
| Tumor vs NAT (matrix) | `cd data && Rscript scripts/DA_tumor_vs_NAT_CPTAC_breast.R` |
| Tumor vs NAT (TMT protein) | `cd data && Rscript scripts/DA_tumor_vs_NAT_CPTAC_breast_MSstatsTMT.R` |
| Consensus + enrichment | `cd data && Rscript scripts/run_consensus_DA_analysis_PDC000120.R` |
| Subtype DA (limma path) | `Rscript --vanilla data/scripts/DA_subtype_MSstats_PDC000120.R` |
| Subtype DA (MSstatsTMT path) | `Rscript --vanilla data/scripts/DA_subtype_MSstatsTMT_PDC000120.R` (set `PDC_SUBTYPE_ANNOT` / `PDC_MSSTATSTMT_OUT_DIR` for mixture subset run — see script env) |
| **CCLE** Luminal vs Basal (**MSstatsTMT**, aligned to CPTAC) | `Rscript --vanilla data/scripts/ccle_DA_luminal_basal_v1.R` |
| **CPTAC** mixture subset build | `python3 data/scripts/build_cptac_basal_luminal_mixture_subset.py` (from `data/` as documented in script) |
| Subtype benchmark Venns | `.venv-venn/bin/python data/scripts/plot_subtype_benchmark_venns.py` (optional venv with matplotlib + matplotlib-venn) |
| Report figures (no Pandoc) | `Rscript --vanilla reports/generate_report_plots.R` |

---

*Last updated to include CPTAC–CCLE subtype benchmarking and presentation deliverables; regenerate large outputs locally after code changes.*
