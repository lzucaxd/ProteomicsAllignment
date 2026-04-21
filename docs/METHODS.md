# Methods (harmonization representations)

This document describes **each representation** used in the benchmark at a level suitable for **Methods** text in a paper or supplement. Statistical evaluation (limma, metrics, nulls, ceilings) is **shared** across methods; see **`scripts/benchmark/README.md`** and **`docs/END_TO_END_TECHNICAL_REPORT.md`**.

---

## Raw (baseline)

**Definition.** No cross-domain harmonization beyond what is already in **MSstatsTMT** study-level processing. Union construction (`scripts/run_preprocessing.py`) concatenates CPTAC and CCLE samples on the **intersection gene set** with prevalence / variance filters.

**Role.** Baseline for both **geometry** and **FC agreement**. CPTAC native limma logFCs for a task are unchanged when adding CCLE samples to the design matrix in downstream **representation-level** DA, because CCLE is a separate “domain” block in the benchmark’s per-domain fits.

---

## Bridge shift

**Goal.** Remove **much of the systematic offset** between CPTAC and CCLE reporter-ion scales using the **TMT bridge / reference (Norm) channel**, then leave biological contrasts **as additive as possible** in log space.

**Where the offset comes from (critical).** **Not** from “domain medians” of sample `gene_matrix` abundances. The pipeline reads **`msstats_input.tsv`** and filters rows with **`Condition == "Norm"`** (see `scripts/benchmark/extract_bridge_summaries.R`). For each protein and each plex, it summarizes bridge intensities, then forms **cross-plex robust location** summaries per protein per domain. The CCLE matrix is shifted by a **per-gene** offset targeting CPTAC’s bridge baseline.

**Application.** **CCLE** columns are shifted on the log₂ scale; **CPTAC** is left unchanged in the stored representation used here.

**Fold changes.** A **constant per gene** added to all CCLE samples cancels in **within-CCLE** contrasts, so **CCLE luminal–basal logFCs are unchanged**. CPTAC logFCs are unchanged because CPTAC was not shifted. **Cross-domain FC correlation** can still differ from raw when genes are filtered or when agreement is computed on matched gene sets with numerical edge cases, but the **intended** invariant is **no per-gene distortion inside CCLE**.

**Robustness.** Genes with **too few** informative bridge plexes may be passed through with weaker or no correction (implementation detail in `bridge_aware_correction.R`).

---

## Bridge shift + scale

**Definition.** Same **bridge-derived location** alignment as **Bridge shift**, plus **MAD-based** rescaling of CCLE per gene to better match CPTAC bridge **dispersion** when S/N compression differs.

**Fold changes.** Scaling is **per gene** and **linear in log space** for CCLE; it **does** change **magnitude** (not only location) of CCLE contrasts unless scale factors are near 1. Treat as **approximately** preserving small effects when scale ≈ 1; always check **`comparison_summary.csv`** rather than assuming exact FC preservation.

---

## Celligner (adapted to proteomics)

**Origin.** Celligner-style pipelines were popularized for **RNA-seq** integration (contrastive PCA + MNN-style alignment). This repository vendors a **Python** implementation under `models/celligner-master/` and applies it to **log₂ gene matrices**.

**Stages (conceptual).**

1. **Contrastive PCA (cPCA):** emphasize variation correlated with a “background” (e.g. broad CPTAC tumor heterogeneity) while retaining CCLE-relevant directions.
2. **Mutual nearest neighbors (MNN)** in a moderate-dimensional PC subspace: learn local correction vectors from anchor pairs.
3. **Back-transform** to gene space via loadings: **all genes** are moved **jointly** — there is **no** per-gene guarantee that a native CPTAC marker such as **FOXA1** retains its luminal–basal logFC.

**Preprocessing interaction.** Celligner is run on matrices after **union** filtering; **high tissue specificity** in CPTAC can reduce prevalence for some markers in pan-tissue unions — interpret marker panels with care (see diagnostics CSVs).

**When it shines / fails.** Often **excellent** domain mixing on PCs with **large** changes to **within-domain** DA relative to native CPTAC inference — exactly the **disconnect** the benchmark quantifies.

---

## Evaluation metrics (pointer)

| Concept | Where computed |
|---------|----------------|
| FC correlation, same-direction | `compute_cross_domain_metrics.R` → `comparison_summary.csv` |
| Permutation null | `run_permutation_nulls.R` → `calibration/null_distribution.csv` |
| Concordance ceiling | `run_concordance_ceilings.R` → `calibration/concordance_ceiling_*.csv` |
| Stratified FC | `compute_stratified_fc.R` → `reports/benchmark_master/diagnostics/` |
| Disconnect | `compute_disconnect_scores.R` → `disconnect_scores.csv` |

---

## References (informal)

- Warren et al., *Nature Communications* (2019) — cPCA.  
- Haghverdi et al. / Batchelor et al. — MNN batch integration (as implemented in `mnnpy`).  
- MSstatsTMT — summarization and design for TMT proteomics.
