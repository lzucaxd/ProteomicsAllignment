# Inference baselines: MSstatsTMT (native TMT) vs limma (gene matrix)

This repo uses **two different statistical “layers”**. Keeping them straight is essential for handoff and for Methods text.

---

## 1. Protein-level, TMT-design-aware (**MSstatsTMT**) — CPTAC native

**Where:** Building **`gene_matrix.csv`** from PSMs — **`data/pdc_psm_to_msstatsTMT_protein_matrix.R`** (and MSstatsTMT steps inside it). Uses run/channel/mixture/bridge structure from **`.sample.txt`**.

**Role:** Correct **TMT** summarization and normalization for CPTAC before genes are collapsed to symbols.

**Further “native” DA on CPTAC** (optional / reporting, not the same as the harmonization benchmark’s 16 limma runs):

- R driver: **`scripts/benchmark/native_domain_da.R`**
- Outputs: **`reports/benchmark_master/native_domain_da/`** (CPTAC often MSstatsTMT-style; CCLE may use limma on the gene matrix — see script headers).

**Deep dive:** [`native_domain_inference_overview.md`](native_domain_inference_overview.md).

---

## 2. Gene-level **limma** on matrices — representation-level benchmark

**Where:** After harmonization methods produce **`transformed_{task}.csv`**, the overnight pipeline runs **`scripts/benchmark/run_all_limma_da.R`** (and related helpers) per **domain** (CPTAC vs CCLE) and **method**.

**Role:** **Fair comparison** across methods using the **same** contrast and **gene-level** logFCs on the **shared** representation (raw, bridge, Celligner, …).

**Config hooks:** `configs/tasks/*.yaml` list `native_domain_inference` vs `representation_level_inference` keys consumed by Python task objects; the **overnight** path is scripted in R/shell as documented in **`scripts/benchmark/README.md`**.

**Deep dive (representation / agreement):** [`representation_level_inference_overview.md`](representation_level_inference_overview.md).

---

## Side-by-side

| Layer | Typical tool | Input feature | Question answered |
|-------|--------------|---------------|-------------------|
| **CPTAC PSM → matrix** | MSstatsTMT + mapping | PSM / peptide → protein → gene | Valid TMT summarization |
| **Native-domain DA** | MSstatsTMT (CPTAC) / limma (CCLE) as implemented | `msstats_input` / gene matrix | “What each domain says alone” |
| **Representation benchmark** | limma on aligned matrices | Harmonized gene × sample | “Do harmonizers agree across domains on the **same** genes?” |

---

## Exploratory subtype scripts (`data/scripts/`)

Files like **`DA_subtype_MSstatsTMT_PDC000120.R`** vs **`ccle_DA_luminal_basal_limma_gene_matrix.R`** implement **separate** exploratory contrasts for slides — **not** the same object as the harmonization benchmark’s `run_all_limma_da.R` block. See **`data/scripts/README.md`**.

---

## What to cite in a paper

- **Matrix construction (CPTAC TMT):** MSstatsTMT workflow in **`data/PIPELINE_README.md`** / **`data/pdc_psm_to_msstatsTMT_protein_matrix.R`**.
- **Cross-method benchmark:** limma protocol + metrics in **`docs/END_TO_END_TECHNICAL_REPORT.md`** and **`scripts/benchmark/README.md`**.
- **Native-only baselines:** **`native_domain_inference_overview.md`** + `native_domain_da/` outputs.
