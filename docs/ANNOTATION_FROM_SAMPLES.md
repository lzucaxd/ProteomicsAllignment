# How MSstatsTMT annotation is built from sample files (CPTAC and CCLE)

MSstatsTMT needs one **annotation** table per study: each row is a **(Run, Channel, Fraction)** with **`BioReplicate`** (who is in that channel), **`Mixture`** (which TMT plex), **`Condition`** (`Norm` for the bridge / reference channel, `Sample` for biology), plus **`TechRepMixture`**. That table must have **exactly one `Norm` per mixture** (the bridge).

This page is the **one place** that explains how that table is produced from **sample design files** on the two domains. Implementation lives in **`data/pdc_psm_to_msstatsTMT_protein_matrix.R`** (CPTAC path) and **`data/ccle_peptide/ccle_to_msstats_input.py`** + the same R script with **`--msstats_input_dir`** (CCLE path).

**Not covered here:** subtype / tissue **benchmark** labels (`data/annotations/`, `data/biospecimen/`, curated `data/ccle/*.csv`) — those attach **gene × sample** metadata for tasks **after** matrices exist; they do not replace the MSstats design.

---

## CPTAC: from `*.sample.txt` to `annotation_filled_corrected.csv`

### What the CPTAC sample file is

CPTAC ships a tab-separated **`*.sample.txt`** per study. It has (at minimum):

- **`FileNameRegEx`** — pattern that matches **raw run** names in the PSM table so each LC–MS **Run** can be assigned to one **plex row**.
- **`AnalyticalSample`** — identifier for that plex (used as **`Mixture`** in the annotation).
- **One column per TMT channel** (e.g. `126`, `127N`, `131C`) — each cell is the **BioReplicate** label in that channel for that plex.
- **Exactly one channel per plex contains `POOL`** (case-insensitive) — that channel is the **bridge**; it becomes **`Condition = Norm`**; all other channels in that plex are **`Condition = Sample`**.

The R function **`load_sample_txt()`** reads this file, detects channel columns dynamically, and records which channel is POOL per plex.

### How the annotation rows are built

1. **PSMs are parsed** → long table with **`Run`**, **`Channel`**, **`Fraction`**, intensities, …
2. **Each unique `(Run, Channel, Fraction)`** from the PSMs is a candidate annotation row.
3. **`match_runs_to_plex()`** assigns each **Run** to the **plex row** whose **`FileNameRegEx`** matches the run string.
4. **`rebuild_annotation_from_sample_txt()`** fills **`BioReplicate`** from the sample file cell at that plex × channel, sets **`Mixture`** from **`AnalyticalSample`**, sets **`Condition`** to `Norm` if that cell is `POOL` else `Sample`, and sets **`TechRepMixture`** (typically 1).
5. Output is written as **`annotation_filled_corrected.csv`** under **`data/results/{study_id}/`** (and optionally overwrites **`annotation_filled.csv`** if you use replace flags). An **`annotation_audit.txt`** records whether the prior annotation matched the sample file or was rebuilt.

If there is **no** `--sample_txt`, the pipeline may use an existing **`--annotation`** CSV, **`--reference_channel`** to auto-fill Norm/Sample from the parsed channel grid only, or stop after writing a **template** — see **`data/PIPELINE_README.md`**.

### What you must supply (CPTAC)

- **`data/sample_files_msstats_tmt.csv`**: `study_id` → **`path`** to the real **`*.sample.txt`** on disk (often mirror-relative + **`CPTAC_LOCAL_MIRROR`**).
- Shell runners (e.g. **`run_msstats_tmt_gene_matrix.sh`**) resolve that path and pass **`--sample_txt`** into R when the file exists.

---

## CCLE: from Excel Sheet2 + peptide TSV to `annotation_filled.csv`

### Idea

CCLE peptide data does not use PDC **`sample.txt`**, but MSstatsTMT still needs the **same column schema**. A **Python converter** builds **`msstats_input.tsv`** and **`annotation_filled.csv`** from:

- **Peptide TSV** — reporter columns like `rq_126_sn`, … (channels auto-detected).
- **Sample sheet** — historically **Table S1 Sheet2** exported to CSV (`sample_info_ccle.csv`); columns map to the same roles as CPTAC.

### Column mapping (conceptual)

Documented in detail in **`data/ccle_peptide/README_CCLE.md`**; in short:

| Sample / design source | → | Annotation field |
|-------------------------|---|-------------------|
| Cell line / sample name | → | **BioReplicate** |
| Plex ID (e.g. Protein 10-Plex ID) | → | **Mixture** |
| TMT label (126, 127N, …) | → | **Channel** (canonicalized) |
| “Bridge line” in notes | → | **Condition = Norm**; others **Sample** |
| Run basename from peptide table | → | **Run** |
| (fixed) | → | **TechRepMixture = 1** |

The converter validates: **one Norm per mixture**, consistent channel sets, each **Run** maps to exactly one **Mixture**, no duplicate **(Mixture, Channel)**.

### What you run

1. Export Sheet2 → **`ccle_peptide/sample_info_ccle.csv`** (see **`export_sample_sheet2_csv.py`**).
2. **`ccle_to_msstats_input.py --tsv … --sample_csv … --outdir results/CCLE`** → writes **`msstats_input.tsv`** + **`annotation_filled.csv`**.
3. **`Rscript pdc_psm_to_msstatsTMT_protein_matrix.R --msstats_input_dir results/CCLE --outdir results/CCLE`** — **same R script** as CPTAC; it **skips PSM parsing** and uses the pre-built MSstats inputs, then runs summarization → **`gene_matrix.csv`**.

---

## Side-by-side

| Step | CPTAC | CCLE (peptide path) |
|------|--------|---------------------|
| Design file | Study **`*.sample.txt`** | **`sample_info_ccle.csv`** (from Sheet2) |
| Who builds MSstats rows | R (parse PSM) + R (annotation from sample.txt) | Python **`ccle_to_msstats_input.py`** |
| Annotation artifact | **`annotation_filled_corrected.csv`** (from sample.txt) | **`annotation_filled.csv`** (from converter) |
| Same R driver for summarization? | **`pdc_psm_to_msstatsTMT_protein_matrix.R`** | Yes, with **`--msstats_input_dir`** |

---

## Further reading (only if you need depth)

- **CPTAC pipeline flags and fallbacks:** [`../data/PIPELINE_README.md`](../data/PIPELINE_README.md)
- **CCLE converter steps:** [`../data/ccle_peptide/README_CCLE.md`](../data/ccle_peptide/README_CCLE.md)
- **File checklist:** [`../data/manifests/EXPECTED_INPUTS.md`](../data/manifests/EXPECTED_INPUTS.md)
