# CCLE peptide data → pipeline

**Same story as CPTAC, in one place:** conceptual overview of sample design → MSstats annotation for **both** domains is **[`../../docs/ANNOTATION_FROM_SAMPLES.md`](../../docs/ANNOTATION_FROM_SAMPLES.md)** (this file is the CCLE-specific runbook).

Use the **same** MSstatsTMT → gene matrix pipeline for CCLE peptide-level data without changing CPTAC (PDC) studies.

The pipeline supports **arbitrary TMT plex sizes** (TMT6, TMT10, TMT11, TMT16/TMTpro, TMT18). Channel columns are **automatically detected** from reporter-ion intensity fields (e.g. `rq_126_sn`, `rq_127n_sn`, `rq_127c_sn`). The annotation file follows the MSstatsTMT schema and requires **exactly one bridge (Condition = Norm) per mixture**.

## Files

- **Table_S1_Sample_Information (1).xlsx** – Sample info; **Sheet2** = channel annotation (same format as CPTAC sample design).
- **ccle_protein_quant_with_peptides_*.tsv** – Peptide-level TMT quantification (ProteinId, PeptideSequence, Charge, RunLoadPath, `rq_*_sn` columns; channels auto-detected).
- **export_sample_sheet2_csv.py** – Exports Sheet2 to `sample_info_ccle.csv` (stdlib only).
- **ccle_to_msstats_input.py** – Converts TSV + sample CSV → `msstats_input.tsv` + `annotation_filled.csv`; runs validation before writing.

## Sample file → MSstatsTMT annotation (same idea as CPTAC)

The annotation tells MSstatsTMT which biological sample is in each TMT channel in each plex (reference normalization, plex correction, protein summarization).

| Sample file column        | → | Annotation / meaning |
|---------------------------|---|------------------------|
| Cell Line                 | → | **BioReplicate** (sample ID) |
| Protein 10-Plex ID        | → | **Mixture** (TMT plex 0, 1, 2, …) |
| Protein TMT Label         | → | **Channel** (126, 127N, 127C, …; canonicalized) |
| Notes “Bridge line”       | → | **Condition** = Norm; else Sample |
| (from peptide table path) | → | **Run** = basename of RunLoadPath (no extension) |
| (fixed)                   | → | **TechRepMixture** = 1 |

- Each **Mixture** must have the **same channel set** (enforced by validation, unless `--allow-inconsistent-channels`).
- Each mixture must have **exactly one bridge (Norm)** (unless `--allow-multiple-norm`); MSstatsTMT uses it for reference normalization.
- Each **Run** maps to exactly one Mixture; (Mixture, Channel) is unique.

### Run → Mixture mapping

**Mixture IDs come from the sample file** (`Protein 10-Plex ID`). They are **not** inferred from run order.

- **Runs** correspond to LC–MS (fraction) files; each `RunLoadPath` basename is one Run.
- **Multiple runs may belong to the same Mixture**: several LC–MS fractions can be part of the same TMT plex.
- The converter assigns Run → Mixture by:
  1. **Run name**: if the run id contains `Prot_NN` (e.g. `Prot_01`, `Prot_11`), that identifies the plex (1-based in file: `Prot_01` → Mixture 0, `Prot_02` → Mixture 1, …).
  2. **Fallback**: if the set of channels with data in that run uniquely matches one plex’s channel set in the sample file, that plex is used.
- **Run order is never used** to assign plexes. Validation enforces that each Run maps to exactly one Mixture.

**Validation** (run before writing outputs): (1) **Channel consistency** — all mixtures use the same channel set (`unique(channel_set_per_mixture) == 1`); otherwise error: inconsistent TMT channel structure across mixtures. (2) Exactly one Norm per mixture. (3) Each Run maps to a single Mixture. (4) **(Mixture, Channel) unique** — no duplicate (Mixture, Channel) pairs (catches annotation corruption).

## Steps

### 1. Export sample sheet (once)

```bash
cd data
python3 ccle_peptide/export_sample_sheet2_csv.py
# -> ccle_peptide/sample_info_ccle.csv
```

### 2. Convert CCLE → MSstats input

```bash
cd data
python3 ccle_peptide/ccle_to_msstats_input.py \
  --tsv ccle_peptide/ccle_protein_quant_with_peptides_14745.tsv \
  --sample_csv ccle_peptide/sample_info_ccle.csv \
  --outdir results/CCLE
```

Writes `results/CCLE/msstats_input.tsv` and `results/CCLE/annotation_filled.csv`.

### 3. Run R pipeline (same script as CPTAC)

```bash
cd data
Rscript --no-init-file pdc_psm_to_msstatsTMT_protein_matrix.R \
  --msstats_input_dir results/CCLE \
  --outdir results/CCLE
```

This skips PSM parsing and annotation; it loads the pre-built input and runs protein summarization → gene matrix. Outputs: `protein_summary.tsv`, `gene_matrix.csv`, `plots/` (MSstatsTMT QC), `qc_summary.txt`.

## CPTAC unchanged

- CPTAC studies still use `--psm_dir`, `--sample_txt`, manifests, and `run_pipeline_per_manifest.sh` as before.
- Only when you pass **`--msstats_input_dir`** does the R script use the pre-built input path; otherwise it requires `--psm_dir`.
