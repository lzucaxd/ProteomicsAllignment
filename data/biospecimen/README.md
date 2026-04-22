# CPTAC biospecimen and clinical adjuncts (PDC000120 / breast)

Files in this folder are **small, curated exports** used by **`data/scripts/build_PDC000120_subtype_mapping.py`** to join PDC aliquot IDs to matrix columns and tumor / NAT / PAM50 context.

| File | Notes |
|------|--------|
| `PDC_study_biospecimen_03162026_190026.csv` | PDC biospecimen export (dated filename; replace with a fresh export if IDs drift). |
| `brca_cptac_2020_clinical_data.tsv` | CPTAC BRCA clinical table used in the mapping script. |
| `S039_*.xlsx` | Prospective / confirmatory study spreadsheets (label ↔ TMT channel mapping, follow-up). |

Do **not** commit patient-level data you are not allowed to share; these are aggregate / mapping tables as used in the lab pipeline.
