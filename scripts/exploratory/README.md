# New exploratory scripts (preferred location)

Put **new** one-off R/Python helpers here instead of under **`data/scripts/`**, so the layout stays:

- **`pipeline/psm_to_gene_matrix/`** + **`data/`** root → reproducible PDC → matrix path  
- **`scripts/benchmark/`** → harmonization benchmark  
- **`data/scripts/`** → frozen legacy / subtype / v1 (still in repo for provenance)

If a script grows into a maintained step, promote it into `scripts/benchmark/` or `src/harmonize/` with tests and docs.
