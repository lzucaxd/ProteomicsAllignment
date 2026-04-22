# `scripts/`

**Manifest → PSM → gene matrix** is documented under **`../pipeline/psm_to_gene_matrix/README.md`** (executables live under **`../data/`**). **`../data/scripts/`** is exploratory / legacy — not that pipeline; see **`../data/scripts/README.md`**. New one-offs: **`exploratory/README.md`**.

| Path | Role |
|------|------|
| **`preprocessing/`** | Documentation index for PSM → `gene_matrix` (code still mostly under `data/`). See [preprocessing/README.md](preprocessing/README.md). |
| **`benchmark/`** | CPTAC–CCLE harmonization benchmark v2 (`run_overnight_v2.sh`). See [benchmark/README.md](benchmark/README.md). |
| **`presentation/`** | Build `presentation_materials/` for slides (`prepare_all.sh`). See [presentation/README.md](presentation/README.md). |
| **`methods/`** | Celligner wrapper (`run_celligner_representation.py`) |

Legacy exploratory / subtype helpers remain under **`../data/scripts/`** (many slide/report notes still cite those paths).

Example: `../data/scripts/notify_when_DA_finishes.sh`
