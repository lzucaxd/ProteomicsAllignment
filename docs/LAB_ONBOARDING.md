# Lab onboarding — clone, data layout, environment variables

This repository is meant to work on **any machine** after clone: no hardcoded home-directory paths in scripts under `scripts/`, `src/`, or portable `data/*.sh` helpers.

**Full pipeline (one doc):** **[`HOW_TO_RUN_EVERYTHING.md`](HOW_TO_RUN_EVERYTHING.md)**. For verify / commit policy only, see **[`CLEAN_CLONE_REPRODUCIBILITY.md`](CLEAN_CLONE_REPRODUCIBILITY.md)**.

---

## 1. Clone and Python

```bash
git clone <your-fork-or-upstream-url> ProteomicsAllignment
cd ProteomicsAllignment
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -U pip
pip install -r requirements.txt
pip install -e .
pip install -e ".[celligner]"   # optional: Celligner extras only
python3 scripts/verify_repro_setup.py
```

Benchmark / preprocessing Python code expects:

```bash
export PYTHONPATH="$PWD/src"
```

The overnight script `scripts/benchmark/run_overnight_v2.sh` sets this automatically.

---

## 2. R

```bash
Rscript install_r_packages.R
```

Use `export R_PROFILE_USER=/dev/null` if a parent directory’s broken `renv` would otherwise break `Rscript` subprocesses (the overnight script does this when appropriate).

---

## 3. CPTAC `.sample.txt` paths (`CPTAC_LOCAL_MIRROR`)

`data/sample_files_msstats_tmt.csv` lists each PDC study’s **path** to the CPTAC `*.sample.txt` design file. Paths in git are **relative** to a directory you choose:

| Resolution order | Meaning |
|------------------|---------|
| Absolute path in CSV | Used as-is if the file exists |
| Relative path | Tried under `data/` first |
| `CPTAC_LOCAL_MIRROR` | If set, `CPTAC_LOCAL_MIRROR/<path-in-CSV>` is tried |

**Typical layout:** your local CPTAC tree has folders `PDC000120/`, `PDC000153/`, … each containing that study’s `.sample.txt`. Point the mirror at the **parent** of those folders:

```bash
export CPTAC_LOCAL_MIRROR=/Volumes/lab/CPTAC/data   # example
```

Then a CSV path like `PDC000120/CPTAC2_Breast_....sample.txt` resolves to  
`/Volumes/lab/CPTAC/data/PDC000120/CPTAC2_Breast_....sample.txt`.

**Checks:**

```bash
cd data
python3 check_studies_sample_file.py
```

---

## 4. R benchmark scripts (`PROTEOMICS_ALIGNMENT_ROOT`)

R scripts under `scripts/benchmark/` resolve the repository root automatically from `Rscript --file=…`. If that ever fails (unusual wrappers), set:

```bash
export PROTEOMICS_ALIGNMENT_ROOT=/absolute/path/to/ProteomicsAllignment
```

Shared helper: `scripts/benchmark/harmonize_paths.R` (`harmonize_repo_root()`).

---

## 5. What not to commit

Large or personal artifacts are usually excluded via `.gitignore` (e.g. `data/results/`, `data/pdc_psm/`, `.venv/`). **Regenerated** benchmark CSVs and figures under `reports/benchmark_master/` are a **team policy** choice; many labs keep them local or on shared storage only.

---

## 6. Documentation map (GitHub renders Markdown)

| Doc | Audience |
|-----|----------|
| [README.md](../README.md) | First-time setup, pipeline entry |
| [docs/README.md](README.md) | Full documentation index (topic notes + legacy reports) |
| [END_TO_END_TECHNICAL_REPORT.md](END_TO_END_TECHNICAL_REPORT.md) | Paper / methods / full benchmark narrative |
| [BENCHMARK_V2_AND_PRESENTATION.md](BENCHMARK_V2_AND_PRESENTATION.md) | Slides + v2 outputs |
| [data/PIPELINE_README.md](../data/PIPELINE_README.md) | PDC → MSstatsTMT → matrix |
| [CONTRIBUTING.md](../CONTRIBUTING.md) | Forks, branches, PRs |

---

## 7. Optional: shared lab fork workflow

- Fork the upstream repo to your lab org; set **default branch** protection and require PRs for `main`.
- Document your lab’s **`CPTAC_LOCAL_MIRROR`** and shared **`data/results/`** location in the team wiki (not in git).
- Use **GitHub Issues** with labels (`data`, `benchmark`, `pipeline`) for reproducibility tickets.
