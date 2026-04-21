# Portable repository root for scripts under scripts/benchmark/
# Used by lab clones, GitHub checkouts, and CI.
#
# Resolution order:
#   1. Environment variable PROTEOMICS_ALIGNMENT_ROOT (must contain pyproject.toml)
#   2. Infer from Rscript --file= path (this file lives in scripts/benchmark/)
#   3. getwd() if it contains pyproject.toml
harmonize_repo_root <- function() {
  env <- Sys.getenv("PROTEOMICS_ALIGNMENT_ROOT", "")
  if (nzchar(env)) {
    root <- normalizePath(env, winslash = "/", mustWork = FALSE)
    if (dir.exists(root) && file.exists(file.path(root, "pyproject.toml"))) {
      return(root)
    }
    warning(
      "PROTEOMICS_ALIGNMENT_ROOT is set but is not this repository root (missing pyproject.toml): ",
      env
    )
  }
  ca <- commandArgs(trailingOnly = FALSE)
  fl <- ca[startsWith(ca, "--file=")]
  if (length(fl)) {
    sp <- normalizePath(sub("^--file=", "", fl[1L]))
    bench <- dirname(sp)
    root <- normalizePath(file.path(bench, "..", ".."), winslash = "/")
    if (file.exists(file.path(root, "pyproject.toml"))) {
      return(root)
    }
  }
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(wd, "pyproject.toml"))) {
    return(wd)
  }
  stop(
    "Cannot find repository root. Either:\n",
    "  export PROTEOMICS_ALIGNMENT_ROOT=/path/to/ProteomicsAllignment\n",
    "  or cd into the repository root before running Rscript.\n",
    call. = FALSE
  )
}
