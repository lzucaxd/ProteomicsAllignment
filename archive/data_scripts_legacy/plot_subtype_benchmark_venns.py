#!/usr/bin/env python3
"""
Venn diagrams for subtype benchmark slides.

- CPTAC full vs subset: FDR < 0.05 genes (Gene_symbol) from MSstatsTMT CSVs.
- CPTAC subset vs CCLE: intersection of significant genes (by symbol).
- Canonical panel: 10 genes — coverage in CPTAC subset vs CCLE protein table.

Usage (repo root; needs matplotlib + matplotlib-venn):
  python3 data/scripts/plot_subtype_benchmark_venns.py [output_dir]

Default output_dir: reports/presentation_subtype_benchmark/03_ccle_subtype/venn_figures/
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path


def repo_root() -> Path:
    p = Path(__file__).resolve().parents[2]
    if (p / "data" / "scripts" / Path(__file__).name).exists():
        return p
    return Path(__file__).resolve().parents[1]


def count_sig_proteins_msstats_csv(path: Path, fdr_max: float = 0.05) -> int:
    n = 0
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                ap = float(row.get("adj.pvalue", "") or "nan")
            except ValueError:
                continue
            if ap < fdr_max:
                n += 1
    return n


def read_sig_genes_from_msstats_csv(path: Path, fdr_max: float = 0.05) -> set[str]:
    """MSstatsTMT groupComparison CSV: adj.pvalue column, Gene_symbol."""
    out: set[str] = set()
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f)
        for row in r:
            sym = (row.get("Gene_symbol") or "").strip()
            if not sym or sym in ('""', "NA", "na"):
                continue
            try:
                ap = float(row.get("adj.pvalue", "") or row.get("adj.P.Val", "") or "nan")
            except ValueError:
                continue
            if ap < fdr_max:
                out.add(sym)
    return out


def read_sig_genes_ccle_gene_csv(path: Path, fdr_max: float = 0.05) -> set[str]:
    """DA_MSstatsTMT_*_gene_symbols.csv with FDR column (may be 0 or scientific)."""
    out: set[str] = set()
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f)
        for row in r:
            sym = (row.get("Gene_symbol") or "").strip().strip('"')
            if not sym:
                continue
            raw = row.get("FDR", "")
            try:
                ap = float(raw) if raw not in ("", "NA") else float("nan")
            except ValueError:
                continue
            if ap < fdr_max:
                out.add(sym)
    return out


def read_canonical_panel_genes() -> list[str]:
    return [
        "ESR1",
        "GATA3",
        "FOXA1",
        "PGR",
        "KRT18",
        "KRT5",
        "KRT14",
        "KRT17",
        "EGFR",
        "FOXC1",
    ]


def genes_in_ccle_any_row(protein_csv: Path) -> set[str]:
    """All Gene_symbol observed in CCLE DA protein-level CSV (tested)."""
    out: set[str] = set()
    with protein_csv.open(newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f)
        for row in r:
            sym = (row.get("Gene_symbol") or "").strip()
            if sym and not sym.startswith("##"):
                out.add(sym)
    return out


def plot_two(ax_title: str, a: set[str], b: set[str], label_a: str, label_b: str, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    ab = a & b
    only_a = a - b
    only_b = b - a
    plt.figure(figsize=(7, 6))
    v = venn2(
        subsets=(len(only_a), len(only_b), len(ab)),
        set_labels=(label_a, label_b),
    )
    if v.get_patch_by_id("10"):
        v.get_patch_by_id("10").set_color("#b3cde3")
    if v.get_patch_by_id("01"):
        v.get_patch_by_id("01").set_color("#fdda95")
    if v.get_patch_by_id("11"):
        v.get_patch_by_id("11").set_color("#ccebc5")
    for text in v.set_labels:
        if text:
            text.set_fontsize(12)
    for text in v.subset_labels:
        if text:
            text.set_fontsize(14)
    plt.title(ax_title, fontsize=13)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_canonical_coverage(
    cptac_sig: set[str],
    ccle_tested: set[str],
    panel: list[str],
    out_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    in_cptac = {g for g in panel if g in cptac_sig}
    in_ccle = {g for g in panel if g in ccle_tested}
    only_cptac = in_cptac - in_ccle
    only_ccle = in_ccle - in_cptac
    both = in_cptac & in_ccle
    plt.figure(figsize=(7.5, 6))
    v = venn2(
        subsets=(len(only_cptac), len(only_ccle), len(both)),
        set_labels=("Sig. in CPTAC\nsubset (panel genes)", "Tested in CCLE\n(panel genes)"),
    )
    if v.get_patch_by_id("10"):
        v.get_patch_by_id("10").set_color("#decbe4")
    if v.get_patch_by_id("01"):
        v.get_patch_by_id("01").set_color("#fed9a6")
    if v.get_patch_by_id("11"):
        v.get_patch_by_id("11").set_color("#b3e2ad")
    plt.title("Canonical marker genes (n=10): CPTAC FDR hit vs CCLE coverage", fontsize=12)
    miss_both = [g for g in panel if g not in cptac_sig and g not in ccle_tested]
    miss_line = f"Panel not sig. CPTAC & absent CCLE: {', '.join(miss_both) if miss_both else 'none'}"
    plt.annotate(miss_line, xy=(0.5, -0.12), xycoords="axes fraction", ha="center", fontsize=9, wrap=True)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def main() -> int:
    root = repo_root()
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "out_dir",
        nargs="?",
        default=str(
            root
            / "reports"
            / "presentation_subtype_benchmark"
            / "03_ccle_subtype"
            / "venn_figures"
        ),
    )
    args = ap.parse_args()
    out_dir = Path(args.out_dir)

    cptac_full = root / "data" / "results" / "PDC000120" / "DA_MSstatsTMT_Basal_vs_Luminal.csv"
    cptac_subset = root / "data" / "results" / "PDC000120" / "DA_subtype_subset_runs" / "DA_MSstatsTMT_Luminal_vs_Basal.csv"
    ccle = root / "data" / "results" / "CCLE_corrected" / "DA_luminal_vs_basal" / "DA_MSstatsTMT_Luminal_vs_Basal.csv"

    missing = [str(p) for p in (cptac_full, cptac_subset, ccle) if not p.exists()]
    if missing:
        print("Missing inputs (run from machine with data/results checked in or restored):", file=sys.stderr)
        for m in missing:
            print(" ", m, file=sys.stderr)
        return 1

    sig_full = read_sig_genes_from_msstats_csv(cptac_full)
    sig_subset = read_sig_genes_from_msstats_csv(cptac_subset)
    sig_ccle = read_sig_genes_from_msstats_csv(ccle)
    ccle_sig_proteins = count_sig_proteins_msstats_csv(ccle)
    ccle_genes_tested = genes_in_ccle_any_row(ccle)

    plot_two(
        "FDR < 0.05 genes: CPTAC full cohort vs mixture-balanced subset\n(same tumor labels; subtype contrast)",
        sig_full,
        sig_subset,
        f"CPTAC full\nn={len(sig_full)}",
        f"CPTAC subset\nn={len(sig_subset)}",
        out_dir / "venn_CPTAC_full_vs_subset_FDR05_genes.png",
    )

    plot_two(
        "FDR < 0.05: CPTAC subset vs CCLE (deduplicated Gene_symbol)\n"
        f"CCLE side = {len(sig_ccle)} genes from {ccle_sig_proteins} significant protein rows",
        sig_subset,
        sig_ccle,
        f"CPTAC subset\nn={len(sig_subset)} genes",
        f"CCLE\nn={len(sig_ccle)} genes\n({ccle_sig_proteins} proteins)",
        out_dir / "venn_CPTAC_subset_vs_CCLE_FDR05_genes.png",
    )

    panel = read_canonical_panel_genes()
    plot_canonical_coverage(
        sig_subset,
        ccle_genes_tested,
        panel,
        out_dir / "venn_canonical_panel_coverage_CPTACsig_vs_CCLEtested.png",
    )

    summary = out_dir / "venn_counts_summary.txt"
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", encoding="utf-8") as f:
        f.write("Subtype benchmark Venn inputs (MSstatsTMT adj.pvalue < 0.05)\n")
        f.write("=" * 60 + "\n")
        f.write(f"CPTAC full unique Gene_symbol: {len(sig_full)}\n")
        f.write(f"CPTAC subset unique Gene_symbol: {len(sig_subset)}\n")
        f.write(f"CCLE unique Gene_symbol (FDR<0.05): {len(sig_ccle)}\n")
        f.write(
            f"CCLE significant protein rows (FDR<0.05): {ccle_sig_proteins}\n"
            "  → Venns use **genes** (empty Gene_symbol dropped; isoforms collapsed).\n"
            "  Example: 103 protein rows vs 100 genes — MAPT x3 + one row w/o symbol.\n"
        )
        f.write(f"Intersection full ∩ subset: {len(sig_full & sig_subset)}\n")
        f.write(f"Intersection subset ∩ CCLE: {len(sig_subset & sig_ccle)}\n")
        f.write("\nCanonical panel (10 genes)\n")
        for g in panel:
            f.write(
                f"  {g}: CPTAC_subset_sig={g in sig_subset}, CCLE_tested={g in ccle_genes_tested}\n"
            )

    print("Wrote:", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
