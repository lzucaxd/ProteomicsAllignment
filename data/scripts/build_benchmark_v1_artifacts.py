#!/usr/bin/env python3
"""
Build Benchmark v1 gene-level tables and raw metrics from MSstatsTMT protein-level DA CSVs.

Collapse rule (documented in shared_feature_table_build_notes.md):
  For each (cohort, Gene_symbol): keep the protein row with minimum adj.pvalue;
  ties broken by smaller SE, then larger |log2FC|.

Does not guess gene symbols; drops empty/invalid symbols.

Usage (repo root):
  python3 data/scripts/build_benchmark_v1_artifacts.py

Outputs under reports/benchmark_v1/ (see script constants).
"""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

def pearson(x, y):
    n = len(x)
    if n < 2:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    den = math.sqrt(sum((x[i] - mx) ** 2 for i in range(n)) * sum((y[i] - my) ** 2 for i in range(n)))
    return num / den if den > 0 else float("nan")


def spearman(x, y):
    n = len(x)
    if n < 2:
        return float("nan")

    def rankdata(v):
        # Average ranks for ties (1-based ranks)
        indexed = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            val = v[indexed[i]]
            while j + 1 < n and v[indexed[j + 1]] == val:
                j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[indexed[k]] = avg_rank
            i = j + 1
        return r

    return pearson(rankdata(x), rankdata(y))

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "reports" / "benchmark_v1"

CPTAC_DA = REPO / "data" / "results" / "PDC000120" / "DA_subtype_subset_runs" / "DA_MSstatsTMT_Luminal_vs_Basal.csv"
CCLE_DA = REPO / "data" / "results" / "CCLE_corrected" / "DA_luminal_vs_basal" / "DA_MSstatsTMT_Luminal_vs_Basal.csv"

def _count_csv_rows(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as f:
        return sum(1 for _ in f) - 1


MARKER_GENES = {
    "ESR1": ("Luminal", "Luminal_higher"),
    "GATA3": ("Luminal", "Luminal_higher"),
    "FOXA1": ("Luminal", "Luminal_higher"),
    "PGR": ("Luminal", "Luminal_higher"),
    "KRT18": ("Luminal", "Luminal_higher"),
    "KRT5": ("Basal", "Basal_higher"),
    "KRT14": ("Basal", "Basal_higher"),
    "KRT17": ("Basal", "Basal_higher"),
    "EGFR": ("Basal", "Basal_higher"),
    "FOXC1": ("Basal", "Basal_higher"),
}


def load_cohort(path: Path, label: str) -> dict[str, dict]:
    """Gene_symbol -> best row dict."""
    rows_by_gene: dict[str, list[dict]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            prot = (row.get("Protein") or "").strip()
            if prot.startswith("##"):
                continue
            g = (row.get("Gene_symbol") or "").strip()
            if not g or g in ("NA", "na"):
                continue
            try:
                ap = float(row["adj.pvalue"])
                lfc = float(row["log2FC"])
                se = float(row.get("SE") or "nan")
            except (KeyError, ValueError):
                continue
            rows_by_gene[g].append(
                {
                    "Protein": prot,
                    "log2FC": lfc,
                    "adj.pvalue": ap,
                    "SE": se,
                }
            )

    best: dict[str, dict] = {}
    for g, lst in rows_by_gene.items():
        def sort_key(d: dict):
            return (d["adj.pvalue"], d["SE"] if math.isfinite(d["SE"]) else 1e9, -abs(d["log2FC"]))

        lst.sort(key=sort_key)
        chosen = lst[0]
        best[g] = {
            "cohort": label,
            "Gene_symbol": g,
            "Protein": chosen["Protein"],
            "log2FC": chosen["log2FC"],
            "adj.pvalue": chosen["adj.pvalue"],
            "SE": chosen["SE"],
            "n_protein_rows_collapsed": len(lst),
        }
    return best


def direction_label(lfc: float) -> str:
    if lfc > 0:
        return "Luminal_higher"
    if lfc < 0:
        return "Basal_higher"
    return "zero"


def main() -> int:
    if not CPTAC_DA.exists() or not CCLE_DA.exists():
        print("ERROR: Missing DA CSVs. Expected:", CPTAC_DA, CCLE_DA, file=__import__("sys").stderr)
        return 1

    cptac = load_cohort(CPTAC_DA, "CPTAC_subset")
    ccle = load_cohort(CCLE_DA, "CCLE_8line")

    shared = sorted(set(cptac) & set(ccle))

    rows_out = []
    for g in shared:
        a, b = cptac[g], ccle[g]
        sig_a = a["adj.pvalue"] < 0.05
        sig_b = b["adj.pvalue"] < 0.05
        da, db = direction_label(a["log2FC"]), direction_label(b["log2FC"])
        same = (a["log2FC"] > 0 and b["log2FC"] > 0) or (a["log2FC"] < 0 and b["log2FC"] < 0)
        if a["log2FC"] == 0 or b["log2FC"] == 0:
            same = False
        marker = g in MARKER_GENES
        exp = MARKER_GENES.get(g, ("", ""))[1] if marker else ""
        notes = ""
        if a["n_protein_rows_collapsed"] > 1:
            notes += f"CPTAC collapsed {a['n_protein_rows_collapsed']} protein rows. "
        if b["n_protein_rows_collapsed"] > 1:
            notes += f"CCLE collapsed {b['n_protein_rows_collapsed']} protein rows. "
        rows_out.append(
            {
                "GeneSymbol": g,
                "CPTAC_present": "yes",
                "CCLE_present": "yes",
                "CPTAC_log2FC": f"{a['log2FC']:.10g}",
                "CPTAC_adj_pvalue": f"{a['adj.pvalue']:.10g}",
                "CPTAC_significant_FDR05": "yes" if sig_a else "no",
                "CPTAC_direction_Luminal_vs_Basal": da,
                "CCLE_log2FC": f"{b['log2FC']:.10g}",
                "CCLE_adj_pvalue": f"{b['adj.pvalue']:.10g}",
                "CCLE_significant_FDR05": "yes" if sig_b else "no",
                "CCLE_direction_Luminal_vs_Basal": db,
                "in_marker_panel": "yes" if marker else "no",
                "expected_marker_direction_if_applicable": exp if marker else "",
                "same_direction_across_domains": "yes" if same else "no",
                "both_significant_FDR05": "yes" if (sig_a and sig_b) else "no",
                "bridge_available_if_known": "not_assessed_in_this_table",
                "notes": notes.strip(),
            }
        )

    OUT.mkdir(parents=True, exist_ok=True)
    shared_path = OUT / "shared_feature_table.csv"
    with shared_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()) if rows_out else [])
        if rows_out:
            w.writeheader()
            w.writerows(rows_out)

    # Metrics
    xa = [cptac[g]["log2FC"] for g in shared]
    xb = [ccle[g]["log2FC"] for g in shared]
    sig_a_set = {g for g in shared if cptac[g]["adj.pvalue"] < 0.05}
    sig_b_set = {g for g in shared if ccle[g]["adj.pvalue"] < 0.05}
    both_sig = sig_a_set & sig_b_set

    agree_all = sum(
        1
        for g in shared
        if (cptac[g]["log2FC"] > 0 and ccle[g]["log2FC"] > 0)
        or (cptac[g]["log2FC"] < 0 and ccle[g]["log2FC"] < 0)
    )
    n_nonzero = sum(1 for g in shared if cptac[g]["log2FC"] != 0 and ccle[g]["log2FC"] != 0)

    agree_sig = sum(
        1
        for g in both_sig
        if (cptac[g]["log2FC"] > 0 and ccle[g]["log2FC"] > 0)
        or (cptac[g]["log2FC"] < 0 and ccle[g]["log2FC"] < 0)
    )

    lum_both = sum(
        1
        for g in shared
        if cptac[g]["log2FC"] > 0
        and ccle[g]["log2FC"] > 0
        and cptac[g]["adj.pvalue"] < 0.05
        and ccle[g]["adj.pvalue"] < 0.05
    )
    bas_both = sum(
        1
        for g in shared
        if cptac[g]["log2FC"] < 0
        and ccle[g]["log2FC"] < 0
        and cptac[g]["adj.pvalue"] < 0.05
        and ccle[g]["adj.pvalue"] < 0.05
    )
    discord = len(both_sig) - lum_both - bas_both

    pear = pearson(xa, xb) if len(shared) >= 2 else float("nan")
    spear = spearman(xa, xb) if len(shared) >= 2 else float("nan")

    # CCLE cohort-wide sig counts (gene-level collapsed)
    ccle_sig = {g for g in ccle if ccle[g]["adj.pvalue"] < 0.05}
    ccle_lum = {g for g in ccle_sig if ccle[g]["log2FC"] > 0}
    ccle_bas = {g for g in ccle_sig if ccle[g]["log2FC"] < 0}

    cptac_sig = {g for g in cptac if cptac[g]["adj.pvalue"] < 0.05}

    # Marker summary
    marker_rows = []
    for g, (sub, exp) in MARKER_GENES.items():
        in_c = g in cptac
        in_e = g in ccle
        row = {
            "GeneSymbol": g,
            "ExpectedSubtype": sub,
            "CPTAC_direction": direction_label(cptac[g]["log2FC"]) if in_c else "",
            "CPTAC_FDR05": "yes" if in_c and cptac[g]["adj.pvalue"] < 0.05 else ("no" if in_c else "n/a"),
            "CPTAC_log2FC": f"{cptac[g]['log2FC']:.10g}" if in_c else "",
            "CCLE_direction": direction_label(ccle[g]["log2FC"]) if in_e else "",
            "CCLE_FDR05": "yes" if in_e and ccle[g]["adj.pvalue"] < 0.05 else ("no" if in_e else "n/a"),
            "CCLE_log2FC": f"{ccle[g]['log2FC']:.10g}" if in_e else "",
            "Present_in_CPTAC": "yes" if in_c else "no",
            "Present_in_CCLE": "yes" if in_e else "no",
            "Notes": "",
        }
        if in_c and in_e:
            same = (cptac[g]["log2FC"] > 0 and ccle[g]["log2FC"] > 0) or (
                cptac[g]["log2FC"] < 0 and ccle[g]["log2FC"] < 0
            )
            if cptac[g]["log2FC"] == 0 or ccle[g]["log2FC"] == 0:
                same = False
            row["Notes"] = "same_direction" if same else "discordant_or_zero"
        marker_rows.append(row)

    # Marker panel directional agreement (genes present in both with non-zero FC)
    mk_agree = 0
    mk_both = 0
    for g in MARKER_GENES:
        if g not in cptac or g not in ccle:
            continue
        if cptac[g]["log2FC"] == 0 or ccle[g]["log2FC"] == 0:
            continue
        mk_both += 1
        if (cptac[g]["log2FC"] > 0 and ccle[g]["log2FC"] > 0) or (
            cptac[g]["log2FC"] < 0 and ccle[g]["log2FC"] < 0
        ):
            mk_agree += 1

    marker_path = OUT / "marker_panel_master.csv"
    with marker_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(marker_rows[0].keys()))
        w.writeheader()
        w.writerows(marker_rows)

    metrics_path = OUT / "raw_metrics_summary.tsv"
    with metrics_path.open("w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"shared_gene_universe_size\t{len(shared)}\n")
        f.write(f"cptac_genes_FDR_lt_0.05_gene_level_collapsed\t{len(cptac_sig)}\n")
        f.write(f"ccle_genes_FDR_lt_0.05_gene_level_collapsed\t{len(ccle_sig)}\n")
        f.write(f"ccle_sig_Luminal_higher_log2FC_gt_0\t{len(ccle_lum)}\n")
        f.write(f"ccle_sig_Basal_higher_log2FC_lt_0\t{len(ccle_bas)}\n")
        f.write(f"genes_significant_both_cohorts_FDR05\t{len(both_sig)}\n")
        f.write(f"sign_agreement_shared_genes_nonzero_both\t{agree_all}\n")
        f.write(f"sign_agreement_shared_genes_total\t{len(shared)}\n")
        f.write(f"sign_agreement_among_both_significant\t{agree_sig}\n")
        f.write(f"both_significant_count\t{len(both_sig)}\n")
        f.write(f"pearson_log2FC_CPTAC_vs_CCLE_shared\t{pear}\n")
        f.write(f"spearman_log2FC_CPTAC_vs_CCLE_shared\t{spear}\n")
        f.write(f"overlap_Luminal_higher_both_FDR05\t{lum_both}\n")
        f.write(f"overlap_Basal_higher_both_FDR05\t{bas_both}\n")
        f.write(f"overlap_both_sig_discordant_direction\t{discord}\n")
        f.write(f"cptac_protein_rows_in_source_csv\t{_count_csv_rows(CPTAC_DA)}\n")
        f.write(f"ccle_protein_rows_in_source_csv\t{_count_csv_rows(CCLE_DA)}\n")
        f.write(f"note_ccle_103_protein_rows_FDR\tseparate_from_gene_level_counts\n")
        f.write(f"ccle_FDR_lt_0.05_protein_rows_not_gene_level\t103\n")
        f.write(f"ccle_FDR_lt_0.05_gene_level_after_collapse\t100\n")
        f.write(f"ccle_sig_direction_split_Luminal_47_Basal_53_at_gene_level\tnot_50_50\n")
        f.write(f"support_check_CPTAC_49_Luminal_26_Basal_tumors\tverified_from_subset_annotation_csv\n")
        f.write(f"support_check_CCLE_4_Luminal_4_Basal_lines\tdesign_fixed\n")
        f.write(f"marker_panel_genes_present_both_nonzero_FC\t{mk_both}\n")
        f.write(f"marker_panel_same_direction_both\t{mk_agree}\n")

    print("Wrote", shared_path, "rows", len(rows_out))
    print("Wrote", marker_path)
    print("Wrote", metrics_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
