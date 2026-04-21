#!/usr/bin/env python3
"""
Convert CCLE peptide-level TSV + sample info (Sheet2) to MSstatsTMT input format.
Outputs msstats_input.tsv and annotation_filled.csv so the main R pipeline can
run with --msstats_input_dir (skips PSM parsing; CPTAC flow unchanged).

Design (corrected):
- Protein 10-Plex ID = 0 rows describe pooled bridge composition only; they are
  NOT a measured TMT mixture and must not become Mixture 0.
- For each real plex (ID >= 1), the sample sheet lists 9 cell-line channels; the
  remaining TMT10 channel (inferred as missing from the full reporter set) is the
  bridge/reference — one annotation row with Condition=Norm and BioReplicate=POOL.
- Run names contain Prot_NN; for this dataset NN matches Protein 10-Plex ID (1..42).

Legacy behavior (pre-fix) is preserved in ccle_to_msstats_input.py.bak_pre_bridge_fix
for diff/audit; the old script used Prot_NN→(NN-1) and treated ID 0 as a mixture.

Usage (from data/ or repo root):
  python3 ccle_peptide/ccle_to_msstats_input.py --tsv ccle_peptide/ccle_protein_quant_with_peptides_14745.tsv \\
    --sample_csv ccle_peptide/sample_info_ccle.csv --outdir results/CCLE_corrected

Prerequisite: export sample sheet2 first:
  python3 ccle_peptide/export_sample_sheet2_csv.py
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

# Reporter-ion columns: rq_126_sn, rq_127n_sn, rq_127c_sn, ... (auto-detected)
RQ_SN_PATTERN = re.compile(r"^rq_(\d+)([nc])?_sn$", re.IGNORECASE)

# e.g. f05449_Prot_09_F09 → fraction 9
FRACTION_RUN_PATTERN = re.compile(r"_F(\d+)(?:_|$)", re.IGNORECASE)

RUN_PLEX_PATTERN = re.compile(r"Prot_(\d+)", re.IGNORECASE)

# Backup of previous revision: ccle_to_msstats_input.py.bak_pre_bridge_fix (repo)


def channel_sort_key(label):
    """Sort key for TMT channels: 126, 127N, 127C, ..."""
    m = re.match(r"^(\d+)([NC])?$", str(label).upper())
    if not m:
        return (0, "Z")
    num = int(m.group(1))
    suffix = m.group(2) or ""
    suf_order = {"": 0, "N": 1, "C": 2}
    return (num, suf_order.get(suffix, 99))


def detect_channel_columns(fieldnames):
    """Return [(col, channel_label), ...] in canonical order."""
    found = []
    for name in fieldnames:
        m = RQ_SN_PATTERN.match(name.strip())
        if not m:
            continue
        num, suf = m.group(1), (m.group(2) or "").upper()
        label = num + suf
        found.append((name, label))
    found.sort(key=lambda x: channel_sort_key(x[1]))
    return found


def normalize_channel(ch):
    """Sample file may have 127n/127c; MSstatsTMT expects 127N, 127C."""
    if not ch:
        return ch
    s = str(ch).strip()
    m = re.match(r"^(\d+)([nc])?$", s, re.IGNORECASE)
    if m:
        return m.group(1) + (m.group(2) or "").upper()
    if re.match(r"^\d+[NC]$", s, re.IGNORECASE):
        return s[: len(s) - 1] + s[-1].upper()
    return s


def load_sample_info(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append({k.strip(): (v.strip() if v else "") for k, v in row.items()})
    return rows


def load_bridge_metadata_rows(sample_rows):
    """Plex ID 0 rows: bridge pool composition (metadata only)."""
    key_col = "Protein 10-Plex ID"
    out = []
    for row in sample_rows:
        try:
            plex = int(row.get(key_col, ""))
        except (ValueError, TypeError):
            continue
        if plex != 0:
            continue
        out.append(row)
    return out


def build_real_plex_channel_annotation(sample_rows, full_channel_set):
    """
    Build (PlexID, Channel) -> (BioReplicate, Condition) for real plexes only (ID >= 1).
    Adds one inferred bridge row per plex: missing TMT channel -> (POOL, Norm).

    full_channel_set: set of canonical channel labels from the peptide TSV (rq_*).
    """
    key_col = "Protein 10-Plex ID"
    ch_col = "Protein TMT Label"
    cell_col = "Cell Line"
    notes_col = "Notes"
    ann = {}
    real_plex_ids = set()

    for row in sample_rows:
        try:
            plex = int(row.get(key_col, ""))
        except (ValueError, TypeError):
            continue
        if plex < 1:
            continue
        ch_raw = (row.get(ch_col) or "").strip()
        ch = normalize_channel(ch_raw)
        if not ch:
            continue
        cell = (row.get(cell_col) or "").strip() or f"Plex{plex}_{ch}"
        notes = (row.get(notes_col) or "").lower()
        if "bridge" in notes:
            cond = "Norm"
        else:
            cond = "Sample"
        ann[(plex, ch)] = (cell, cond)
        real_plex_ids.add(plex)

    inconsistencies = []
    bridge_labels_added = {}

    for plex in sorted(real_plex_ids):
        listed = {ch for (pp, ch) in ann.keys() if pp == plex}
        missing = full_channel_set - listed
        if len(missing) == 0:
            inconsistencies.append(
                f"Plex {plex}: no missing channel (listed covers full reporter set); bridge not added."
            )
            continue
        if len(missing) > 1:
            inconsistencies.append(
                f"Plex {plex}: multiple missing channels {sorted(missing, key=channel_sort_key)} — check sample sheet vs TSV."
            )
            # Conservative: do not guess; skip adding bridge for this plex
            continue
        (bridge_ch,) = tuple(missing)
        key = (plex, bridge_ch)
        if key in ann:
            inconsistencies.append(f"Plex {plex}: bridge channel {bridge_ch} already present.")
            continue
        ann[key] = ("POOL", "Norm")
        bridge_labels_added[plex] = bridge_ch

    return ann, sorted(real_plex_ids), inconsistencies, bridge_labels_added


def run_id_from_path(path):
    base = os.path.basename(path)
    return base.replace(".mzXML", "").replace(".raw", "").replace(".mzML", "").strip()


def mixture_from_run_id(run_id, valid_plex_ids):
    """
    Infer Mixture (= Protein 10-Plex ID) from run_id.
    CCLE files use Prot_NN with NN matching the sample sheet plex ID (1..42 here).
    """
    m = RUN_PLEX_PATTERN.search(run_id)
    if not m:
        return None
    num = int(m.group(1))
    if num in valid_plex_ids:
        return num
    return None


def parse_fraction_from_run_id(run_id):
    """Return positive int if _F## found, else None."""
    m = FRACTION_RUN_PATTERN.search(run_id)
    if not m:
        return None
    return int(m.group(1))


def build_plex_metadata(sample_rows, full_channel_set):
    plex_ch_ann, real_plex_ids, inconsistencies, bridge_labels = build_real_plex_channel_annotation(
        sample_rows, full_channel_set
    )
    plex_ids = set(real_plex_ids)
    channels_per_plex = {}
    norm_count_per_plex = {}
    for (p, ch), (bio, cond) in plex_ch_ann.items():
        channels_per_plex.setdefault(p, set()).add(ch)
        if cond == "Norm":
            norm_count_per_plex[p] = norm_count_per_plex.get(p, 0) + 1
    return plex_ch_ann, plex_ids, channels_per_plex, norm_count_per_plex, inconsistencies, bridge_labels


def validate_annotation(annotation_rows, plex_ch_ann, run_id_to_mixture, allow_inconsistent_channels=False, allow_multiple_norm=False):
    if not annotation_rows:
        raise SystemExit("Validation failed: annotation is empty.")

    seen = set()
    for r in annotation_rows:
        key = (r["Run"], r["Channel"])
        if key in seen:
            raise SystemExit(
                f"Validation failed: (Run, Channel) must be unique. Duplicate: Run {key[0]}, Channel {key[1]}."
            )
        seen.add(key)

    by_mix = {}
    for r in annotation_rows:
        m = r["Mixture"]
        if m not in by_mix:
            by_mix[m] = {"channels": set(), "norm_channels": set()}
        by_mix[m]["channels"].add(r["Channel"])
        if (r.get("Condition") or "").strip().lower() == "norm":
            by_mix[m]["norm_channels"].add(r["Channel"])

    ch_sets = [frozenset(v["channels"]) for v in by_mix.values()]
    if len(set(ch_sets)) != 1:
        if not allow_inconsistent_channels:
            raise SystemExit(
                "Validation failed: inconsistent TMT channel structure across mixtures. "
                "Use --allow-inconsistent-channels only for exceptional debugging."
            )
        print("Warning: Inconsistent TMT channel structure across mixtures (--allow-inconsistent-channels).")

    # One reference *channel label* per mixture (not one row: many runs repeat the same channel).
    for m, v in by_mix.items():
        nref = len(v["norm_channels"])
        if nref != 1 and not allow_multiple_norm:
            raise SystemExit(
                f"Validation failed: Mixture {m} has Norm on {nref} distinct channel label(s): {sorted(v['norm_channels'])}. "
                "MSstatsTMT expects exactly one bridge channel per mixture."
            )
    if allow_multiple_norm and any(len(by_mix[m]["norm_channels"]) != 1 for m in by_mix):
        print("Warning: At least one mixture has multiple Norm channel labels (--allow-multiple-norm).")

    run_to_mix = {}
    for r in annotation_rows:
        run = r["Run"]
        mix = r["Mixture"]
        if run in run_to_mix and run_to_mix[run] != mix:
            raise SystemExit(
                f"Validation failed: Run must map to exactly one Mixture. Run '{run}' maps to both {run_to_mix[run]} and {mix}."
            )
        run_to_mix[run] = mix
    print("Validation passed: channels consistent, one Norm per mixture, Run->single Mixture, (Run,Channel) unique.")


def main():
    ap = argparse.ArgumentParser(description="CCLE peptide TSV + sample sheet -> MSstatsTMT input (bridge-corrected)")
    ap.add_argument("--tsv", required=True, help="CCLE peptide TSV")
    ap.add_argument("--sample_csv", required=True, help="Sample info CSV from Sheet2")
    ap.add_argument("--outdir", default="results/CCLE", help="Output directory")
    ap.add_argument("--allow-inconsistent-channels", action="store_true", help="Allow different channel sets per mixture (debug only)")
    ap.add_argument("--allow-multiple-norm", action="store_true", help="Allow more than one Norm per mixture (debug only)")
    ap.add_argument("--summary-file", default="converter_summary.txt", help="Written under outdir")
    args = ap.parse_args()

    if not os.path.isfile(args.tsv):
        raise SystemExit(f"Not found: {args.tsv}")
    if not os.path.isfile(args.sample_csv):
        raise SystemExit(f"Not found: {args.sample_csv}")

    sample_rows = load_sample_info(args.sample_csv)
    bridge_meta = load_bridge_metadata_rows(sample_rows)

    with open(args.tsv, newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = list(r.fieldnames or [])
    channel_cols = detect_channel_columns(fieldnames)
    if not channel_cols:
        raise SystemExit(
            "No reporter-ion columns found. Expected rq_*_sn (e.g. rq_126_sn, rq_127n_sn)."
        )
    full_channel_set = {ch for _, ch in channel_cols}
    channel_col_names = [c for c, _ in channel_cols]
    required = ["ProteinId", "PeptideSequence", "Charge", "RunLoadPath"] + channel_col_names
    for col in required:
        if col not in fieldnames:
            raise SystemExit(f"TSV missing column: {col}")

    plex_ch_ann, plex_ids, channels_per_plex, norm_count_per_plex, inconsistencies, bridge_labels = build_plex_metadata(
        sample_rows, full_channel_set
    )
    if not plex_ids:
        raise SystemExit("Sample file has no valid real plex rows (Protein 10-Plex ID >= 1 with channels).")

    channel_set_to_plexes = defaultdict(list)
    for p in plex_ids:
        channel_set_to_plexes[frozenset(channels_per_plex[p])].append(p)

    os.makedirs(args.outdir, exist_ok=True)

    n_plexes = len(plex_ids)
    size_to_plexes = defaultdict(list)
    for p in plex_ids:
        size_to_plexes[len(channels_per_plex[p])].append(p)
    print("Real plexes (ID >= 1):", n_plexes)
    if len(size_to_plexes) == 1:
        (n_ch,) = size_to_plexes.keys()
        print("Channels per plex (after bridge inference):", n_ch)
    else:
        for s, plist in sorted(size_to_plexes.items()):
            print(f"  {s} channels: plexes {sorted(plist)[:8]}{'...' if len(plist) > 8 else ''}")
    print("Bridge metadata rows (ID=0, not used as mixture):", len(bridge_meta))
    print("Inferred bridge rows added (one per plex where possible):", len(bridge_labels))
    print("Detected", len(channel_cols), "reporter columns in TSV:", [ch for _, ch in channel_cols])

    def resolve_mixture(run_id, row):
        if run_id in run_id_to_mixture:
            return run_id_to_mixture[run_id]
        mixture = mixture_from_run_id(run_id, plex_ids)
        if mixture is None:
            ch_with_data = frozenset(
                ch for col, ch in channel_cols
                if float(row.get(col) or 0) > 0
            )
            candidates = channel_set_to_plexes.get(ch_with_data, [])
            if len(candidates) == 1:
                mixture = candidates[0]
            else:
                return None
        run_id_to_mixture[run_id] = mixture
        return mixture

    out_path = os.path.join(args.outdir, "msstats_input.tsv")
    run_id_to_mixture = {}
    rows_out = []
    runs_seen = set()

    with open(args.tsv, newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            path = (row.get("RunLoadPath") or "").strip()
            if not path:
                continue
            run_id = run_id_from_path(path)
            mixture = resolve_mixture(run_id, row)
            if mixture is None:
                raise SystemExit(
                    f"Run '{run_id}' could not be assigned to a Mixture. "
                    "Ensure run names contain Prot_NN matching Protein 10-Plex ID."
                )
            runs_seen.add(run_id)
            protein = (row.get("ProteinId") or "").strip()
            peptide = (row.get("PeptideSequence") or "").strip()
            charge = (row.get("Charge") or "").strip()
            if not protein or not peptide:
                continue
            psm = f"{peptide}_{charge}"
            for col, ch in channel_cols:
                try:
                    val = float(row.get(col) or 0)
                except (ValueError, TypeError):
                    continue
                if val <= 0:
                    continue
                key = (mixture, ch)
                if key not in plex_ch_ann:
                    continue
                bio, cond = plex_ch_ann[key]
                rows_out.append({
                    "ProteinName": protein,
                    "PeptideSequence": peptide,
                    "Charge": charge,
                    "PSM": psm,
                    "Mixture": mixture,
                    "TechRepMixture": 1,
                    "Run": run_id,
                    "Channel": ch,
                    "Condition": cond,
                    "BioReplicate": bio,
                    "Intensity": val,
                })

    # Fix counters: we incremented per row; report unique runs with/without fraction
    runs_fraction_ok = sum(1 for rid in runs_seen if parse_fraction_from_run_id(rid) is not None)
    runs_fraction_bad = len(runs_seen) - runs_fraction_ok

    print("Detected runs:", len(runs_seen))
    print("Runs with _F## fraction tag:", runs_fraction_ok, "; without:", runs_fraction_bad)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", "TechRepMixture",
                "Run", "Channel", "Condition", "BioReplicate", "Intensity",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows_out)
    print("Wrote", out_path, "(", len(rows_out), "rows)")

    annotation_rows = []
    for run_id, mix in run_id_to_mixture.items():
        frac = parse_fraction_from_run_id(run_id)
        frac_out = frac if frac is not None else 1
        for (p, ch), (bio, cond) in plex_ch_ann.items():
            if p != mix:
                continue
            annotation_rows.append({
                "Run": run_id,
                "Channel": ch,
                "Condition": cond,
                "BioReplicate": bio,
                "Mixture": mix,
                "Fraction": frac_out,
                "TechRepMixture": 1,
            })

    validate_annotation(
        annotation_rows, plex_ch_ann, run_id_to_mixture,
        args.allow_inconsistent_channels, args.allow_multiple_norm,
    )

    ann_path = os.path.join(args.outdir, "annotation_filled.csv")
    with open(ann_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["Run", "Channel", "Condition", "BioReplicate", "Mixture", "Fraction", "TechRepMixture"],
        )
        w.writeheader()
        w.writerows(annotation_rows)
    print("Wrote", ann_path, "(", len(annotation_rows), "rows)")

    all_ten = all(len(channels_per_plex[p]) == len(full_channel_set) for p in plex_ids)
    summary_lines = [
        "CCLE converter summary (bridge-corrected)",
        f"Real plexes (Protein 10-Plex ID >= 1): {n_plexes}",
        f"Bridge metadata rows (ID=0) in sample sheet: {len(bridge_meta)} (not emitted as a mixture)",
        f"Inferred bridge annotation rows (one per plex): {len(bridge_labels)}",
        f"Reporter channels in TSV: {len(full_channel_set)}",
        f"All real plexes have full channel annotation: {all_ten}",
        f"Unique MS runs: {len(runs_seen)}",
        f"Runs with parsed Fraction from _F##: {runs_fraction_ok}",
        f"Runs without _F## (Fraction defaulted to 1 in CSV): {runs_fraction_bad}",
        f"msstats_input long rows: {len(rows_out)}",
        f"annotation_filled rows: {len(annotation_rows)}",
        "Inconsistencies / warnings:",
    ]
    if inconsistencies:
        summary_lines.extend("  - " + x for x in inconsistencies)
    else:
        summary_lines.append("  (none)")
    summary_lines.append("Inferred bridge channel per plex (from set difference):")
    for p in sorted(bridge_labels.keys()):
        summary_lines.append(f"  plex {p}: {bridge_labels[p]}")
    summary_path = os.path.join(args.outdir, args.summary_file)
    with open(summary_path, "w", encoding="utf-8") as sf:
        sf.write("\n".join(summary_lines) + "\n")
    print("Wrote", summary_path)
    print("Next: run R pipeline with --msstats_input_dir", os.path.abspath(args.outdir))


if __name__ == "__main__":
    main()
