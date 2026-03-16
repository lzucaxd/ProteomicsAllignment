#!/usr/bin/env python3
"""
Convert CCLE peptide-level TSV + sample info (Sheet2) to MSstatsTMT input format.
Outputs msstats_input.tsv and annotation_filled.csv so the main R pipeline can
run with --msstats_input_dir (skips PSM parsing; CPTAC flow unchanged).

Usage (from data/ or repo root):
  python3 ccle_peptide/ccle_to_msstats_input.py --tsv ccle_peptide/ccle_protein_quant_with_peptides_14745.tsv --sample_csv ccle_peptide/sample_info_ccle.csv --outdir results/CCLE

Prerequisite: export sample sheet2 first:
  python3 ccle_peptide/export_sample_sheet2_csv.py
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

# Reporter-ion columns: rq_126_sn, rq_127n_sn, rq_127c_sn, ... (auto-detected; supports TMT6/10/11/16/TMTpro/18)
RQ_SN_PATTERN = re.compile(r"^rq_(\d+)([nc])?_sn$", re.IGNORECASE)


def channel_sort_key(label):
    """Sort key for TMT channels: 126, 127N, 127C, 128N, ... 131, 131N, 131C, ... 134N."""
    m = re.match(r"^(\d+)([NC])?$", str(label).upper())
    if not m:
        return (0, "Z")
    num = int(m.group(1))
    suffix = m.group(2) or ""
    # '' < 'N' < 'C' for same number
    suf_order = {"": 0, "N": 1, "C": 2}
    return (num, suf_order.get(suffix, 99))


def detect_channel_columns(fieldnames):
    """Detect rq_*_sn columns and return [(col, channel_label), ...] in canonical order.
    Channel label: 126, 127N, 127C, ... (number + optional N/C, capitalized).
    """
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
    """Sample file may have 127n/127c; MSstatsTMT expects 127N, 127C, etc. Works for any TMT plex."""
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
    """Load sample_info_ccle.csv (Sheet2). Returns list of dicts with keys:
    Cell Line, CCLE Code, Tissue of Origin, Protein 10-Plex ID, Protein TMT Label, Notes.
    """
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append({k.strip(): (v.strip() if v else "") for k, v in row.items()})
    return rows


def build_plex_channel_annotation(sample_rows):
    """(PlexID, Channel) -> (BioReplicate, Condition). PlexID = Mixture (Protein 10-Plex ID). Channel = MSstatsTMT canonical (127N, 127C, ...). Condition=Norm if Notes contains 'Bridge' else Sample."""
    key_col = "Protein 10-Plex ID"
    ch_col = "Protein TMT Label"
    cell_col = "Cell Line"
    notes_col = "Notes"
    ann = {}
    for row in sample_rows:
        try:
            plex = int(row.get(key_col, ""))
        except (ValueError, TypeError):
            continue
        ch_raw = (row.get(ch_col) or "").strip()
        ch = normalize_channel(ch_raw)
        if not ch:
            continue
        cell = (row.get(cell_col) or "").strip() or f"Plex{plex}_{ch}"
        notes = (row.get(notes_col) or "").lower()
        cond = "Norm" if "bridge" in notes else "Sample"
        ann[(plex, ch)] = (cell, cond)
    return ann


def run_id_from_path(path):
    """RunLoadPath -> run_id (basename without extension)."""
    base = os.path.basename(path)
    return base.replace(".mzXML", "").replace(".raw", "").replace(".mzML", "").strip()


# Run name pattern: Prot_NN (e.g. Prot_01, Prot_11) -> plex number. Sample file uses 0-based Protein 10-Plex ID;
# CCLE run names typically use 1-based (Prot_01 = first plex = 0). So mixture = int(NN) - 1 for NN in 1..43.
RUN_PLEX_PATTERN = re.compile(r"Prot_(\d+)", re.IGNORECASE)


def mixture_from_run_id(run_id, valid_plex_ids):
    """Infer Mixture (plex ID) from run_id. Uses Prot_NN in run name (1-based: Prot_01 -> 0). Returns int or None if not determinable."""
    m = RUN_PLEX_PATTERN.search(run_id)
    if not m:
        return None
    num = int(m.group(1))
    # 1-based run naming: Prot_01 = plex 0, Prot_02 = plex 1, ... Prot_43 = plex 42
    mixture = num - 1
    if mixture in valid_plex_ids:
        return mixture
    # 0-based: Prot_00 -> 0, Prot_01 -> 1
    if num in valid_plex_ids:
        return num
    return None


def build_plex_metadata(sample_rows):
    """From sample file build: plex_ch_ann (Plex,Channel)->(BioReplicate,Condition), plex_ids set, channels_per_plex dict, norm_count_per_plex."""
    plex_ch_ann = build_plex_channel_annotation(sample_rows)
    plex_ids = {p for (p, _) in plex_ch_ann}
    channels_per_plex = {}
    norm_count_per_plex = {}
    for (p, ch), (bio, cond) in plex_ch_ann.items():
        channels_per_plex.setdefault(p, set()).add(ch)
        if cond == "Norm":
            norm_count_per_plex[p] = norm_count_per_plex.get(p, 0) + 1
    return plex_ch_ann, plex_ids, channels_per_plex, norm_count_per_plex


def validate_annotation(annotation_rows, plex_ch_ann, run_id_to_mixture, allow_inconsistent_channels=False, allow_multiple_norm=False):
    """Enforce MSstatsTMT requirements: same channel set per mixture, one Norm per mixture, Run->single Mixture, (Mixture,Channel) unique. Raises SystemExit on failure unless allow_* flags set."""
    if not annotation_rows:
        raise SystemExit("Validation failed: annotation is empty.")

    # Check: (Run, Channel) unique — no duplicate (Run, Channel) (annotation corruption; multiple runs per mixture are allowed)
    seen = set()
    for r in annotation_rows:
        key = (r["Run"], r["Channel"])
        if key in seen:
            raise SystemExit(
                f"Validation failed: (Run, Channel) must be unique. Duplicate: Run {key[0]}, Channel {key[1]}."
            )
        seen.add(key)

    # Group by Mixture for channel consistency and Norm count
    by_mix = {}
    for r in annotation_rows:
        m = r["Mixture"]
        if m not in by_mix:
            by_mix[m] = {"channels": set(), "norm_count": 0}
        by_mix[m]["channels"].add(r["Channel"])
        if (r.get("Condition") or "").strip().lower() == "norm":
            by_mix[m]["norm_count"] += 1

    # Check: channel consistency across mixtures — unique(channel_set_per_mixture) == 1 (would break MSstats normalization)
    ch_sets = [frozenset(v["channels"]) for v in by_mix.values()]
    if len(set(ch_sets)) != 1:
        if not allow_inconsistent_channels:
            raise SystemExit(
                "Validation failed: inconsistent TMT channel structure across mixtures. "
                "All plexes must use the same channel set (e.g. Mixture 0 and Mixture 1 must have the same channels). "
                "Use --allow-inconsistent-channels to run anyway (e.g. CCLE where plex 0 has 10 channels, others 9)."
            )
        print("Warning: Inconsistent TMT channel structure across mixtures (e.g. plex 0 has 10 channels, others 9). Proceeding with --allow-inconsistent-channels.")

    # Check: exactly one Norm per mixture (unless allow_multiple_norm, e.g. all-bridge plex)
    for m, v in by_mix.items():
        if v["norm_count"] != 1 and not allow_multiple_norm:
            raise SystemExit(
                f"Validation failed: Mixture {m} has {v['norm_count']} Norm channel(s). "
                "MSstatsTMT expects exactly one bridge (Condition=Norm) per mixture. "
                "Use --allow-multiple-norm to run anyway (e.g. CCLE all-bridge reference plex)."
            )
    if allow_multiple_norm and any(by_mix[m]["norm_count"] != 1 for m in by_mix):
        print("Warning: At least one mixture has multiple Norm channels (e.g. all-bridge reference plex). Proceeding with --allow-multiple-norm.")

    # Check 3: Run must map to exactly one Mixture
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
    ap = argparse.ArgumentParser(description="CCLE peptide TSV + sample sheet -> MSstatsTMT input")
    ap.add_argument("--tsv", required=True, help="CCLE peptide TSV (e.g. ccle_protein_quant_with_peptides_14745.tsv)")
    ap.add_argument("--sample_csv", required=True, help="Sample info CSV from Sheet2 (sample_info_ccle.csv)")
    ap.add_argument("--outdir", default="results/CCLE", help="Output directory (msstats_input.tsv, annotation_filled.csv)")
    ap.add_argument("--allow-inconsistent-channels", action="store_true", help="Allow different channel sets per mixture (e.g. CCLE plex 0 has 10 ch, others 9)")
    ap.add_argument("--allow-multiple-norm", action="store_true", help="Allow more than one Norm per mixture (e.g. all-bridge reference plex)")
    args = ap.parse_args()

    if not os.path.isfile(args.tsv):
        raise SystemExit(f"Not found: {args.tsv}")
    if not os.path.isfile(args.sample_csv):
        raise SystemExit(f"Not found: {args.sample_csv}")

    sample_rows = load_sample_info(args.sample_csv)
    plex_ch_ann, plex_ids, channels_per_plex, norm_count_per_plex = build_plex_metadata(sample_rows)
    if not plex_ids:
        raise SystemExit("Sample file has no valid plex rows (Protein 10-Plex ID + Protein TMT Label).")

    # Channel-set -> list of plex IDs (for fallback when run name has no Prot_NN)
    channel_set_to_plexes = defaultdict(list)
    for p in plex_ids:
        channel_set_to_plexes[frozenset(channels_per_plex[p])].append(p)

    os.makedirs(args.outdir, exist_ok=True)

    # Auto-detect channel columns (rq_*_sn) — supports TMT6, TMT10, TMT11, TMT16/TMTpro, TMT18
    with open(args.tsv, newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = list(r.fieldnames or [])
    channel_cols = detect_channel_columns(fieldnames)
    if not channel_cols:
        raise SystemExit(
            "No reporter-ion columns found. Expected TSV columns matching rq_*_sn (e.g. rq_126_sn, rq_127n_sn, rq_127c_sn)."
        )
    channel_col_names = [c for c, _ in channel_cols]
    required = ["ProteinId", "PeptideSequence", "Charge", "RunLoadPath"] + channel_col_names
    for col in required:
        if col not in fieldnames:
            raise SystemExit(f"TSV missing column: {col}")

    # Logging: plex and channel summary from sample file
    n_plexes = len(plex_ids)
    size_to_plexes = defaultdict(list)
    for p in plex_ids:
        size_to_plexes[len(channels_per_plex[p])].append(p)
    n_bridge = sum(1 for p in plex_ids if norm_count_per_plex.get(p, 0) >= 1)
    print("Detected plexes:", n_plexes)
    if len(size_to_plexes) == 1:
        (n_ch,) = size_to_plexes.keys()
        print("Channels per plex:", n_ch)
    else:
        main_size, main_plexes = max(size_to_plexes.items(), key=lambda x: len(x[1]))
        others = [(s, plist) for s, plist in size_to_plexes.items() if s != main_size]
        except_str = ", ".join(f"plex {p}: {s}" for s, plist in others for p in plist)
        print("Channels per plex:", main_size, "(except", except_str + ")")
    print("Plexes with Norm (bridge) channel:", n_bridge)
    print("Detected", len(channel_cols), "channels in TSV:", [ch for _, ch in channel_cols])

    def resolve_mixture(run_id, row):
        """Resolve Mixture for this run: from run name (Prot_NN) or from channel structure in row. Cached per run."""
        if run_id in run_id_to_mixture:
            return run_id_to_mixture[run_id]
        mixture = mixture_from_run_id(run_id, plex_ids)
        if mixture is None:
            # Fallback: which channels have data in this row?
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
                    "Ensure run names contain plex id (e.g. Prot_01, Prot_02) or that channel structure uniquely identifies one plex."
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

    print("Detected runs:", len(runs_seen))

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", "TechRepMixture", "Run", "Channel", "Condition", "BioReplicate", "Intensity"], delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)
    print("Wrote", out_path, "(", len(rows_out), "rows)")

    # Annotation: one row per (Run, Channel). Mixture = Protein 10-Plex ID from sample file; Run = MS run (fraction).
    annotation_rows = []
    for run_id, mix in run_id_to_mixture.items():
        for (p, ch), (bio, cond) in plex_ch_ann.items():
            if p != mix:
                continue
            annotation_rows.append({
                "Run": run_id,
                "Channel": ch,
                "Condition": cond,
                "BioReplicate": bio,
                "Mixture": mix,
                "Fraction": 1,
                "TechRepMixture": 1,
            })

    validate_annotation(annotation_rows, plex_ch_ann, run_id_to_mixture, args.allow_inconsistent_channels, args.allow_multiple_norm)

    ann_path = os.path.join(args.outdir, "annotation_filled.csv")
    with open(ann_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Run", "Channel", "Condition", "BioReplicate", "Mixture", "Fraction", "TechRepMixture"])
        w.writeheader()
        w.writerows(annotation_rows)
    print("Wrote", ann_path, "(", len(annotation_rows), "rows)")
    print("Next: run R pipeline with --msstats_input_dir", os.path.abspath(args.outdir))


if __name__ == "__main__":
    main()
