#!/usr/bin/env python3
"""
After downloading: check that each study (from manifests or pdc_psm/) has a
corresponding row in sample_files_msstats_tmt.csv.

Usage:
  python check_studies_sample_file.py
  python check_studies_sample_file.py --manifests manifests/
  python check_studies_sample_file.py --psm_dir pdc_psm
"""
import csv
import os
import sys
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent
SAMPLE_CSV = DATA_DIR / "sample_files_msstats_tmt.csv"
MANIFESTS_DIR = DATA_DIR / "manifests"
PSM_ROOT = DATA_DIR / "pdc_psm"


def get_study_ids_from_manifests(manifests_dir: Path) -> set:
    studies = set()
    for path in sorted(manifests_dir.glob("PDC_file_manifest_*.csv")):
        with open(path, newline="", encoding="utf-8", errors="replace") as f:
            r = csv.DictReader(f)
            for row in r:
                sid = (row.get("PDC Study ID") or "").strip().strip('"')
                if sid:
                    studies.add(sid)
                    break
    return studies


def get_study_ids_from_psm_dir(psm_root: Path) -> set:
    if not psm_root.is_dir():
        return set()
    return {d.name for d in psm_root.iterdir() if d.is_dir() and d.name.startswith("PDC")}


def load_sample_csv(sample_path: Path) -> dict:
    """Return dict study_id -> row (path, file_name, etc.)."""
    out = {}
    if not sample_path.is_file():
        return out
    with open(sample_path, newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f)
        for row in r:
            sid = (row.get("study_id") or "").strip().strip('"')
            if sid:
                out[sid] = row
    return out


def main():
    import argparse
    p = argparse.ArgumentParser(description="Check studies have entry in sample_files_msstats_tmt.csv")
    p.add_argument("--manifests", default=str(MANIFESTS_DIR), help="Manifests directory")
    p.add_argument("--psm_dir", default=str(PSM_ROOT), help="pdc_psm root (study subdirs)")
    p.add_argument("--sample_csv", default=str(SAMPLE_CSV), help="sample_files_msstats_tmt.csv path")
    p.add_argument("--study", default=None, help="Check only this study_id; exit 0 if present, 1 if missing")
    args = p.parse_args()

    manifests_dir = Path(args.manifests)
    psm_root = Path(args.psm_dir)
    sample_path = Path(args.sample_csv)
    only_study = (args.study or "").strip()

    studies_manifests = get_study_ids_from_manifests(manifests_dir)
    studies_psm = get_study_ids_from_psm_dir(psm_root)
    studies = studies_manifests or studies_psm
    if only_study:
        studies = {only_study}
    elif not studies:
        print("No studies found (check --manifests and --psm_dir).")
        sys.exit(0)

    sample_rows = load_sample_csv(sample_path)
    if not sample_rows and not only_study:
        print(f"Warning: no rows in {sample_path}")

    def resolve_sample_path(path_str: str) -> Path:
        raw = path_str.strip().strip('"').strip("'")
        if not raw:
            return Path("")
        p = Path(raw)
        if p.is_file():
            return p
        rel_data = DATA_DIR / raw
        if rel_data.is_file():
            return rel_data
        mirror = (os.environ.get("CPTAC_LOCAL_MIRROR") or "").strip()
        if mirror:
            under_mirror = Path(mirror) / raw
            if under_mirror.is_file():
                return under_mirror
        return rel_data

    have = []
    missing = []
    for sid in sorted(studies):
        row = sample_rows.get(sid)
        if row:
            path = (row.get("path") or "").strip()
            file_name = (row.get("file_name") or "").strip()
            path_obj = resolve_sample_path(path) if path else Path("")
            path_exists = "yes" if path and path_obj.is_file() else "no (path missing)"
            have.append((sid, path_exists, file_name))
        else:
            missing.append(sid)

    if only_study:
        if missing:
            print(f"  {only_study}: no entry in sample_files_msstats_tmt.csv")
            sys.exit(1)
        print(f"  {only_study}: has entry (path_exists={have[0][1]}, file={have[0][2]})")
        sys.exit(0)

    print("Studies vs sample_files_msstats_tmt.csv")
    print("=" * 60)
    print(f"\nHave sample file entry: {len(have)}")
    for sid, path_ok, fname in have:
        print(f"  {sid}  path_exists={path_ok}  file={fname}")

    print(f"\nMissing sample file entry: {len(missing)}")
    for sid in missing:
        print(f"  {sid}")

    if missing:
        print("\n→ Add rows for the missing studies to sample_files_msstats_tmt.csv (study_id, path, file_name, format, reference_channel, etc.).")
    sys.exit(1 if missing else 0)


if __name__ == "__main__":
    main()
