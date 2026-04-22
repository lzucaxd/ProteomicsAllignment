#!/usr/bin/env python3
"""Check manifests folder: list studies per file and find duplicate File IDs / File Names."""
import csv
from collections import defaultdict
from pathlib import Path

MANIFESTS_DIR = Path(__file__).resolve().parent / "manifests"

def norm(s):
    return (s or "").strip().strip('"')

def main():
    manifests = sorted(MANIFESTS_DIR.glob("*.csv"))
    if not manifests:
        print("No CSV files in manifests/")
        return

    # Per-manifest: studies, file_ids, file_names
    all_file_ids = defaultdict(list)   # file_id -> [(manifest, row)]
    all_file_names = defaultdict(list) # file_name -> [(manifest, row)]
    studies_per_manifest = {}
    study_to_manifests = defaultdict(list)

    for path in manifests:
        if path.name == "example_pdc_file_manifest.csv":
            continue
        studies = set()
        with open(path, newline="", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            if not rows:
                studies_per_manifest[path.name] = []
                continue
            for row in rows:
                study = norm(row.get("PDC Study ID", ""))
                if study:
                    studies.add(study)
                    study_to_manifests[study].append(path.name)
                fid = norm(row.get("File ID", ""))
                fname = norm(row.get("File Name", ""))
                if fid:
                    all_file_ids[fid].append((path.name, row))
                if fname:
                    all_file_names[fname].append((path.name, row))
            studies_per_manifest[path.name] = sorted(studies)

    # Report studies per manifest
    print("=" * 60)
    print("MANIFESTS FOLDER: studies per file")
    print("=" * 60)
    for name in sorted(studies_per_manifest.keys()):
        studies = studies_per_manifest[name]
        print(f"\n{name}")
        print(f"  Studies ({len(studies)}): {', '.join(studies)}")

    # All unique studies
    all_studies = sorted(set().union(*[set(s) for s in studies_per_manifest.values()]))
    print("\n" + "=" * 60)
    print("ALL UNIQUE STUDIES (across all manifests)")
    print("=" * 60)
    print(", ".join(all_studies))
    print(f"\nTotal: {len(all_studies)} studies")

    # Duplicates: File ID or File Name appearing in more than one manifest
    dup_ids = {k: v for k, v in all_file_ids.items() if len(v) > 1}
    dup_names = {k: v for k, v in all_file_names.items() if len(v) > 1}

    print("\n" + "=" * 60)
    print("DUPLICATES")
    print("=" * 60)
    print(f"\nFile IDs appearing in more than one manifest: {len(dup_ids)}")
    if dup_ids:
        for fid, occurrences in list(dup_ids.items())[:10]:
            manifests_seen = [m for m, _ in occurrences]
            print(f"  {fid[:40]}... in {manifests_seen}")
        if len(dup_ids) > 10:
            print(f"  ... and {len(dup_ids) - 10} more")

    print(f"\nFile names appearing in more than one manifest: {len(dup_names)}")
    if dup_names:
        for fname, occurrences in list(dup_names.items())[:10]:
            manifests_seen = [m for m, _ in occurrences]
            print(f"  {fname} in {manifests_seen}")
        if len(dup_names) > 10:
            print(f"  ... and {len(dup_names) - 10} more")

    # Studies that appear in multiple manifests (unique manifest filenames per study)
    study_to_unique_manifests = {s: sorted(set(study_to_manifests[s])) for s in all_studies}
    multi_manifest_studies = {s: study_to_unique_manifests[s] for s in all_studies if len(study_to_unique_manifests[s]) > 1}
    print(f"\nStudies that appear in more than one manifest file: {len(multi_manifest_studies)}")
    for s, mans in sorted(multi_manifest_studies.items()):
        print(f"  {s}: {mans}")

if __name__ == "__main__":
    main()
