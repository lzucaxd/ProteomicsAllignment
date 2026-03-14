#!/usr/bin/env python3
"""
Generic PDC manifest downloader (works with PDC/PDC-like CSV manifests).

Features:
- Auto-detect columns for: filename, download URL, size, md5, category
- Filter by regex / contains / extension
- Polite, rate-limit aware downloads (429 Retry-After)
- Retries with exponential backoff for 5xx/network errors
- Resume via Range requests + .part files
- Optional MD5 verification
- Skips files already downloaded with correct size (and md5 if enabled)

Usage examples:
  # Download PSM files only (.psm) from a study ID
  python pdc_manifest_downloader.py --manifest manifest.csv --outdir ./dl \
    --include-category "Peptide Spectral Matches" --ext .psm

  # Download peptide summary files (.peptides.tsv) for anything peptide-related
  python pdc_manifest_downloader.py --manifest manifest.csv --outdir ./dl \
    --include-any "peptide" --ext .tsv --include-name "peptides"

  # Dry-run to see how many files match
  python pdc_manifest_downloader.py --manifest manifest.csv --outdir ./dl \
    --include-any "peptide" --dry-run
"""

import argparse
import csv
import hashlib
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import requests


# ---------- Column detection helpers ----------

CANDIDATE_URL_KEYS = [
    "File Download Link", "Download Link", "Download URL", "URL", "Url", "link", "download_url"
]
CANDIDATE_NAME_KEYS = [
    "File Name", "Filename", "File", "Name", "file_name"
]
CANDIDATE_SIZE_KEYS = [
    "File Size (in bytes)", "Size (bytes)", "Size", "Bytes", "file_size"
]
CANDIDATE_MD5_KEYS = [
    "Md5sum", "MD5", "md5", "Checksum", "checksum"
]
CANDIDATE_CATEGORY_KEYS = [
    "Data Category", "Category", "File Category", "data_category"
]
CANDIDATE_STUDY_ID_KEYS = [
    "PDC Study ID", "Study ID", "Study", "pdc_study_id"
]
CANDIDATE_DOWNLOADABLE_KEYS = [
    "Downloadable", "Is Downloadable", "downloadable"
]


def pick_first_existing(header: List[str], candidates: List[str]) -> Optional[str]:
    header_set = {h.strip(): h for h in header}
    for c in candidates:
        if c in header_set:
            return header_set[c]
    # also try case-insensitive match
    lower = {h.strip().lower(): h for h in header}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


# ---------- Download logic ----------

@dataclass
class FileRow:
    name: str
    url: str
    size: Optional[int] = None
    md5: Optional[str] = None
    category: Optional[str] = None
    study_id: Optional[str] = None


def md5_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def safe_int(x: Optional[str]) -> Optional[int]:
    if x is None:
        return None
    x = str(x).strip()
    if not x:
        return None
    try:
        return int(float(x))
    except ValueError:
        return None


def should_skip_existing(out_path: Path, expected_size: Optional[int], expected_md5: Optional[str]) -> bool:
    if not out_path.exists():
        return False
    if expected_size is not None and out_path.stat().st_size != expected_size:
        return False
    if expected_md5:
        got = md5_file(out_path).lower()
        return got == expected_md5.lower()
    return True


def download_with_resume(
    session: requests.Session,
    url: str,
    out_path: Path,
    expected_size: Optional[int],
    expected_md5: Optional[str],
    *,
    polite_delay: float = 0.35,
    timeout: int = 60,
    max_retries: int = 10,
    base_backoff: float = 1.0,
    max_backoff: float = 90.0,
) -> bool:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if should_skip_existing(out_path, expected_size, expected_md5):
        return True

    part_path = out_path.with_suffix(out_path.suffix + ".part")
    existing = part_path.stat().st_size if part_path.exists() else 0

    headers: Dict[str, str] = {}
    if existing > 0:
        headers["Range"] = f"bytes={existing}-"

    for attempt in range(max_retries):
        time.sleep(polite_delay)
        try:
            r = session.get(url, stream=True, headers=headers, timeout=timeout)
        except requests.RequestException as e:
            wait = min(max_backoff, base_backoff * (2 ** attempt))
            print(f"[retry] network error: {e} | wait {wait:.1f}s")
            time.sleep(wait)
            continue

        if r.status_code == 429 or 500 <= r.status_code < 600:
            retry_after = r.headers.get("Retry-After")
            if retry_after:
                try:
                    wait = float(retry_after)
                except ValueError:
                    wait = min(max_backoff, base_backoff * (2 ** attempt))
            else:
                wait = min(max_backoff, base_backoff * (2 ** attempt))
            print(f"[retry] HTTP {r.status_code} | wait {wait:.1f}s")
            time.sleep(wait)
            continue

        if r.status_code not in (200, 206):
            print(f"[fail] HTTP {r.status_code} for {out_path.name}")
            return False

        # If server ignored Range, restart
        if r.status_code == 200 and existing > 0:
            existing = 0
            headers.pop("Range", None)
            if part_path.exists():
                part_path.unlink()

        mode = "ab" if r.status_code == 206 else "wb"
        with part_path.open(mode) as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

        # If size known, ensure complete; else accept as-is
        if expected_size is not None:
            got_size = part_path.stat().st_size
            if got_size != expected_size:
                # Continue with Range on next attempt
                headers["Range"] = f"bytes={got_size}-"
                wait = min(max_backoff, base_backoff * (2 ** attempt))
                print(f"[retry] incomplete size {got_size}/{expected_size} | wait {wait:.1f}s")
                time.sleep(wait)
                continue

        part_path.replace(out_path)

        if expected_md5:
            got = md5_file(out_path).lower()
            if got != expected_md5.lower():
                print(f"[retry] MD5 mismatch: {out_path.name}")
                out_path.unlink(missing_ok=True)
                # start over
                existing = 0
                headers.pop("Range", None)
                wait = min(max_backoff, base_backoff * (2 ** attempt))
                time.sleep(wait)
                continue

        return True

    print(f"[fail] exceeded retries: {out_path.name}")
    return False


# ---------- Filtering ----------

def ci_contains(hay: Optional[str], needle: str) -> bool:
    return needle.lower() in (hay or "").lower()


def compile_optional_regex(pattern: Optional[str]) -> Optional[re.Pattern]:
    if not pattern:
        return None
    return re.compile(pattern)


def filter_rows(
    rows: List[FileRow],
    include_any: Optional[str],
    include_category: Optional[str],
    include_name: Optional[str],
    include_regex: Optional[re.Pattern],
    ext: Optional[str],
    study_id: Optional[str],
) -> List[FileRow]:
    out: List[FileRow] = []
    for r in rows:
        if study_id and (r.study_id != study_id):
            continue
        if include_category and not ci_contains(r.category, include_category):
            continue
        if include_name and not ci_contains(r.name, include_name):
            continue
        if include_any:
            blob = " | ".join([r.name or "", r.category or "", r.study_id or ""])
            if not ci_contains(blob, include_any):
                continue
        if include_regex and not include_regex.search(r.name):
            continue
        if ext and not r.name.endswith(ext):
            continue
        out.append(r)
    return out


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="Path to manifest CSV")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--study-id", default=None, help="Exact Study ID filter (if present in manifest)")
    ap.add_argument("--include-any", default=None, help='Case-insensitive contains across name/category/study (e.g., "peptide")')
    ap.add_argument("--include-category", default=None, help='Case-insensitive contains in category (e.g., "Peptide Spectral Matches")')
    ap.add_argument("--include-name", default=None, help='Case-insensitive contains in file name (e.g., "peptides")')
    ap.add_argument("--include-regex", default=None, help='Regex applied to file name (e.g., ".*plex_03.*")')
    ap.add_argument("--ext", default=None, help="File extension filter (e.g., .psm, .tsv, .gz)")
    ap.add_argument("--max-files", type=int, default=0, help="0 = no limit")
    ap.add_argument("--verify-md5", action="store_true")
    ap.add_argument("--polite-delay", type=float, default=0.35)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--flat", action="store_true", help="Store all files directly under outdir (no subfolders)")
    args = ap.parse_args()

    manifest = Path(args.manifest)
    outdir = Path(args.outdir)

    with manifest.open(newline="") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []
        url_key = pick_first_existing(header, CANDIDATE_URL_KEYS)
        name_key = pick_first_existing(header, CANDIDATE_NAME_KEYS)
        size_key = pick_first_existing(header, CANDIDATE_SIZE_KEYS)
        md5_key = pick_first_existing(header, CANDIDATE_MD5_KEYS)
        cat_key = pick_first_existing(header, CANDIDATE_CATEGORY_KEYS)
        study_key = pick_first_existing(header, CANDIDATE_STUDY_ID_KEYS)
        downloadable_key = pick_first_existing(header, CANDIDATE_DOWNLOADABLE_KEYS)

        if not url_key or not name_key:
            print("[error] Could not detect URL and/or filename columns.")
            print("Detected header columns:", header)
            sys.exit(2)

        rows: List[FileRow] = []
        for raw in reader:
            # If a "Downloadable" column exists, honor it
            if downloadable_key:
                if str(raw.get(downloadable_key, "")).strip().lower() not in ("yes", "true", "1", "open"):
                    continue

            url = (raw.get(url_key) or "").strip()
            name = (raw.get(name_key) or "").strip()
            if not url or not name:
                continue

            rows.append(
                FileRow(
                    name=name,
                    url=url,
                    size=safe_int(raw.get(size_key)) if size_key else None,
                    md5=(raw.get(md5_key) or "").strip() if md5_key else None,
                    category=(raw.get(cat_key) or "").strip() if cat_key else None,
                    study_id=(raw.get(study_key) or "").strip() if study_key else None,
                )
            )

    include_regex = compile_optional_regex(args.include_regex)
    filtered = filter_rows(
        rows,
        include_any=args.include_any,
        include_category=args.include_category,
        include_name=args.include_name,
        include_regex=include_regex,
        ext=args.ext,
        study_id=args.study_id,
    )

    print(f"Manifest rows loaded: {len(rows)}")
    print(f"Matched after filters: {len(filtered)}")

    if args.max_files and args.max_files > 0:
        filtered = filtered[: args.max_files]
        print(f"Limiting to first {len(filtered)} files.")

    if args.dry_run:
        # show a few examples
        for ex in filtered[:10]:
            print(f"  - {ex.study_id or 'NA'} | {ex.category or 'NA'} | {ex.name}")
        return

    session = requests.Session()
    session.headers.update({"User-Agent": "pdc-manifest-downloader/1.0"})

    ok = 0
    for i, r in enumerate(filtered, 1):
        # Output structure: outdir/{study_id}/{category_sanitized}/filename
        if args.flat:
            target = outdir / r.name
        else:
            study = r.study_id or "UNKNOWN_STUDY"
            cat = re.sub(r"[^A-Za-z0-9_.-]+", "_", (r.category or "UNCATEGORIZED"))[:120]
            target = outdir / study / cat / r.name

        expected_md5 = r.md5 if args.verify_md5 and r.md5 else None

        print(f"[{i}/{len(filtered)}] {r.name}")
        if download_with_resume(
            session,
            r.url,
            target,
            r.size,
            expected_md5,
            polite_delay=args.polite_delay,
        ):
            ok += 1

    print(f"Done: {ok}/{len(filtered)} succeeded.")


if __name__ == "__main__":
    main()
