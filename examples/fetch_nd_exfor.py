#!/usr/bin/env python3
"""
Purpose:
  Fetch and parse nucleon-deuteron elastic scattering experimental data from
  the IAEA-NDS EXFOR archive (https://nds.iaea.org/nrdc/exfor-master/).

Design choice:
  EXFOR's web search UI is JavaScript-heavy and unstable for headless use.
  The raw per-entry text files at nds.iaea.org/nrdc/exfor-master/entry/<N>/<ENTRY>.txt
  are stable and directly parseable. This script targets that path.

What this does:
  - Given one or more EXFOR entry IDs (e.g. 22151), download the raw text.
  - Split into subentries (SUBENT ... ENDSUBENT).
  - Identify subentries whose REACTION line matches a user-provided pattern
    (default: "1-H-2(N,EL)" for n+d elastic) and whose COMMON/DATA columns
    contain angle + cross-section (or analyzing power).
  - Parse the DATA block into a tab-separated table and write to disk.

Usage:
  python3 examples/fetch_nd_exfor.py \
    --entries 22151,22172,22176 \
    --reaction "1-H-2(N,EL)" \
    --out-dir data/nd_scattering/

Known useful entries (n+d elastic, cross section or analyzing power):
  — Curate your own list by scanning the EXFOR catalog at
    https://nds.iaea.org/exfor/ (reaction = 1-H-2(N,EL), quantity = DA or POL).
  — Shimizu et al., Nucl. Phys. A 382, 242 (1982) — 65 MeV nd elastic
  — Rühl et al., Nucl. Phys. A 524, 377 (1991) — 67 MeV nd analyzing power
  — Ermisch et al., PRC 68, 051001 (2003) — 135-190 MeV nd elastic

This script does not attempt to be a general EXFOR parser; it is tailored to
the common DATA layout used for DA-type angular distributions.
"""

from __future__ import annotations

import argparse
import json
import re
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ARCHIVE_TEMPLATE = "https://nds.iaea.org/nrdc/exfor-master/entry/{prefix}/{entry}.txt"


@dataclass
class SubEntry:
    entry_id: str
    subent_id: str
    reaction: str = ""
    title: str = ""
    authors: str = ""
    reference: str = ""
    fields: List[str] = field(default_factory=list)          # column names
    units: List[str] = field(default_factory=list)           # column units
    common: Dict[str, Tuple[float, str]] = field(default_factory=dict)
    data_rows: List[List[float]] = field(default_factory=list)


def fetch_entry(entry_id: str, cache_dir: Optional[Path] = None) -> str:
    entry_id = entry_id.strip()
    if len(entry_id) < 4 or not re.match(r"^[0-9A-Za-z]{5}$", entry_id):
        raise ValueError(f"entry id must be a 5-char alphanumeric EXFOR ID (e.g. 22151, O1563); got {entry_id!r}")
    # EXFOR master archive uses the first character as directory prefix (case-sensitive).
    prefix = entry_id[0].lower()
    entry_lower = entry_id.lower()
    url = ARCHIVE_TEMPLATE.format(prefix=prefix, entry=entry_lower)

    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / f"{entry_id}.txt"
        if cache_path.exists():
            return cache_path.read_text(encoding="utf-8", errors="replace")

    req = urllib.request.Request(url, headers={"User-Agent": "Tic-tac-fetch/0.1"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    if "Not Found" in raw[:200] or "Access forbidden" in raw[:200]:
        raise RuntimeError(f"EXFOR entry {entry_id} not found at {url}")
    if cache_dir:
        (cache_dir / f"{entry_id}.txt").write_text(raw, encoding="utf-8")
    return raw


def _strip_exfor_cols(line: str) -> str:
    # EXFOR lines are 66 columns wide; after col 66 there is a sequence id
    # that we don't want to treat as data.
    return line[:66] if len(line) > 66 else line


def _parse_float_tokens(line: str) -> List[float]:
    # Tokens are fixed-width 11 chars, right-aligned. Blank -> NaN.
    body = _strip_exfor_cols(line)
    # Strip the first column (which is the line-type indicator).
    body = body[11:] if len(body) > 11 else ""
    values: List[float] = []
    for i in range(0, len(body), 11):
        tok = body[i:i+11].strip()
        if not tok:
            values.append(float("nan"))
            continue
        # EXFOR occasionally uses formats like "1.234+3" meaning 1.234e+3
        if re.match(r"^[\-+]?\d*\.\d+[+\-]\d+$", tok) or re.match(r"^[\-+]?\d+[+\-]\d+$", tok):
            # Insert E before sign
            m = re.match(r"^([+\-]?\d*\.?\d+)([+\-]\d+)$", tok)
            if m:
                tok = f"{m.group(1)}E{m.group(2)}"
        try:
            values.append(float(tok))
        except ValueError:
            values.append(float("nan"))
    return values


def _parse_header_text(line: str, prefix_len: int = 11) -> List[str]:
    body = _strip_exfor_cols(line)
    body = body[prefix_len:] if len(body) > prefix_len else ""
    tokens: List[str] = []
    for i in range(0, len(body), 11):
        tokens.append(body[i:i+11].strip())
    return tokens


def parse_subentries(raw: str, entry_id: str) -> List[SubEntry]:
    lines = raw.splitlines()
    subentries: List[SubEntry] = []
    current: Optional[SubEntry] = None
    in_bib = False
    in_common = False
    in_data = False
    last_bib_key: Optional[str] = None
    bib_buffer: Dict[str, List[str]] = {}

    common_fields: List[str] = []
    common_units: List[str] = []
    common_values: List[float] = []

    data_fields: List[str] = []
    data_units: List[str] = []

    for raw_line in lines:
        line = raw_line.rstrip("\r\n")
        stripped = line.strip()
        head6 = line[:11].strip()

        if head6 == "SUBENT" and stripped.startswith("SUBENT"):
            current = SubEntry(entry_id=entry_id, subent_id=stripped.split()[1])
            subentries.append(current)
            in_bib = False
            in_common = False
            in_data = False
            bib_buffer = {}
            continue

        if head6 == "ENDSUBENT":
            if current is not None:
                # Flush bib_buffer to structured fields on current
                current.title = " ".join(bib_buffer.get("TITLE", []))
                current.authors = " ".join(bib_buffer.get("AUTHOR", []))
                current.reference = " ".join(bib_buffer.get("REFERENCE", []))
                current.reaction = " ".join(bib_buffer.get("REACTION", []))
            current = None
            continue

        if current is None:
            continue

        if head6 == "BIB":
            in_bib = True
            in_common = False
            in_data = False
            last_bib_key = None
            continue
        if head6 == "ENDBIB":
            in_bib = False
            last_bib_key = None
            continue

        if head6 == "COMMON":
            in_common = True
            in_data = False
            common_fields = []
            common_units = []
            common_values = []
            continue
        if head6 == "ENDCOMMON":
            in_common = False
            # pack into dict
            for name, unit, val in zip(common_fields, common_units, common_values):
                current.common[name] = (val, unit)
            continue

        if head6 == "DATA":
            in_data = True
            in_common = False
            continue
        if head6 == "ENDDATA":
            in_data = False
            current.fields = list(data_fields)
            current.units = list(data_units)
            continue

        if in_bib:
            key = line[:11].strip()
            if key:
                last_bib_key = key
                bib_buffer.setdefault(last_bib_key, []).append(line[11:66].strip())
            elif last_bib_key is not None:
                bib_buffer[last_bib_key].append(line[11:66].strip())
            continue

        if in_common:
            # EXFOR COMMON layout: line 1 header (field names), line 2 units, line 3 values.
            # Unlike BIB keys, COMMON/DATA headers do NOT have a line-type prefix —
            # column 0 is the first field name.
            if not common_fields:
                common_fields = _parse_header_text(line, prefix_len=0)
                continue
            if not common_units:
                common_units = _parse_header_text(line, prefix_len=0)
                continue
            common_values = _parse_float_tokens("XXXXXXXXXXX" + line)
            continue

        if in_data:
            if not data_fields:
                data_fields = _parse_header_text(line, prefix_len=0)
                continue
            if not data_units:
                data_units = _parse_header_text(line, prefix_len=0)
                continue
            row = _parse_float_tokens("XXXXXXXXXXX" + line)
            if any(not (v != v) for v in row):  # at least one non-NaN
                current.data_rows.append(row)
            continue

    # Propagate common-BIB (TITLE / AUTHORS / REFERENCE) from SUBENT *001 to the
    # sibling subentries, which usually only declare REACTION/METHOD/etc locally.
    shared: Optional[SubEntry] = None
    for sub in subentries:
        if sub.subent_id.endswith("001"):
            shared = sub
            break
    if shared is not None:
        for sub in subentries:
            if sub is shared:
                continue
            if not sub.title:
                sub.title = shared.title
            if not sub.authors:
                sub.authors = shared.authors
            if not sub.reference:
                sub.reference = shared.reference

    return subentries


def _normalize_reaction(s: str) -> str:
    return re.sub(r"\s+", "", s.upper())


def filter_subentries(subentries: List[SubEntry], reaction_pattern: str) -> List[SubEntry]:
    pattern = _normalize_reaction(reaction_pattern)
    kept: List[SubEntry] = []
    for sub in subentries:
        react = _normalize_reaction(sub.reaction)
        if pattern in react:
            kept.append(sub)
    return kept


def format_output_table(sub: SubEntry) -> str:
    out_lines: List[str] = []
    out_lines.append(f"# EXFOR entry:   {sub.entry_id}")
    out_lines.append(f"# SUBENT:        {sub.subent_id}")
    out_lines.append(f"# TITLE:         {sub.title}")
    out_lines.append(f"# AUTHORS:       {sub.authors}")
    out_lines.append(f"# REFERENCE:     {sub.reference}")
    out_lines.append(f"# REACTION:      {sub.reaction}")
    if sub.common:
        out_lines.append("# COMMON:")
        for name, (val, unit) in sub.common.items():
            out_lines.append(f"#   {name}={val}  [{unit}]")
    out_lines.append("# DATA columns: " + "\t".join(f"{f}[{u}]" for f, u in zip(sub.fields, sub.units)))
    for row in sub.data_rows:
        out_lines.append("\t".join(f"{v:.6e}" if v == v else "NaN" for v in row))
    return "\n".join(out_lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Fetch and parse EXFOR entries for n+d elastic data")
    p.add_argument("--entries", required=True,
                   help="Comma-separated EXFOR entry IDs (e.g. '22151,22172')")
    p.add_argument("--reaction", default="1-H-2(N,EL)",
                   help="EXFOR REACTION substring to keep (case-insensitive, whitespace ignored)")
    p.add_argument("--out-dir", default="data/nd_scattering",
                   help="Directory to write per-subentry .txt tables")
    p.add_argument("--cache-dir", default="data/nd_scattering/_raw_cache",
                   help="Directory to cache raw EXFOR entry files")
    return p


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    out_dir = (root / args.out_dir).resolve()
    cache_dir = (root / args.cache_dir).resolve() if args.cache_dir else None
    out_dir.mkdir(parents=True, exist_ok=True)

    entry_ids = [e.strip() for e in args.entries.split(",") if e.strip()]
    manifest: List[Dict[str, object]] = []
    for eid in entry_ids:
        try:
            raw = fetch_entry(eid, cache_dir)
        except Exception as exc:
            print(f"[fetch_nd_exfor] {eid}: FETCH FAILED: {exc}")
            manifest.append({"entry": eid, "status": "fetch_failed", "error": str(exc)})
            continue
        subs = parse_subentries(raw, eid)
        matches = filter_subentries(subs, args.reaction)
        print(f"[fetch_nd_exfor] {eid}: {len(subs)} subentries, {len(matches)} match '{args.reaction}'")
        for sub in matches:
            if not sub.data_rows:
                continue
            fname = f"exfor_{sub.entry_id}_{sub.subent_id[-3:]}.tsv"
            fpath = out_dir / fname
            fpath.write_text(format_output_table(sub), encoding="utf-8")
            manifest.append({
                "entry": sub.entry_id,
                "subent": sub.subent_id,
                "reaction": sub.reaction,
                "reference": sub.reference,
                "rows": len(sub.data_rows),
                "path": str(fpath),
            })

    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[fetch_nd_exfor] manifest -> {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
