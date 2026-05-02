#!/usr/bin/env python3
"""
migrate_legacy_cache.py — Walk legacy P123_sparse_JP_*.h5 files, reconstruct
the new P123Key from the filename, compute the canonical SHA-256 hash, and
(in --apply mode) copy/move them into the new hash-keyed cache layout with
proper root attributes and manifest entries.

Default mode: dry-run (no files are written).
Pass --apply to execute writes.

Usage:
    python3 migrate_legacy_cache.py \\
        --legacy-roots CPP/Output [CPP/Output2 ...] \\
        --target-root /path/to/cache \\
        [--apply] [--copy]

Options:
    --legacy-roots  One or more root directories to scan recursively.
    --target-root   Destination cache root (will be created if needed).
    --apply         Execute moves/copies; without this flag, dry-run only.
    --copy          Copy files instead of moving them (default: move).
    --verbose       Print extra debug info.
"""

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import sys
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# h5py is the only non-stdlib dependency.
# ---------------------------------------------------------------------------
try:
    import h5py
except ImportError:
    print("ERROR: h5py is not installed.  Run:  pip install h5py", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Regex that matches the legacy filename format.
# ---------------------------------------------------------------------------
LEGACY_P123_RE = re.compile(
    r"^P123_sparse_JP_(?P<two_J_3N>-?\d+)_(?P<P_3N>-?\d+)"
    r"_Np_(?P<Np>\d+)_Nq_(?P<Nq>\d+)_J2max_(?P<J2max>\d+)\.h5$"
)

# ---------------------------------------------------------------------------
# Historic defaults that were hard-wired in the legacy code (pre-cache era).
# ---------------------------------------------------------------------------
LEGACY_DEFAULTS = dict(
    Nphi=24,
    Nx=24,
    tensor_force=True,
    isospin_breaking_1S0=False,
    chebyshev_s=1.5,
    chebyshev_t=1.0,
)

P123_SCHEMA_VERSION = 1
WRITER_VERSION = "legacy_imported_v1"


# ---------------------------------------------------------------------------
# Canonical JSON — must mirror cache_keys.cpp: canonical_json(const P123Key&)
# Keys are in lexicographic order (same as sorted()).
# ---------------------------------------------------------------------------
def _quantize_1e9(x: float) -> float:
    """Quantize to 1e-9 precision, matching C++ floor(x*1e9+0.5)/1e9."""
    return math.floor(x * 1e9 + 0.5) / 1e9


def _fmt_double_q(x: float) -> str:
    """%.9f formatting of quantized value, matching C++ fmt_double_q."""
    return "%.9f" % _quantize_1e9(x)


def canonical_p123_json(
    schema_version: int,
    Np_WP: int,
    Nq_WP: int,
    J_2N_max: int,
    two_J_3N: int,
    P_3N: int,
    Nphi: int,
    Nx: int,
    tensor_force: bool,
    isospin_breaking_1S0: bool,
    chebyshev_s: float,
    chebyshev_t: float,
) -> str:
    """
    Produce canonical JSON matching C++ canonical_json(P123Key).
    Keys are sorted lexicographically; doubles use %.9f after 1e-9 quantization;
    booleans are lowercase true/false.
    """
    # Manually assemble in lex-sorted key order (same as C++ explicit order):
    # J_2N_max, Nphi, Np_WP, Nq_WP, Nx, P_3N,
    # chebyshev_s, chebyshev_t, isospin_breaking_1S0,
    # schema_version, tensor_force, two_J_3N
    pairs = [
        ("J_2N_max", str(J_2N_max)),
        ("Nphi", str(Nphi)),
        ("Np_WP", str(Np_WP)),
        ("Nq_WP", str(Nq_WP)),
        ("Nx", str(Nx)),
        ("P_3N", str(P_3N)),
        ("chebyshev_s", _fmt_double_q(chebyshev_s)),
        ("chebyshev_t", _fmt_double_q(chebyshev_t)),
        ("isospin_breaking_1S0", "true" if isospin_breaking_1S0 else "false"),
        ("schema_version", str(schema_version)),
        ("tensor_force", "true" if tensor_force else "false"),
        ("two_J_3N", str(two_J_3N)),
    ]
    inner = ",".join(f'"{k}":{v}' for k, v in pairs)
    return "{" + inner + "}"


def sha256_hex(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def hash_short(full_hex: str) -> str:
    return full_hex[:6]


def filename_prefix(
    Np_WP: int, Nq_WP: int, J_2N_max: int, two_J_3N: int, P_3N: int
) -> str:
    """Human-readable prefix matching C++ filename_prefix(P123Key)."""
    parity = "+1" if P_3N > 0 else "-1"
    return f"Np{Np_WP}_Nq{Nq_WP}_J2max{J_2N_max}_JP{two_J_3N}{parity}"


def now_utc_iso8601() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# Root attribute writers
# ---------------------------------------------------------------------------
def write_root_attrs(h5file: h5py.File, key_hash_full: str, key_json: str) -> None:
    """Tag an HDF5 file with the new cache root attributes."""
    h5file.attrs["tictac_cache_kind"] = "p123"
    h5file.attrs["tictac_schema_version"] = P123_SCHEMA_VERSION
    h5file.attrs["tictac_key_hash_full"] = key_hash_full
    h5file.attrs["tictac_key_json"] = key_json
    h5file.attrs["tictac_created_utc"] = now_utc_iso8601()
    h5file.attrs["tictac_writer_version"] = WRITER_VERSION


# ---------------------------------------------------------------------------
# Manifest helpers
# ---------------------------------------------------------------------------
MANIFEST_VERSION = 1


def load_manifest(cache_root: str) -> dict:
    path = os.path.join(cache_root, "manifest.json")
    if not os.path.exists(path):
        return {"manifest_version": MANIFEST_VERSION, "entries": []}
    with open(path) as f:
        try:
            return json.load(f)
        except json.JSONDecodeError:
            return {"manifest_version": MANIFEST_VERSION, "entries": []}


def save_manifest(cache_root: str, manifest: dict) -> None:
    path = os.path.join(cache_root, "manifest.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(manifest, f, indent=2)
    os.replace(tmp, path)


def upsert_manifest_entry(manifest: dict, entry: dict) -> bool:
    """Insert or update by key_hash_full. Returns True if inserted."""
    for i, e in enumerate(manifest["entries"]):
        if e.get("key_hash_full") == entry["key_hash_full"]:
            manifest["entries"][i] = entry
            return False
    manifest["entries"].append(entry)
    return True


# ---------------------------------------------------------------------------
# Core migration logic for a single file
# ---------------------------------------------------------------------------
def migrate_file(
    src_path: str,
    target_root: str,
    apply: bool,
    copy_mode: bool,
    manifest: dict,
    verbose: bool,
) -> str:
    """
    Process one legacy h5 file.

    Returns a one-line summary string describing the action taken (or planned).
    """
    fname = os.path.basename(src_path)
    m = LEGACY_P123_RE.match(fname)
    if m is None:
        return f"  SKIP  {src_path}  (filename does not match legacy pattern)"

    two_J_3N = int(m.group("two_J_3N"))
    P_3N = int(m.group("P_3N"))
    Np = int(m.group("Np"))
    Nq = int(m.group("Nq"))
    J2max = int(m.group("J2max"))

    Nphi = LEGACY_DEFAULTS["Nphi"]
    Nx = LEGACY_DEFAULTS["Nx"]
    tensor_force = LEGACY_DEFAULTS["tensor_force"]
    isospin_breaking_1S0 = LEGACY_DEFAULTS["isospin_breaking_1S0"]
    chebyshev_s = LEGACY_DEFAULTS["chebyshev_s"]
    chebyshev_t = LEGACY_DEFAULTS["chebyshev_t"]

    key_json = canonical_p123_json(
        schema_version=P123_SCHEMA_VERSION,
        Np_WP=Np,
        Nq_WP=Nq,
        J_2N_max=J2max,
        two_J_3N=two_J_3N,
        P_3N=P_3N,
        Nphi=Nphi,
        Nx=Nx,
        tensor_force=tensor_force,
        isospin_breaking_1S0=isospin_breaking_1S0,
        chebyshev_s=chebyshev_s,
        chebyshev_t=chebyshev_t,
    )

    key_hash_full = sha256_hex(key_json)
    h_short = hash_short(key_hash_full)
    fp = filename_prefix(Np, Nq, J2max, two_J_3N, P_3N)
    dest_fname = f"{fp}__{h_short}.h5"
    dest_rel = os.path.join("p123", dest_fname)
    dest_path = os.path.join(target_root, dest_rel)

    action = "COPY" if copy_mode else "MOVE"
    summary = f"  {action}  {src_path}\n        -> {dest_path}\n        hash={key_hash_full[:12]}..."

    if verbose:
        print(f"  key_json: {key_json}")

    if apply:
        os.makedirs(os.path.join(target_root, "p123"), exist_ok=True)

        if os.path.exists(dest_path):
            try:
                with h5py.File(dest_path, "r") as hf:
                    raw = hf.attrs.get("tictac_key_hash_full", "")
                existing_hash = raw.decode() if isinstance(raw, bytes) else str(raw)
            except Exception:
                existing_hash = ""
            if existing_hash == key_hash_full:
                summary += "\n        SKIP (target already exists with matching hash)"
                return summary
            summary += (f"\n        ERROR collision: target exists with different hash"
                        f" ({existing_hash[:12]}... != {key_hash_full[:12]}...); skipped")
            return summary

        shutil.copy2(src_path, dest_path)

        # Write new root attributes
        with h5py.File(dest_path, "a") as hf:
            write_root_attrs(hf, key_hash_full, key_json)

        # Remove source if move mode
        if not copy_mode:
            os.remove(src_path)

        # Update manifest
        size_bytes = os.path.getsize(dest_path)
        utc_now = now_utc_iso8601()
        entry = {
            "kind": "p123",
            "filename": dest_rel,
            "key_hash_full": key_hash_full,
            "key": json.loads(key_json),
            "key_json": key_json,
            "size_bytes": size_bytes,
            "schema_version": P123_SCHEMA_VERSION,
            "created_utc": utc_now,
            "last_used_utc": utc_now,
            "legacy_imported": True,
            "writer_version": WRITER_VERSION,
        }
        inserted = upsert_manifest_entry(manifest, entry)
        status = "inserted" if inserted else "updated"
        summary += f"\n        manifest {status}"

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Migrate legacy P123_sparse_JP_*.h5 files into the new cache layout."
    )
    p.add_argument(
        "--legacy-roots",
        nargs="+",
        required=True,
        metavar="DIR",
        help="Root directories to scan recursively for legacy h5 files.",
    )
    p.add_argument(
        "--target-root",
        required=True,
        metavar="DIR",
        help="Destination cache root directory.",
    )
    p.add_argument(
        "--apply",
        action="store_true",
        default=False,
        help="Execute writes (default: dry-run only).",
    )
    p.add_argument(
        "--copy",
        action="store_true",
        default=False,
        help="Copy files instead of moving them.",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Print extra debug information.",
    )
    return p.parse_args()


def collect_legacy_files(roots: list) -> list:
    found = []
    for root in roots:
        if not os.path.isdir(root):
            print(f"WARNING: legacy root not found, skipping: {root}", file=sys.stderr)
            continue
        for dirpath, _dirnames, filenames in os.walk(root):
            for fname in filenames:
                if LEGACY_P123_RE.match(fname):
                    found.append(os.path.join(dirpath, fname))
    return sorted(found)


def main() -> None:
    args = parse_args()
    mode_label = "APPLY" if args.apply else "DRY-RUN"
    action_label = "COPY" if args.copy else "MOVE"
    print(f"migrate_legacy_cache  [{mode_label}]  action={action_label}")
    print(f"  legacy roots : {args.legacy_roots}")
    print(f"  target root  : {args.target_root}")
    print()

    legacy_files = collect_legacy_files(args.legacy_roots)
    if not legacy_files:
        print("No legacy P123_sparse_JP_*.h5 files found.")
        if not args.apply:
            print("(dry-run; pass --apply to execute)")
        return

    print(f"Found {len(legacy_files)} legacy file(s):")

    # Load or initialise manifest
    if args.apply:
        os.makedirs(args.target_root, exist_ok=True)
    manifest = load_manifest(args.target_root)

    n_ok = 0
    n_skip = 0
    for src in legacy_files:
        summary = migrate_file(
            src_path=src,
            target_root=args.target_root,
            apply=args.apply,
            copy_mode=args.copy,
            manifest=manifest,
            verbose=args.verbose,
        )
        print(summary)
        if "SKIP" in summary:
            n_skip += 1
        else:
            n_ok += 1

    if args.apply:
        save_manifest(args.target_root, manifest)
        print(f"\nDone: {n_ok} migrated, {n_skip} skipped.")
    else:
        print(f"\n{n_ok} candidate(s), {n_skip} skipped.")
        print("(dry-run; pass --apply to execute)")


if __name__ == "__main__":
    main()
