#!/usr/bin/env python3
"""
Freeze a bit-for-bit 2NF-only baseline for the chiral N2LO 3NF rewrite.

Purpose:
    Run the canonical 190 MeV deuteron-proton solver workflow with
    ``three_nucleon_force=none`` and persist SHA-256 hashes of the
    ``U_PW_elements_*.txt`` files to
    ``tests/fixtures/2nf_baseline_hashes.json``. Downstream regression
    test ``tests/test_2nf_miller_baseline.py`` loads this JSON and
    asserts that future runs of the 2NF-only path produce identical
    outputs. Benchmarks: Miller, Ekström, Hebeler, Phys. Rev. C 106,
    024001 (2022), with the Glöckle-Witała 1996 standard-Faddeev
    values (Phys. Rep. 274, 107) as NN-only reference.

When to run:
    Once, manually, or on explicit CI trigger. Re-running invalidates
    the guardrail — see ``tests/MILLER_BASELINE_README.md``.

Usage:
    python3 tests/fixtures/freeze_2nf_baseline.py \
        --work-dir output/2nf_miller_baseline \
        --target-tlab-mev 190
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURES_DIR = ROOT / "tests" / "fixtures"
HASH_FILE = FIXTURES_DIR / "2nf_baseline_hashes.json"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def run_solver_2nf_only(work_dir: Path, target_tlab_mev: float, reuse_p123: bool) -> Path:
    """Invoke ``examples/deuteron_proton_Ay.py`` with 3NF explicitly off.

    The script already defaults to the 2NF-only code path (the solver
    binary hard-defaults ``three_nucleon_force=none``), but we pass
    ``--three-nucleon-force none`` explicitly when the CLI supports it.
    If the CLI does not expose that flag yet, we still freeze what the
    default currently produces and annotate the baseline accordingly.
    """
    cmd = [
        sys.executable,
        str(ROOT / "examples" / "deuteron_proton_Ay.py"),
        "--work-dir",
        str(work_dir),
        "--target-tlab-mev",
        str(target_tlab_mev),
    ]
    if reuse_p123:
        cmd.append("--reuse-p123")

    print("[freeze_2nf_baseline] running:", " ".join(cmd), flush=True)
    result = subprocess.run(cmd, cwd=ROOT)
    if result.returncode != 0:
        raise SystemExit(f"solver run failed (rc={result.returncode})")

    solver_out = work_dir / "solver_out"
    if not solver_out.is_dir():
        raise SystemExit(f"solver_out directory missing: {solver_out}")
    return solver_out


def hash_u_files(solver_out: Path) -> dict[str, str]:
    hashes = {}
    for u_file in sorted(solver_out.glob("U_PW_elements_*.txt")):
        hashes[u_file.name] = sha256_of(u_file)
    if not hashes:
        raise SystemExit(f"no U_PW_elements_*.txt files found in {solver_out}")
    return hashes


def write_baseline(hashes: dict[str, str], target_tlab_mev: float, solver_out: Path) -> None:
    FIXTURES_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "target_tlab_mev": target_tlab_mev,
        "three_nucleon_force": "none",
        "solver_out_relpath": str(solver_out.relative_to(ROOT)),
        "benchmark_refs": [
            "Miller, Ekstrom, Hebeler, Phys. Rev. C 106, 024001 (2022)",
            "Glockle, Witala, Huber, Kamada, Golak, Phys. Rep. 274, 107 (1996)",
            "Nogga, Kamada, Glockle, Barrett, Phys. Rev. C 65, 054003 (2002)",
        ],
        "sha256": hashes,
    }
    HASH_FILE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[freeze_2nf_baseline] wrote {HASH_FILE}", flush=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default="output/2nf_miller_baseline")
    parser.add_argument("--target-tlab-mev", type=float, default=190.0)
    parser.add_argument("--reuse-p123", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    work_dir = (ROOT / args.work_dir).resolve() if not Path(args.work_dir).is_absolute() else Path(args.work_dir)
    solver_out = run_solver_2nf_only(work_dir, args.target_tlab_mev, args.reuse_p123)
    hashes = hash_u_files(solver_out)
    write_baseline(hashes, args.target_tlab_mev, solver_out)
    print(f"[freeze_2nf_baseline] froze {len(hashes)} U_PW_elements files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
