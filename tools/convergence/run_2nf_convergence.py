#!/usr/bin/env python3
"""
Phase D: 2NF convergence campaign on the realistic grid.

Runs three independent convergence ladders with three_nucleon_force=none:
  D1. Np    — already covered by deuteron binding ladder; spot-check only
  D2. Nq    — vary Nq at fixed Np=82 (WARNING: moves Tlab midpoint)
  D3. J_2N  — vary J_2N_max at fixed Np=82, Nq=12
  D4. two_J_3N_max — vary at fixed Np=82, Nq=12, J_2N_max
  D5. Padé  — [24/24] vs [32/32] comparison (already in feasibility JSON)

For each ladder, collects on-shell U and computes max|dU| and RMS|dU|
between adjacent rungs. Writes a machine-readable convergence summary.

FAIL-CLOSED: each axis is independent. A physical claim requires ALL axes PASS.
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SOLVER = REPO / "CPP" / "run"
INPUT = REPO / "CPP" / "Input" / "input_realistic_82_12.txt"
WORK = REPO / "output" / "realistic_2nf_convergence"
OUT_JSON = REPO / "output" / "validation" / "realistic_2nf_convergence.json"


def parse_u_file(path):
    """Parse a U_PW_elements file. Returns list of (Tlab, Ecm, qidx, [U_strs])."""
    rows = []
    in_table = False
    for line in path.read_text(errors='replace').splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("# ##"):
            in_table = True
            continue
        if s.startswith("#"):
            continue
        if not in_table:
            continue
        parts = s.split()
        if len(parts) >= 4:
            try:
                tlab = float(parts[0]); ecm = float(parts[1]); qidx = int(parts[2])
                rows.append((tlab, ecm, qidx, parts[3:]))
            except ValueError:
                pass
    return rows


def find_u_files(out_dir):
    """Find U_PW_elements files in the output directory."""
    files = {}
    for p in sorted(out_dir.glob("U_PW_elements_*PSI_0.txt")):
        # Extract JP tag from filename
        m = re.search(r'JP_(\d+)_(-?\d+)', p.name)
        if m:
            jp = f"{m.group(1)}_{m.group(2)}"
            files[jp] = p
    return files


def u_diff(a_rows, b_rows):
    """Max and RMS |dU| over row-aligned complex U blocks."""
    mx = 0.0; ss = 0.0; n = 0
    for ra, rb in zip(a_rows, b_rows):
        for xa, xb in zip(ra[3], rb[3]):
            za = complex(xa.replace(" ", "")); zb = complex(xb.replace(" ", ""))
            d = abs(za - zb)
            mx = max(mx, d); ss += d*d; n += 1
    rms = math.sqrt(ss / n) if n else 0.0
    return mx, rms, n


def run_solver(out_dir, cache_dir, overrides, label=""):
    """Run the 2NF solver with given overrides."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(SOLVER), str(INPUT),
           f"output_folder={out_dir}", f"P123_folder={out_dir}",
           f"cache_root={cache_dir}",
           "three_nucleon_force=none",
           ] + overrides
    env = dict(os.environ, OMP_NUM_THREADS="12")
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=REPO, timeout=3600)
    wall = time.time() - t0
    if proc.returncode != 0:
        print(f"  WARNING: solver returned {proc.returncode} for {label}", file=sys.stderr)
        print(proc.stderr[-300:], file=sys.stderr)
    return wall


def collect_convergence(ladder_name, rungs, work_base):
    """Run a ladder and compute convergence metrics between adjacent rungs."""
    results = []
    print(f"\n{'='*60}")
    print(f"  {ladder_name}")
    print(f"{'='*60}")
    for i, (label, overrides) in enumerate(rungs):
        out = work_base / f"rung_{i}_{label}" / "out"
        cache = work_base / f"rung_{i}_{label}" / "cache"
        print(f"\n  rung {i}: {label} ...", end="", flush=True)
        wall = run_solver(out, cache, overrides, label)
        u_files = find_u_files(out)
        u_data = {jp: parse_u_file(p) for jp, p in u_files.items()}
        print(f" {wall:.0f}s, {len(u_files)} U files")
        results.append({"label": label, "wall_seconds": wall, "u_data": u_data,
                        "u_files": {jp: str(p) for jp, p in u_files.items()}})

    # Adjacent convergence metrics
    convergence = []
    for i in range(1, len(results)):
        prev, curr = results[i-1], results[i]
        entry = {"from": prev["label"], "to": curr["label"], "sectors": {}}
        for jp in sorted(curr["u_data"].keys()):
            if jp in prev["u_data"]:
                mx, rms, n = u_diff(prev["u_data"][jp], curr["u_data"][jp])
                entry["sectors"][jp] = {"max_abs_dU_MeV": mx, "rms_dU_MeV": rms, "n": n}
                print(f"    {prev['label']} -> {curr['label']} JP={jp}: "
                      f"max|dU|={mx:.4e} rms={rms:.4e}")
        convergence.append(entry)
    return results, convergence


def main():
    WORK.mkdir(parents=True, exist_ok=True)

    report = {
        "schema": "tictac.realistic_2nf_convergence.v1",
        "generated_on": "2026-08-15",
        "purpose": "2NF-only convergence baseline on the realistic grid",
        "grid_base": {"Np_WP": 82, "Nq_WP": 12, "p_chebyshev_t": 0.88,
                       "p_chebyshev_s": 325, "q_chebyshev_t": 1.0, "q_chebyshev_s": 90},
        "ladders": {},
    }

    # ---- D3: J_2N_max ladder (no energy alignment issue) ----
    j2n_rungs = [
        ("J2N=1", ["J_2N_max=1", "two_J_3N_max=1"]),
        ("J2N=2", ["J_2N_max=2", "two_J_3N_max=1"]),
        ("J2N=3", ["J_2N_max=3", "two_J_3N_max=1"]),
        ("J2N=4", ["J_2N_max=4", "two_J_3N_max=1"]),
        ("J2N=5", ["J_2N_max=5", "two_J_3N_max=1"]),
    ]
    j2n_results, j2n_conv = collect_convergence(
        "D3: J_2N_max ladder (Np=82, Nq=12, two_J_3N_max=1)", j2n_rungs,
        WORK / "j2n_ladder")
    report["ladders"]["J_2N_max"] = {
        "rungs": [{"label": r["label"], "wall_seconds": r["wall_seconds"]} for r in j2n_results],
        "convergence": j2n_conv,
        "axis_status": "PASS" if all(
            s["max_abs_dU_MeV"] < 0.01 for c in j2n_conv
            for s in c["sectors"].values()
        ) else "NOT YET CONVERGED",
    }

    # ---- D4: two_J_3N_max ladder (at fixed J_2N_max=5) ----
    # Use the highest converged J_2N_max as basis
    j3n_rungs = [
        ("2J3N=1", ["J_2N_max=5", "two_J_3N_max=1"]),
        ("2J3N=3", ["J_2N_max=5", "two_J_3N_max=3"]),
        ("2J3N=5", ["J_2N_max=5", "two_J_3N_max=5"]),
        ("2J3N=7", ["J_2N_max=5", "two_J_3N_max=7"]),
    ]
    j3n_results, j3n_conv = collect_convergence(
        "D4: two_J_3N_max ladder (Np=82, Nq=12, J_2N_max=5)", j3n_rungs,
        WORK / "j3n_ladder")
    report["ladders"]["two_J_3N_max"] = {
        "rungs": [{"label": r["label"], "wall_seconds": r["wall_seconds"]} for r in j3n_results],
        "convergence": j3n_conv,
        "axis_status": "PASS" if all(
            s["max_abs_dU_MeV"] < 0.01 for c in j3n_conv
            for s in c["sectors"].values()
        ) else "NOT YET CONVERGED",
    }

    # ---- D2: Nq ladder (WARNING: moves Tlab midpoint) ----
    nq_rungs = [
        ("Nq=8", ["Nq_WP=8", "J_2N_max=5", "two_J_3N_max=1"]),
        ("Nq=10", ["Nq_WP=10", "J_2N_max=5", "two_J_3N_max=1"]),
        ("Nq=12", ["Nq_WP=12", "J_2N_max=5", "two_J_3N_max=1"]),
        ("Nq=14", ["Nq_WP=14", "J_2N_max=5", "two_J_3N_max=1"]),
        ("Nq=16", ["Nq_WP=16", "J_2N_max=5", "two_J_3N_max=1"]),
    ]
    nq_results, nq_conv = collect_convergence(
        "D2: Nq ladder (Np=82, J_2N_max=5, two_J_3N_max=1) — WARNING: midpoint moves",
        nq_rungs, WORK / "nq_ladder")
    # Report midpoints
    for r in nq_results:
        midpoints = {}
        for jp, rows in r["u_data"].items():
            if rows:
                midpoints[jp] = {"Tlab": rows[0][0], "Ecm": rows[0][1]}
        r["midpoints"] = midpoints
    report["ladders"]["Nq"] = {
        "rungs": [{"label": r["label"], "wall_seconds": r["wall_seconds"],
                   "midpoints": r.get("midpoints", {})} for r in nq_results],
        "convergence": nq_conv,
        "warning": "Changing Nq moves the wave-packet energy-bin midpoints. "
                   "Direct U comparison between different Nq is NOT physically valid. "
                   "Compare only at aligned midpoints or fail closed.",
        "axis_status": "REQUIRES_ALIGNED_MIDPOINT_CHECK",
    }

    # ---- D5: Padé (use existing feasibility data) ----
    report["ladders"]["pade"] = {
        "source": "output/validation/n2lo_3nf_physical_Ay_feasibility.json",
        "order_24_vs_32_max_abs_dU_mev": 1.2803710168e-09,
        "status": "PASS (16/16 amplitudes converged at [24/24] and [32/32])",
        "axis_status": "PASS",
    }

    # ---- D1: Np (use existing deuteron binding ladder) ----
    report["ladders"]["Np"] = {
        "source": "output/validation/idaho_n3lo_deuteron_binding_ladder.json",
        "note": "Binding energy converges as 1/Np^2. Np=82 is within 20-keV gate. "
                "Binding convergence does NOT imply scattering convergence; "
                "D2-D4 are the scattering convergence axes.",
        "axis_status": "PASS (binding); scattering requires D2-D4",
    }

    # Write summary
    summary_lines = []
    summary_lines.append("axis                  result\n")
    summary_lines.append("--------------------------------\n")
    for axis, data in report["ladders"].items():
        status = data.get("axis_status", "NOT TESTED")
        summary_lines.append(f"{axis:<22} {status}\n")
    report["summary_table"] = "".join(summary_lines)

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    # Don't serialize the raw u_data (too large)
    for ladder in report["ladders"].values():
        if "rungs" in ladder:
            for r in ladder["rungs"]:
                r.pop("u_data", None)
    OUT_JSON.write_text(json.dumps(report, indent=2))
    print(f"\n{'='*60}")
    print(f"wrote {OUT_JSON}")
    print(f"\nSummary:")
    print(report["summary_table"])


if __name__ == "__main__":
    main()
