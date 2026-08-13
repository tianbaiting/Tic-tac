#!/usr/bin/env python3
"""
J_3NF convergence ladder driver (Phase C).

Runs the 3NF-active-J ladder

    Jc = (0), 1/2, 3/2, ...   ->   two_J_3NF_force_max = 2*Jc

on a FIXED grid/basis/W1-quadrature, reusing already-built W1 sectors, and
reports the incremental observable change

    delta_J U  =  U[2NF + 3NF active for J <= Jc] - U[2NF + 3NF active for J <= Jc_prev]

per on-shell elastic U element and per J^pi sector. Lower-J W1 sectors are cache
hits on every rung above the first, so the marginal cost is one sector's W1
build + a cheap solve.

FAIL-CLOSED: the J_3NF axis is independent of the Np/Nq / W1-quadrature / Pade /
J_2N axes. On toy grids every result is labelled DIAGNOSTIC, not physical; a
physical claim requires those other gates to be met first.

This is a thin orchestrator over the solver binary + the w1_worker; the physics
lives in the C++ core and the tictac observable package.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "python"))
from tictac.io import parse_u_file_meta  # noqa: E402

# On-shell elastic U values live in the data block of U_PW_elements_*.txt after
# the header lines. parse_u_file_meta gives metadata; the numeric rows are parsed
# here (same format examples/pw_amplitudes.py consumes).
def parse_u_values(path: Path):
    """Return list of (Tlab, Ecm, q_idx, [U00..U22 complex strings]) rows."""
    rows = []
    in_table = False
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            in_table = s.startswith("# ##") or in_table
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

def u_abs_diff_block(a_rows, b_rows):
    """Max |Delta U| and RMS over the complex U block, row-aligned."""
    mx = 0.0; ss = 0.0; n = 0
    for ra, rb in zip(a_rows, b_rows):
        for xa, xb in zip(ra[3], rb[3]):
            za = complex(xa.replace(" ", "")); zb = complex(xb.replace(" ", ""))
            d = abs(za - zb)
            mx = max(mx, d); ss += d*d; n += 1
    rms = math.sqrt(ss / n) if n else 0.0
    return mx, rms, n

def run_solver(input_txt: Path, out: Path, cache: Path, two_j_force_max: int, exe: str):
    out.mkdir(parents=True, exist_ok=True)
    subprocess.run([exe, str(input_txt),
                    f"output_folder={out}", f"P123_folder={out}", f"cache_root={cache}",
                    f"two_J_3NF_force_max={two_j_force_max}"],
                   cwd=REPO, check=True, capture_output=True)

def collect_u(out: Path):
    res = {}
    for p in sorted(out.glob("U_PW_elements_*_PSI_0.txt")):
        meta = parse_u_file_meta(p)
        if meta:
            # JP tag from filename: ..._JP_<twoJ>_<P>_Jmax_...
            tag = p.name
            jp = tag.split("JP_")[1].split("_Jmax")[0]  # "1_1" or "1_-1"
            res[jp] = (meta, parse_u_values(p))
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="solver input file (multi-sector grid)")
    ap.add_argument("--rungs", default="0,1,3",
                    help="comma list of two_J_3NF_force_max values (0=2NF baseline)")
    ap.add_argument("--exe", default=str(REPO / "CPP" / "run"))
    ap.add_argument("--work", default=str(REPO / "output" / "j3nf_ladder"))
    args = ap.parse_args()

    inp = Path(args.input)
    rungs = [int(x) for x in args.rungs.split(",")]
    work = Path(args.work)
    cache = work / "w1_db"
    per_rung = {}
    print(f"# J_3NF convergence ladder; grid fixed; reusing W1 db at {cache}")
    for r in rungs:
        out = work / f"cutoff_{r}"
        run_solver(inp, out, cache, r, args.exe)
        per_rung[r] = collect_u(out)
        sectors = sorted(per_rung[r].keys())
        print(f"  rung two_J_3NF_force_max={r:2d} (J<= {r/2:.1f}): sectors {sectors}")

    # Baseline = smallest rung (typically 0 = all 2NF).
    base = min(rungs)
    print("\n# Delta vs 2NF baseline (two_J_3NF_force_max=%d):" % base)
    report = {"baseline_cutoff": base, "rungs": [], "diagnostic": True,
              "note": "toy/diagnostic grid; physical claim needs Np/Nq/W1/Pade/J2N gates met"}
    for r in rungs:
        if r == base:
            continue
        entry = {"cutoff": r, "J_max_active": r / 2.0, "sectors": {}}
        for jp, (meta, rows) in sorted(per_rung[r].items()):
            if jp in per_rung[base]:
                mx, rms, n = u_abs_diff_block(rows, per_rung[base][jp][1])
                entry["sectors"][jp] = {"max_abs_dU_MeV": mx, "rms_dU_MeV": rms, "n": n}
                print(f"  cutoff={r} JP={jp:>5}: max|dU|={mx:.4e} MeV  rms={rms:.4e} (n={n})")
        # incremental vs previous rung
        prev = max(x for x in rungs if x < r) if any(x < r for x in rungs) else None
        if prev is not None:
            entry["delta_vs_prev"] = {}
            print(f"  # incremental vs previous rung (cutoff {prev}):")
            for jp in sorted(per_rung[r].keys()):
                if jp in per_rung[prev]:
                    mx, rms, n = u_abs_diff_block(per_rung[r][jp][1], per_rung[prev][jp][1])
                    entry["delta_vs_prev"][jp] = {"max_abs_dU_MeV": mx, "rms_dU_MeV": rms}
                    print(f"    cutoff {prev}->{r} JP={jp:>5}: max|dU|={mx:.4e} MeV  rms={rms:.4e}")
        report["rungs"].append(entry)

    out_json = work / "j3nf_convergence.json"
    out_json.write_text(json.dumps(report, indent=2))
    print(f"\n# wrote {out_json}")
    print("# DIAGNOSTIC ONLY: J_3NF axis isolated; other convergence gates not asserted here.")

if __name__ == "__main__":
    main()
