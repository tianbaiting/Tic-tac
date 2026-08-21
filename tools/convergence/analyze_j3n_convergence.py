#!/usr/bin/env python3
"""Analyze two_J_3N_max convergence from the Sean-compatible campaign."""
import json, math, re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BASE = REPO / "output" / "sean_2nf_convergence" / "j3n_ladder"
OUT_JSON = REPO / "output" / "validation" / "sean_j3n_convergence.json"

def parse_u(path):
    rows = []
    in_table = False
    for line in Path(path).read_text(errors='replace').splitlines():
        s = line.strip()
        if not s: continue
        if s.startswith("# ##"): in_table = True; continue
        if s.startswith("#"): continue
        if not in_table: continue
        parts = s.split()
        if len(parts) >= 4:
            try: rows.append((float(parts[0]), float(parts[1]), int(parts[2]), parts[3:]))
            except ValueError: pass
    return rows

def u_diff(a, b):
    mx = 0.0; ss = 0.0; n = 0; max_rel = 0.0
    for ra, rb in zip(a, b):
        for xa, xb in zip(ra[3], rb[3]):
            za = complex(xa.replace(" ", "")); zb = complex(xb.replace(" ", ""))
            d = abs(za - zb)
            mx = max(mx, d); ss += d*d; n += 1
            if abs(za) > 1e-10: max_rel = max(max_rel, d / abs(za))
    return {"max_abs_dU_MeV": mx, "rms_dU_MeV": math.sqrt(ss/n) if n else 0, "n": n, "max_relative": max_rel}

rungs = {}
for tj3 in [1, 3, 5, 7, 9, 11, 13, 15, 17]:
    d = BASE / f"J3N_{tj3}" / "out"
    if not d.exists(): continue
    u_files = {}
    for p in sorted(d.glob("U_PW_elements_*PSI_0.txt")):
        m = re.search(r'JP_(\d+)_(-?\d+)', p.name)
        if m:
            jp = f"{m.group(1)}_{m.group(2)}"
            u_files[jp] = parse_u(p)
    if u_files:
        rungs[tj3] = u_files

# Compare adjacent rungs
report = {"schema": "tictac.sean_j3n_convergence.v1", "generated_on": "2026-08-16",
          "grid": {"Np_WP": 82, "Nq_WP": 12, "J_2N_max": 3, "force": "none"},
          "comparisons": [], "summary_at_10MeV": []}

sorted_tj3 = sorted(rungs.keys())
for i in range(1, len(sorted_tj3)):
    prev_tj3 = sorted_tj3[i-1]
    curr_tj3 = sorted_tj3[i]
    print(f"\n=== 2J3N={prev_tj3} -> 2J3N={curr_tj3} ===")
    comp = {"from": f"2J3N={prev_tj3}", "to": f"2J3N={curr_tj3}", "by_energy": []}
    for jp in sorted(rungs[curr_tj3].keys()):
        if jp not in rungs[prev_tj3]: continue
        prev_rows = rungs[prev_tj3][jp]
        curr_rows = rungs[curr_tj3][jp]
        for prow in prev_rows:
            tlab = prow[0]
            best = min(curr_rows, key=lambda r: abs(r[0] - tlab))
            if abs(best[0] - tlab) / tlab > 0.01: continue
            d = u_diff([prow], [best])
            d["Tlab_MeV"] = tlab; d["JP"] = jp
            comp["by_energy"].append(d)
            if abs(tlab - 10.0) < 2.0:
                print(f"  JP={jp} Tlab={tlab:.4f}: max|dU|={d['max_abs_dU_MeV']:.4e} rms={d['rms_dU_MeV']:.4e}")
                report["summary_at_10MeV"].append({
                    "from": f"2J3N={prev_tj3}", "to": f"2J3N={curr_tj3}",
                    "JP": jp, "Tlab_MeV": tlab,
                    "max_abs_dU_MeV": d["max_abs_dU_MeV"], "rms_dU_MeV": d["rms_dU_MeV"]
                })
    report["comparisons"].append(comp)

print("\n=== Summary at Tlab~10 MeV ===")
print(f"{'transition':<25} {'JP':>8} {'max|dU| (MeV)':>15} {'rms (MeV)':>12}")
print("-" * 65)
for s in report["summary_at_10MeV"]:
    print(f"{s['from']} -> {s['to']:<10} {s['JP']:>8} {s['max_abs_dU_MeV']:>15.4e} {s['rms_dU_MeV']:>12.4e}")

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
OUT_JSON.write_text(json.dumps(report, indent=2))
print(f"\nwrote {OUT_JSON}")
