#!/usr/bin/env python3
"""
Analyze J_2N_max convergence from the Sean-compatible campaign.
Compares on-shell U elements between adjacent rungs at Tlab≈10 MeV.
"""
import json, math, re, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO / "output" / "sean_2nf_convergence" / "j2n_ladder"
OUT_JSON = REPO / "output" / "validation" / "sean_j2n_convergence.json"

def parse_u_file(path):
    """Parse U file. Returns list of (Tlab, Ecm, q_idx, [U_strs])."""
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
            try:
                rows.append((float(parts[0]), float(parts[1]), int(parts[2]), parts[3:]))
            except ValueError: pass
    return rows

def u_diff(a, b):
    mx = 0.0; ss = 0.0; n = 0; max_rel = 0.0
    for ra, rb in zip(a, b):
        for xa, xb in zip(ra[3], rb[3]):
            za = complex(xa.replace(" ", "")); zb = complex(xb.replace(" ", ""))
            d = abs(za - zb)
            mx = max(mx, d); ss += d*d; n += 1
            if abs(za) > 1e-10:
                max_rel = max(max_rel, d / abs(za))
    rms = math.sqrt(ss / n) if n else 0.0
    return {"max_abs_dU_MeV": mx, "rms_dU_MeV": rms, "n": n,
            "max_relative": max_rel}

def main():
    # Collect data for each J2N rung
    rungs = {}
    for j2n in [1, 2, 3, 4]:
        d = CAMPAIGN / f"J2N_{j2n}" / "out"
        if not d.exists(): continue
        u_files = {}
        for p in sorted(d.glob("U_PW_elements_*PSI_0.txt")):
            m = re.search(r'JP_(\d+)_(-?\d+)', p.name)
            if m:
                jp = f"{m.group(1)}_{m.group(2)}"
                u_files[jp] = parse_u_file(p)
        if u_files:
            rungs[j2n] = u_files
            # Report energies
            for jp, rows in u_files.items():
                if rows:
                    energies = [r[0] for r in rows]
                    print(f"J2N={j2n} JP={jp}: {len(rows)} energies, "
                          f"Tlab range [{min(energies):.2f}, {max(energies):.2f}] MeV")

    # Compare adjacent rungs at each energy
    report = {"schema": "tictac.sean_j2n_convergence.v1", "generated_on": "2026-08-15",
              "grid": {"Np_WP": 82, "Nq_WP": 12, "two_J_3N_max": 1, "force": "none"},
              "comparisons": []}

    sorted_j2n = sorted(rungs.keys())
    for i in range(1, len(sorted_j2n)):
        prev_j2n = sorted_j2n[i-1]
        curr_j2n = sorted_j2n[i]
        print(f"\n=== J2N={prev_j2n} -> J2N={curr_j2n} ===")
        comp = {"from": f"J2N={prev_j2n}", "to": f"J2N={curr_j2n}", "by_energy": []}
        for jp in sorted(rungs[curr_j2n].keys()):
            if jp not in rungs[prev_j2n]: continue
            prev_rows = rungs[prev_j2n][jp]
            curr_rows = rungs[curr_j2n][jp]
            # Match by Tlab (closest)
            for prow in prev_rows:
                tlab_p = prow[0]
                # Find closest in curr
                best = min(curr_rows, key=lambda r: abs(r[0] - tlab_p))
                if abs(best[0] - tlab_p) / tlab_p > 0.01:  # >1% mismatch
                    continue
                d = u_diff([prow], [best])
                d["Tlab_MeV"] = tlab_p
                d["JP"] = jp
                comp["by_energy"].append(d)
                print(f"  JP={jp} Tlab={tlab_p:.4f}: max|dU|={d['max_abs_dU_MeV']:.4e} "
                      f"rms={d['rms_dU_MeV']:.4e} max_rel={d['max_relative']:.4e}")
        report["comparisons"].append(comp)

    # Summary at Tlab≈10 MeV
    print("\n=== Summary at Tlab≈10 MeV ===")
    for comp in report["comparisons"]:
        for d in comp["by_energy"]:
            if abs(d["Tlab_MeV"] - 10.0) < 2.0:
                print(f"  {comp['from']} -> {comp['to']} JP={d['JP']} Tlab={d['Tlab_MeV']:.2f}: "
                      f"max|dU|={d['max_abs_dU_MeV']:.4e} MeV")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {OUT_JSON}")

if __name__ == "__main__":
    main()
