#!/usr/bin/env python3
"""
Proper J3N convergence analysis: sum U across all J^pi sectors at each rung,
then compare the cumulative sum between adjacent rungs.

The elastic nd scattering amplitude is the sum over all J^pi sectors.
Adding a new J sector changes the total by exactly that sector's contribution.
"""
import json, math, re
from pathlib import Path
from collections import defaultdict

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

def sum_u_across_sectors(rung_dir):
    """Sum U elements across all J^pi sectors, matching by (Tlab, Ecm, q_idx, element_idx)."""
    u_files = sorted(rung_dir.glob("U_PW_elements_*PSI_0.txt"))
    if not u_files:
        return None

    # Collect all U values keyed by (Tlab, element_index)
    # Each file has the same structure: rows of (Tlab, Ecm, q_idx, [U00, U01, ...])
    # The total elastic U = sum of U across all J^pi sectors
    cumulative = defaultdict(lambda: defaultdict(complex))  # Tlab -> elem_idx -> sum
    meta = {}  # Tlab -> (Ecm, q_idx)
    sectors = []

    for p in u_files:
        m = re.search(r'JP_(\d+)_(-?\d+)', p.name)
        if not m: continue
        jp = f"{m.group(1)}_{m.group(2)}"
        rows = parse_u(p)
        sectors.append(jp)
        for tlab, ecm, qidx, u_strs in rows:
            meta[tlab] = (ecm, qidx)
            for i, us in enumerate(u_strs):
                cumulative[tlab][i] += complex(us.replace(" ", ""))

    # Convert to sorted list
    result = []
    for tlab in sorted(meta.keys()):
        ecm, qidx = meta[tlab]
        u_vals = [cumulative[tlab][i] for i in range(len(cumulative[tlab]))]
        result.append((tlab, ecm, qidx, u_vals))
    return result, sectors

def cumulative_diff(a, b):
    mx = 0.0; ss = 0.0; n = 0; max_rel = 0.0
    for (ta, ea, qa, ua), (tb, eb, qb, ub) in zip(a, b):
        assert abs(ta - tb) / ta < 0.01, f"Tlab mismatch: {ta} vs {tb}"
        for za, zb in zip(ua, ub):
            d = abs(za - zb)
            mx = max(mx, d); ss += d*d; n += 1
            denom = max(abs(za), 1e-10)
            max_rel = max(max_rel, d / denom)
    return {"max_abs_dU_MeV": mx, "rms_dU_MeV": math.sqrt(ss/n) if n else 0,
            "n": n, "max_relative": max_rel}

def main():
    rungs = {}
    for tj3 in [1, 3, 5, 7, 9, 11, 13, 15, 17]:
        d = BASE / f"J3N_{tj3}" / "out"
        if not d.exists(): continue
        result = sum_u_across_sectors(d)
        if result:
            cum_u, sectors = result
            rungs[tj3] = (cum_u, sectors)
            print(f"2J3N={tj3:2d}: {len(sectors)} sectors ({', '.join(sectors)}), {len(cum_u)} energies")

    # Compare cumulative U between adjacent rungs
    report = {"schema": "tictac.sean_j3n_convergence.v1", "generated_on": "2026-08-26",
              "grid": {"Np_WP": 82, "Nq_WP": 12, "J_2N_max": 3, "force": "none"},
              "method": "Sum U across all J^pi sectors; compare cumulative total between rungs",
              "comparisons": [], "summary_at_10MeV": []}

    sorted_tj3 = sorted(rungs.keys())
    for i in range(1, len(sorted_tj3)):
        prev_tj3 = sorted_tj3[i-1]
        curr_tj3 = sorted_tj3[i]
        prev_u, prev_sectors = rungs[prev_tj3]
        curr_u, curr_sectors = rungs[curr_tj3]
        new_sectors = [s for s in curr_sectors if s not in prev_sectors]

        print(f"\n=== 2J3N={prev_tj3} -> 2J3N={curr_tj3} (new sectors: {new_sectors}) ===")
        comp = {"from": f"2J3N={prev_tj3}", "to": f"2J3N={curr_tj3}",
                "new_sectors": new_sectors, "by_energy": []}

        d = cumulative_diff(prev_u, curr_u)
        print(f"  Cumulative change across all energies:")
        print(f"    max|dU|={d['max_abs_dU_MeV']:.4e} rms={d['rms_dU_MeV']:.4e} max_rel={d['max_relative']:.4e}")

        # Per-energy breakdown
        for (tlab, ecm, qidx, ua), (_, _, _, ub) in zip(prev_u, curr_u):
            mx = max(abs(a-b) for a, b in zip(ua, ub))
            rms = math.sqrt(sum(abs(a-b)**2 for a, b in zip(ua, ub)) / len(ua))
            entry = {"Tlab_MeV": tlab, "Ecm_MeV": ecm,
                     "max_abs_dU_MeV": mx, "rms_dU_MeV": rms}
            comp["by_energy"].append(entry)
            if abs(tlab - 10.0) < 2.0:
                report["summary_at_10MeV"].append({
                    "from": f"2J3N={prev_tj3}", "to": f"2J3N={curr_tj3}",
                    "new_sectors": new_sectors,
                    "Tlab_MeV": tlab,
                    "max_abs_dU_MeV": mx, "rms_dU_MeV": rms
                })

        report["comparisons"].append(comp)

    # Summary table
    print(f"\n{'='*80}")
    print(f"Cumulative J3N convergence at Tlab~10 MeV")
    print(f"{'='*80}")
    print(f"{'transition':<25} {'new sectors':<20} {'max|dU| (MeV)':>15} {'rms (MeV)':>12}")
    print(f"{'-'*75}")
    for s in report["summary_at_10MeV"]:
        ns = ",".join(s["new_sectors"])
        print(f"{s['from']} -> {s['to']:<10} {ns:<20} {s['max_abs_dU_MeV']:>15.4e} {s['rms_dU_MeV']:>12.4e}")

    # Also print the full energy table for the last few transitions
    print(f"\n{'='*80}")
    print(f"Full energy dependence: 2J3N=15 -> 2J3N=17")
    print(f"{'='*80}")
    print(f"{'Tlab (MeV)':>12} {'max|dU| (MeV)':>15} {'rms (MeV)':>12}")
    print(f"{'-'*40}")
    last_comp = report["comparisons"][-1]
    for e in last_comp["by_energy"]:
        print(f"{e['Tlab_MeV']:>12.4f} {e['max_abs_dU_MeV']:>15.4e} {e['rms_dU_MeV']:>12.4e}")

    # Timing summary
    print(f"\n{'='*80}")
    print(f"Resource summary")
    print(f"{'='*80}")
    timings = {1:3296, 3:18771, 5:39727, 7:65674, 9:89549, 11:113663, 13:140951, 15:164852, 17:187685}
    total = sum(timings.values())
    for tj3 in sorted(timings.keys()):
        n_sectors = rungs.get(tj3, ([], []))[1]
        print(f"  2J3N={tj3:2d}: {timings[tj3]/3600:5.1f}h  {len(n_sectors):2d} sectors  cum={sum(timings[k] for k in timings if k<=tj3)/3600:6.1f}h")
    print(f"  total: {total/3600:.1f}h ({total/3600/24:.1f} days)")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {OUT_JSON}")

if __name__ == "__main__":
    main()
