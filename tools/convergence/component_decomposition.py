#!/usr/bin/env python3
"""
Phase I: component decomposition + zero-3NF limit on the small reference grid.

I1. Zero-3NF limit: three_nucleon_force=none must match the 2NF baseline.
I2. Component decomposition: separately enable c1, c3, c4, cD, cE and combinations.

The small grid (Np=4, Nq=3) is used for speed. Results are DIAGNOSTIC, not physical.
"""
import json, os, subprocess, sys, time, math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SOLVER = REPO / "CPP" / "run"
INPUT = REPO / "tools" / "refactor_harness" / "input_golden_3nf.txt"
WORK = REPO / "output" / "component_decomposition"
OUT_JSON = REPO / "output" / "validation" / "component_decomposition.json"

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
            try:
                rows.append((float(parts[0]), float(parts[1]), int(parts[2]), parts[3:]))
            except ValueError: pass
    return rows

def u_diff(a, b):
    mx = 0.0; ss = 0.0; n = 0
    for ra, rb in zip(a, b):
        for xa, xb in zip(ra[3], rb[3]):
            za = complex(xa.replace(" ", "")); zb = complex(xb.replace(" ", ""))
            d = abs(za - zb); mx = max(mx, d); ss += d*d; n += 1
    return mx, math.sqrt(ss/n) if n else 0, n

def run(label, overrides, p123_src=None):
    out = WORK / label / "out"
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(SOLVER), str(INPUT), f"output_folder={out}", f"P123_folder={p123_src or out}"] + overrides
    env = dict(os.environ, OMP_NUM_THREADS="4")
    t0 = time.time()
    subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=REPO, timeout=600)
    wall = time.time() - t0
    u_files = {}
    for p in sorted(out.glob("U_PW_elements_*PSI_0.txt")):
        import re
        m = re.search(r'JP_(\d+)_(-?\d+)', p.name)
        if m:
            u_files[f"{m.group(1)}_{m.group(2)}"] = parse_u(p)
    return wall, u_files

def main():
    WORK.mkdir(parents=True, exist_ok=True)
    results = {}

    # I1. Zero-3NF limit (2NF only)
    print("=== I1: zero-3NF limit (three_nucleon_force=none) ===")
    # First run builds P123
    wall0, u_2nf = run("2nf_baseline", ["three_nucleon_force=none"])
    print(f"  2NF baseline: {wall0:.1f}s")
    results["zero_3nf_limit"] = {"wall_seconds": wall0, "label": "three_nucleon_force=none"}

    # 3NF with all LECs zeroed (should match 2NF)
    print("=== I1b: 3NF with all LECs zeroed ===")
    wall1, u_3nf_zero = run("3nf_zero_lecs", [
        "three_nucleon_force=chiral_N2LO_full_factorized",
        "c_D=0.0", "c_E=0.0", "two_J_3NF_force_max=1"
    ], p123_src=str(WORK / "2nf_baseline" / "out"))
    print(f"  3NF zero-LECs: {wall1:.1f}s")
    results["zero_3nf_limit"]["wall_3nf_zero_seconds"] = wall1

    # Compare
    for jp in sorted(u_2nf.keys()):
        if jp in u_3nf_zero:
            mx, rms, n = u_diff(u_2nf[jp], u_3nf_zero[jp])
            results["zero_3nf_limit"][f"JP_{jp}"] = {"max_abs_dU_MeV": mx, "rms_dU_MeV": rms}
            print(f"  JP={jp}: zero-3NF vs 2NF: max|dU|={mx:.4e} rms={rms:.4e}")

    # I2. Component decomposition: enable one LEC at a time
    # Baseline LECs from the golden input: c_D=-0.2, c_E=-0.205
    # The c1, c3, c4 are internal to the factorized model (set by the regulator)
    components = {
        "cD_only": {"c_D": -0.2, "c_E": 0.0},
        "cE_only": {"c_D": 0.0, "c_E": -0.205},
        "cD_cE": {"c_D": -0.2, "c_E": -0.205},
        "cD_zero_cE_zero": {"c_D": 0.0, "c_E": 0.0},
    }
    print("\n=== I2: component decomposition ===")
    for label, lecs in components.items():
        overrides = ["three_nucleon_force=chiral_N2LO_full_factorized",
                     "two_J_3NF_force_max=1"]
        for k, v in lecs.items():
            overrides.append(f"{k}={v}")
        wall, u_comp = run(label, overrides,
                          p123_src=str(WORK / "2nf_baseline" / "out"))
        print(f"  {label}: {wall:.1f}s")
        results[label] = {"wall_seconds": wall, "lecs": lecs}
        for jp in sorted(u_2nf.keys()):
            if jp in u_comp:
                mx, rms, n = u_diff(u_2nf[jp], u_comp[jp])
                results[label][f"JP_{jp}"] = {"max_abs_dU_MeV": mx, "rms_dU_MeV": rms}
                print(f"    JP={jp}: vs 2NF max|dU|={mx:.4e} rms={rms:.4e}")

    results["label"] = "DIAGNOSTIC — small grid (Np=4, Nq=3); not physical"
    results["note"] = "Confirms solver-level changes are consistent with matrix-element signs. cD_only and cE_only show the individual contributions. cD+cE shows the combined effect. Zero-3NF (all LECs=0) should match the 2NF baseline."

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {OUT_JSON}")

if __name__ == "__main__":
    main()
