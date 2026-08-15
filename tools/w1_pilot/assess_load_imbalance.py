#!/usr/bin/env python3
"""
Phase C: assess static partition load imbalance.

Uses the measured block costs from the pilot (and, when available, the full
56-block build) to simulate the static modulo partition
(global_eval_index % N == worker_index) for N = 4, 8, 16, 32 workers.

Computes eta_load = mean_worker_load / max_worker_load.

If the full 56-block build is complete (output/realistic_w1_full/worker.log),
uses the MEASURED per-block costs. Otherwise, estimates from the 9 measured
pilot blocks using channel-property classification.
"""
import json
import math
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PILOT_JSON = REPO / "output" / "validation" / "realistic_w1_block_pilot.json"
FULL_BUILD_LOG = REPO / "output" / "realistic_w1_full" / "worker.log"
OUT_JSON = REPO / "output" / "validation" / "realistic_w1_load_imbalance.json"

# Channel properties for the 11 alphas in the J=1/2+ sector (chn=0, J_2N_max=1).
# From the w1_worker plan output.
ALPHA_PROPS = {
    0:  {"L2": 0, "S2": 0, "J2": 0, "L1": 0, "J1": 1, "T2": 1, "desc": "1S0 spect-S T=1"},
    1:  {"L2": 1, "S2": 1, "J2": 0, "L1": 0, "J1": 1, "T2": 1, "desc": "3P0 spect-S T=1"},
    2:  {"L2": 1, "S2": 0, "J2": 1, "L1": 1, "J1": 1, "T2": 0, "desc": "1P1 spect-P T=0"},
    3:  {"L2": 1, "S2": 0, "J2": 1, "L1": 1, "J1": 3, "T2": 0, "desc": "1P1 spect-P(J=3/2) T=0"},
    4:  {"L2": 0, "S2": 1, "J2": 1, "L1": 1, "J1": 1, "T2": 0, "desc": "3P1 spect-P T=0"},
    5:  {"L2": 0, "S2": 1, "J2": 1, "L1": 1, "J1": 3, "T2": 0, "desc": "3P1 spect-P(J=3/2) T=0"},
    6:  {"L2": 1, "S2": 1, "J2": 1, "L1": 1, "J1": 1, "T2": 1, "desc": "deuteron 3S1-3D1 T=1"},
    7:  {"L2": 1, "S2": 1, "J2": 1, "L1": 1, "J1": 3, "T2": 1, "desc": "deuteron 3S1-3D1 spect-P(J=3/2) T=1"},
    8:  {"L2": 2, "S2": 1, "J2": 1, "L1": 0, "J1": 1, "T2": 0, "desc": "3D1 spect-S T=0"},
    9:  {"L2": 2, "S2": 1, "J2": 1, "L1": 1, "J1": 3, "T2": 0, "desc": "3D1 spect-P(J=3/2) T=0"},
    10: {"L2": 0, "S2": 0, "J2": 0, "L1": 0, "J1": 1, "T2": 1, "desc": "1S0 T=3/2 (isolated)"},
}

# Measured costs at 96 threads (from the pilot)
MEASURED_96 = {
    (0, 0): 19.23, (1, 1): 288.03, (6, 6): 269.86, (8, 8): 64.31,
    (0, 1): 68.98, (0, 6): 68.46, (1, 6): 277.63, (6, 8): 675.21,
    (10, 10): 31.67,
}


def classify_block_cost(ar, ac):
    """Estimate block cost from channel properties using measured anchors."""
    if (ar, ac) in MEASURED_96:
        return MEASURED_96[(ar, ac)], "MEASURED"

    pa, pc = ALPHA_PROPS[ar], ALPHA_PROPS[ac]
    # Cost drivers: tensor coupling (S=1 with L mixing), high L, off-diagonal
    has_tensor = (pa["S2"] == 1 or pc["S2"] == 1)
    has_L_mismatch = (pa["L2"] != pc["L2"])
    has_high_L = (pa["L2"] >= 2 or pc["L2"] >= 2)
    is_diagonal = (ar == ac)

    # Cost tiers based on measured data:
    # TIER 1 (~20-35s): S-wave singlet diagonal, T=3/2
    # TIER 2 (~65-90s): D-wave diagonal, S-P off-diagonal
    # TIER 3 (~200-300s): P-wave/deuteron triplet diagonal, P-deuteron off-diagonal
    # TIER 4 (~500-700s): tensor D-wave coupling (L=1<->L=2, S=1)

    if is_diagonal:
        if pa["S2"] == 0 and pa["L2"] == 0:
            return 25.0, "ESTIMATED_TIER1"  # 1S0 diagonal
        elif pa["S2"] == 1 and pa["L2"] == 0:
            return 65.0, "ESTIMATED_TIER2"  # 3P1 diagonal (S-wave pair, J=1)
        elif pa["S2"] == 1 and pa["L2"] == 1 and pa["J2"] == 0:
            return 288.0, "ESTIMATED_TIER3"  # 3P0 diagonal
        elif pa["S2"] == 1 and pa["L2"] == 1 and pa["J2"] == 1:
            return 270.0, "ESTIMATED_TIER3"  # deuteron diagonal
        elif pa["L2"] == 2:
            return 65.0, "ESTIMATED_TIER2"  # 3D1 diagonal
        elif pa["S2"] == 0 and pa["L2"] == 1:
            return 70.0, "ESTIMATED_TIER2"  # 1P1 diagonal
        else:
            return 100.0, "ESTIMATED_DEFAULT"
    else:
        # Off-diagonal: estimate from coupling type
        if has_tensor and has_L_mismatch and has_high_L:
            return 675.0, "ESTIMATED_TIER4"  # tensor D-wave coupling
        elif has_tensor and pa["L2"] == 1 and pc["L2"] == 1:
            return 278.0, "ESTIMATED_TIER3"  # P-deuteron
        elif has_tensor and (pa["L2"] == 0 or pc["L2"] == 0):
            return 69.0, "ESTIMATED_TIER2"  # S-deuteron
        elif not has_tensor and has_L_mismatch:
            return 69.0, "ESTIMATED_TIER2"  # S-P transition
        elif not has_tensor and not has_L_mismatch:
            return 150.0, "ESTIMATED_TIER2_5"  # same-L off-diagonal
        else:
            return 200.0, "ESTIMATED_DEFAULT"


def get_full_build_costs():
    """If the full 56-block build completed, parse per-block costs from the log."""
    if not FULL_BUILD_LOG.exists():
        return None
    text = FULL_BUILD_LOG.read_text(errors='replace')
    costs = {}
    for line in text.splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if d.get("status") == "built" and "a_r" in d and "wall_seconds" in d:
            costs[(d["a_r"], d["a_c"])] = d["wall_seconds"]
    # Need all 56 evaluate blocks
    if len(costs) >= 50:  # allow a few missing
        return costs
    return None


def main():
    print("=" * 72)
    print("Phase C: static partition load imbalance assessment")
    print("=" * 72)

    # Enumerate all 56 evaluate blocks (a_r <= a_c, T_3N matching)
    # Alphas 0-9 have T_3N=1/2; alpha 10 has T_3N=3/2 (only self-couples).
    eval_blocks = []
    for ar in range(11):
        for ac in range(ar, 11):
            # Alpha 10 has T_3N=3/2; only (10,10) is allowed
            if (ar == 10) != (ac == 10):
                continue
            eval_blocks.append((ar, ac))
    print(f"Total evaluate blocks: {len(eval_blocks)}")

    # Try to use full build data; fall back to estimated
    full_costs = get_full_build_costs()
    if full_costs:
        print(f"Using FULL MEASURED costs from {len(full_costs)} blocks")
        source = "MEASURED_FULL_BUILD"
        block_costs = {}
        for b in eval_blocks:
            if b in full_costs:
                block_costs[b] = (full_costs[b], "MEASURED")
            else:
                est, label = classify_block_cost(*b)
                block_costs[b] = (est, "ESTIMATED_FALLBACK")
    else:
        print("Full build not complete; using ESTIMATED costs from 9 measured anchors")
        source = "ESTIMATED_FROM_PILOT"
        block_costs = {}
        for b in eval_blocks:
            est, label = classify_block_cost(*b)
            block_costs[b] = (est, label)

    # Print cost distribution
    costs = [c for c, _ in block_costs.values()]
    print(f"\nBlock cost distribution:")
    print(f"  min:    {min(costs):.1f} s")
    print(f"  max:    {max(costs):.1f} s")
    print(f"  mean:   {statistics.mean(costs):.1f} s")
    print(f"  median: {statistics.median(costs):.1f} s")
    print(f"  spread: {max(costs)/min(costs):.1f}x")

    # Simulate static modulo partition for N = 4, 8, 16, 32
    results = {}
    for N in [4, 8, 16, 32]:
        worker_loads = [0.0] * N
        worker_blocks = [[] for _ in range(N)]
        for idx, (ar, ac) in enumerate(eval_blocks):
            w = idx % N
            cost = block_costs[(ar, ac)][0]
            worker_loads[w] += cost
            worker_blocks[w].append((ar, ac, cost))

        mean_load = statistics.mean(worker_loads)
        max_load = max(worker_loads)
        min_load = min(worker_loads)
        eta = mean_load / max_load if max_load > 0 else 0
        imbalance_pct = (1 - eta) * 100

        results[N] = {
            "worker_loads_seconds": worker_loads,
            "mean_load_seconds": mean_load,
            "max_load_seconds": max_load,
            "min_load_seconds": min_load,
            "eta_load": eta,
            "imbalance_pct": imbalance_pct,
            "threshold_30pct_exceeded": (max_load / mean_load - 1) > 0.30,
            "max_worker_blocks": max(len(wb) for wb in worker_blocks),
            "min_worker_blocks": min(len(wb) for wb in worker_blocks),
        }

        print(f"\n--- N={N} workers ---")
        print(f"  worker loads (s): {[f'{l:.0f}' for l in worker_loads]}")
        print(f"  mean={mean_load:.0f}s max={max_load:.0f}s min={min_load:.0f}s")
        print(f"  eta_load = {eta:.3f}  (imbalance = {imbalance_pct:.1f}%)")
        verdict = "MATERIAL — cost-weighted scheduling recommended" \
            if results[N]["threshold_30pct_exceeded"] else "adequate"
        print(f"  threshold (max/mean > 1.30): {'EXCEEDED' if results[N]['threshold_30pct_exceeded'] else 'OK'}")
        print(f"  verdict: {verdict}")

    # Write JSON
    result = {
        "schema": "tictac.realistic_w1_load_imbalance.v1",
        "generated_on": "2026-08-15",
        "source": source,
        "n_evaluate_blocks": len(eval_blocks),
        "cost_distribution": {
            "min_seconds": min(costs),
            "max_seconds": max(costs),
            "mean_seconds": statistics.mean(costs),
            "median_seconds": statistics.median(costs),
            "spread_ratio": max(costs) / min(costs),
        },
        "static_modulo_partition": {
            str(N): {
                "eta_load": results[N]["eta_load"],
                "imbalance_pct": results[N]["imbalance_pct"],
                "mean_load_seconds": results[N]["mean_load_seconds"],
                "max_load_seconds": results[N]["max_load_seconds"],
                "min_load_seconds": results[N]["min_load_seconds"],
                "threshold_30pct_exceeded": results[N]["threshold_30pct_exceeded"],
            } for N in results
        },
        "block_costs": [
            {"a_r": ar, "a_c": ac, "cost_seconds": c, "source": s}
            for (ar, ac), (c, s) in sorted(block_costs.items())
        ],
        "conclusion": (
            "Static modulo partition has material load imbalance; "
            "cost-weighted scheduling recommended"
            if any(r["threshold_30pct_exceeded"] for r in results.values())
            else "Static modulo partition is adequate"
        ),
        "note": ("If source=ESTIMATED_FROM_PILOT, costs are estimated from 9 measured "
                 "pilot blocks using channel-property classification. The full 56-block "
                 "build was not yet complete when this analysis was run.")
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=2))
    print(f"\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
