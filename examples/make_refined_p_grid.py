#!/usr/bin/env python3
"""Generate a Chebyshev p-grid locally refined in the P-wave sensitive band.

The default Tic-tac p-grid uses a Chebyshev distribution
    b_i = s * tan((2i-1) pi / (4 N))^t,   b_0 = 0
which concentrates points near p=0 but is coarse around p ~ 50-200 MeV/c —
the band that drives the doublet ²P phase shifts (Ay puzzle). This script
takes the same Chebyshev seed and subdivides every bin whose midpoint falls
in [p_min_refine, p_max_refine] into `refine_factor` equal pieces. The
resulting boundary list is consumed by Tic-tac via:
    p_grid_type=custom
    p_grid_filename=<output path>

Output is one boundary per line (Np+1 values), the format expected by
read_WP_boundaries_from_txt. b_0 = 0 is preserved.

Example:
    python examples/make_refined_p_grid.py \\
        --Np 20 --chebyshev-s 100 --chebyshev-t 1 \\
        --p-min-refine 30 --p-max-refine 200 --refine-factor 2 \\
        --out CPP/Input/p_grid_np20_refined.txt
"""

import argparse
import math
from pathlib import Path


def chebyshev_boundaries(N: int, s: float, t: float) -> list[float]:
    """Reproduce make_chebyshev_distribution() from src/core/state_space/make_wp_states.cpp."""
    boundaries = [0.0]
    for i in range(1, N + 1):
        tan_term = math.tan((2 * i - 1) * math.pi / (4 * N))
        boundaries.append(s * tan_term**t)
    return boundaries


def refine_band(boundaries: list[float],
                p_min: float,
                p_max: float,
                factor: int) -> list[float]:
    """Subdivide every bin whose midpoint is in [p_min, p_max] into `factor` pieces."""
    if factor < 1:
        raise ValueError("refine_factor must be >= 1")
    if factor == 1:
        return list(boundaries)

    refined: list[float] = []
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        refined.append(lo)
        mid = 0.5 * (lo + hi)
        if p_min <= mid <= p_max:
            for k in range(1, factor):
                refined.append(lo + (hi - lo) * k / factor)
    refined.append(boundaries[-1])
    return refined


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Np", type=int, required=True,
                    help="Chebyshev seed bin count (matches Np_WP in input.txt)")
    ap.add_argument("--chebyshev-s", type=float, default=100.0,
                    help="Chebyshev scale s (MeV)")
    ap.add_argument("--chebyshev-t", type=float, default=1.0,
                    help="Chebyshev sparseness exponent t")
    ap.add_argument("--p-min-refine", type=float, default=30.0,
                    help="Lower edge of P-wave sensitive band (MeV)")
    ap.add_argument("--p-max-refine", type=float, default=200.0,
                    help="Upper edge of P-wave sensitive band (MeV)")
    ap.add_argument("--refine-factor", type=int, default=2,
                    help="Subdivision factor inside the band (k bins -> k*factor bins)")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output boundary file (one value per line, Np+1 entries)")
    args = ap.parse_args()

    if args.p_min_refine >= args.p_max_refine:
        raise SystemExit("p_min_refine must be < p_max_refine")

    seed = chebyshev_boundaries(args.Np, args.chebyshev_s, args.chebyshev_t)
    refined = refine_band(seed, args.p_min_refine, args.p_max_refine, args.refine_factor)

    n_added = len(refined) - len(seed)
    n_bins = len(refined) - 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        for b in refined:
            f.write(f"{b:.12g}\n")

    print(f"seed Chebyshev grid:    Np={args.Np} (s={args.chebyshev_s}, t={args.chebyshev_t})")
    print(f"refine band [MeV]:      [{args.p_min_refine}, {args.p_max_refine}] x{args.refine_factor}")
    print(f"refined bins added:     {n_added}")
    print(f"final Np_WP:            {n_bins}")
    print(f"wrote:                  {args.out}")
    print()
    print("Use in input.txt:")
    print(f"  Np_WP={n_bins}")
    print( "  p_grid_type=custom")
    print(f"  p_grid_filename={args.out}")


if __name__ == "__main__":
    main()
