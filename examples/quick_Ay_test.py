#!/usr/bin/env python3
"""
Purpose:
  Fast smoke test for the maintained Tlab benchmark pipeline.

Data flow (left -> right):
  CLI args
    -> run `deuteron_proton_Ay.py`
    -> run `compare_Ay_experiment.py`
    -> optionally run `plot_validation_curves.py`
    -> print return codes and summary JSON fields

Calls / dependencies:
  - External process wrapper around scripts in `examples/`.
  - Reads `solver_validation_tlab_*.json` when present.

Usage:
  python3 examples/quick_Ay_test.py --work-dir output/quick_Ay_test --plot
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import List

from tlab_utils import format_tlab_label


def run_command(cmd: List[str], cwd: Path, timeout_s: int) -> int:
    print(f"[RUN] {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            text=True,
            capture_output=True,
            timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        print(f"[FAIL] Timeout after {timeout_s} s")
        return 124

    if result.stdout.strip():
        print("[STDOUT tail]")
        print("\n".join(result.stdout.strip().splitlines()[-20:]))
    if result.stderr.strip():
        print("[STDERR tail]")
        print("\n".join(result.stderr.strip().splitlines()[-20:]))

    print(f"[RC] {result.returncode}")
    return result.returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Quick Tic-tac Tlab smoke validation")
    parser.add_argument("--work-dir", default="output/quick_Ay_test", help="Working directory")
    parser.add_argument("--target-tlab-mev", type=float, default=190.0)
    parser.add_argument("--two-j-3n-max", type=int, default=1)
    parser.add_argument("--j-2n-max", type=int, default=2)
    parser.add_argument("--np", type=int, default=20)
    parser.add_argument("--nq", type=int, default=20)
    parser.add_argument("--nphi", type=int, default=24)
    parser.add_argument("--nx", type=int, default=24)
    parser.add_argument("--np-per-wp", type=int, default=6)
    parser.add_argument("--nq-per-wp", type=int, default=6)
    parser.add_argument("--chebyshev-s", type=float, default=150.0)
    parser.add_argument("--chebyshev-t", type=float, default=1.0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--potential-model", default="LO_internal")
    parser.add_argument("--timeout-solver", type=int, default=3600)
    parser.add_argument("--reuse-p123", action="store_true")
    parser.add_argument("--plot", action="store_true", help="Also run plotting script")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    work_dir = (root / args.work_dir).resolve()
    python = sys.executable

    # Step 1: produce solver outputs for the requested settings.
    run_solver_cmd = [
        python,
        str(root / "examples" / "deuteron_proton_Ay.py"),
        "--work-dir",
        str(work_dir),
        "--target-tlab-mev",
        str(args.target_tlab_mev),
        "--two-j-3n-max",
        str(args.two_j_3n_max),
        "--j-2n-max",
        str(args.j_2n_max),
        "--np",
        str(args.np),
        "--nq",
        str(args.nq),
        "--nphi",
        str(args.nphi),
        "--nx",
        str(args.nx),
        "--np-per-wp",
        str(args.np_per_wp),
        "--nq-per-wp",
        str(args.nq_per_wp),
        "--chebyshev-s",
        str(args.chebyshev_s),
        "--chebyshev-t",
        str(args.chebyshev_t),
        "--threads",
        str(args.threads),
        "--potential-model",
        args.potential_model,
        "--timeout",
        str(args.timeout_solver),
    ]
    if args.reuse_p123:
        run_solver_cmd.append("--reuse-p123")

    rc = run_command(run_solver_cmd, root, timeout_s=args.timeout_solver + 120)
    if rc != 0:
        return rc

    # Step 2: compare model outputs against experiment and write validation artifacts.
    compare_cmd = [
        python,
        str(root / "examples" / "compare_Ay_experiment.py"),
        "--work-dir",
        str(work_dir),
        "--solver-out-dir",
        str(work_dir / "solver_out"),
        "--target-tlab-mev",
        str(args.target_tlab_mev),
    ]
    rc = run_command(compare_cmd, root, timeout_s=180)

    # Step 3: optional lightweight summary extraction for CI/smoke visibility.
    summary_label = format_tlab_label(args.target_tlab_mev)
    summary_path = work_dir / f"solver_validation_tlab_{summary_label}.json"
    if summary_path.exists():
        try:
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            best = payload.get("best_energy", {})
            print("[SUMMARY]")
            print(
                f"best_tlab_mev={best.get('tlab_mev', 'NA')} "
                f"target_tlab_mev={payload.get('target_tlab_mev', 'NA')} "
                f"abs_delta_tlab_mev={best.get('abs_delta_tlab_mev', 'NA')}"
            )
            print(f"status={payload.get('status', 'NA')}")
        except Exception as exc:
            print(f"[WARN] Failed to parse summary JSON: {exc}")

    if args.plot:
        plot_cmd = [
            python,
            str(root / "examples" / "plot_validation_curves.py"),
            "--work-dir",
            str(work_dir),
            "--summary-json",
            summary_path.name,
        ]
        _ = run_command(plot_cmd, root, timeout_s=180)

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
