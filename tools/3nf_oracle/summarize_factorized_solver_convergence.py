#!/usr/bin/env python3
"""Summarize complete-N2LO-3NF reduced-solver convergence runs.

The report keeps three distinct diagnostics separate:

* wave-packet grid resolution and the resulting deuteron binding energy;
* W1 radial-cell quadrature at identical grid/kinematic points;
* Padé honesty codes from the solver sidecar.

Different Chebyshev q grids generally produce different on-shell bin midpoints,
so amplitudes are compared only when every kinematic row matches exactly.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from compare_reduced_solver_outputs import (
    convergence_path,
    parse_convergence,
    parse_run_parameters,
    parse_u_rows,
    sha256,
)


def parse_run_spec(text: str) -> tuple[str, Path]:
    try:
        label, directory = text.split("=", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError("run must be LABEL=DIRECTORY") from error
    return label, Path(directory)


def parse_timing_spec(text: str) -> tuple[str, dict[str, float]]:
    try:
        label, w1_seconds, total_seconds = text.split(",", 2)
        return label, {
            "w1_build_seconds": float(w1_seconds),
            "total_run_seconds": float(total_seconds),
        }
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "timing must be LABEL,W1_BUILD_SECONDS,TOTAL_RUN_SECONDS"
        ) from error


def parse_compare_spec(text: str) -> tuple[str, str]:
    try:
        left, right = text.split(",", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError("comparison must be LEFT_LABEL,RIGHT_LABEL") from error
    return left, right


def summarize_run(label: str, directory: Path, timing: dict[str, float] | None) -> dict:
    parameters = parse_run_parameters(directory)
    u_files = sorted(directory.glob("U_PW_elements_*.txt"))
    if len(u_files) != 1:
        raise ValueError(f"expected one U file under {directory}, found {len(u_files)}")
    u_path = u_files[0]
    binding, rows = parse_u_rows(u_path)
    conv_path = convergence_path(u_path)
    convergence = parse_convergence(conv_path)
    codes = [record["code"] for record in convergence]
    result = {
        "label": label,
        "directory": str(directory),
        "grid": {
            "Np_WP": int(parameters["Np_WP"]),
            "Nq_WP": int(parameters["Nq_WP"]),
            "J_2N_max": int(parameters["J_2N_max"]),
            "two_J_3N_max": int(parameters["two_J_3N_max"]),
            "Np_per_WP_W1": int(parameters["Np_per_WP_W1"]),
            "Nq_per_WP_W1": int(parameters["Nq_per_WP_W1"]),
            "Nangle_3NF": int(parameters["Nangle_3NF"]),
        },
        "deuteron_binding_mev": binding,
        "u_sha256": sha256(u_path),
        "convergence_sha256": sha256(conv_path),
        "num_kinematic_rows": len(rows),
        "num_complex_amplitudes": sum(len(row[3]) for row in rows),
        "kinematics": [
            {"Tlab_mev": row[0], "Ecm_mev": row[1], "q_index": row[2]}
            for row in rows
        ],
        "pade_status": {
            "num_truly_converged": codes.count(1),
            "num_max_order_truncated": codes.count(2),
            "all_truly_converged": bool(codes) and all(code == 1 for code in codes),
        },
        "_rows": rows,
    }
    if timing is not None:
        result["timing"] = timing
    return result


def compare_runs(left: dict, right: dict) -> dict:
    left_rows = left["_rows"]
    right_rows = right["_rows"]
    same_kinematics = (
        len(left_rows) == len(right_rows)
        and all(left_row[:3] == right_row[:3] for left_row, right_row in zip(left_rows, right_rows))
    )
    result = {
        "left": left["label"],
        "right": right["label"],
        "same_kinematics": same_kinematics,
    }
    if not same_kinematics:
        result["amplitude_comparison"] = "not_performed"
        return result

    differences = []
    for left_row, right_row in zip(left_rows, right_rows):
        if len(left_row[3]) != len(right_row[3]):
            raise ValueError(f"amplitude-count mismatch: {left['label']} vs {right['label']}")
        for left_value, right_value in zip(left_row[3], right_row[3]):
            absolute = abs(left_value - right_value)
            relative = absolute / max(abs(right_value), 1.0e-300)
            differences.append((absolute, relative))
    result.update({
        "num_complex_amplitudes": len(differences),
        "max_absolute_difference_mev": max(item[0] for item in differences),
        "max_relative_difference": max(item[1] for item in differences),
        "rms_absolute_difference_mev": math.sqrt(
            sum(item[0] ** 2 for item in differences) / len(differences)
        ),
    })
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", type=parse_run_spec, required=True)
    parser.add_argument("--timing", action="append", type=parse_timing_spec, default=[])
    parser.add_argument("--compare", action="append", type=parse_compare_spec, default=[])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    timings = dict(args.timing)
    runs = [summarize_run(label, directory, timings.get(label)) for label, directory in args.run]
    by_label = {run["label"]: run for run in runs}
    comparisons = []
    for left_label, right_label in args.compare:
        comparisons.append(compare_runs(by_label[left_label], by_label[right_label]))
    for run in runs:
        run.pop("_rows")

    report = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat()
                         .replace("+00:00", "Z"),
        "purpose": "reduced complete-N2LO-3NF WPCD convergence diagnostics",
        "runs": runs,
        "comparisons": comparisons,
        "limitations": [
            "J_2N_max=1 and two_J_3N_max=1 are reduced diagnostic truncations.",
            "Nangle_3NF=2 is not an angular-convergence certificate.",
            "Radial order one is the unconverged midpoint diagnostic.",
            "Chebyshev q grids at different Nq have different on-shell bin midpoints; their amplitudes are not compared directly.",
            "Padé maximum-order truncation codes remain a separate failed convergence gate.",
        ],
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
