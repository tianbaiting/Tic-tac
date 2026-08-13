#!/usr/bin/env python3
"""Fail-closed audit and comparison for the low-energy nd Ay benchmark.

The script computes paired 2NF and 2NF+3NF Ay curves, but calls the result an
acceptance benchmark only when the stored solver artifacts prove that the two
runs differ solely by the 3NF switch, use the complete factorized N2LO force,
have a physical deuteron, cover every requested J^pi block, contain no
max-order-truncated Pade amplitudes, and include converged Ay ladders for every
required numerical dimension.

Convergence references are supplied as
``--convergence-pair DIMENSION=TWO_NF_DIR,THREE_NF_DIR``.  The required
dimensions are Np, Nq, W1, J2N, J3N, and Pade.  Missing ladders remain explicit
failed gates; no command-line assertion can turn them into evidence.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable, Sequence

import numpy as np


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "examples"))

import pw_amplitudes as pw  # noqa: E402


COMPLETE_3NF_MODEL = "chiral_N2LO_full_factorized"
NO_3NF_MODELS = {"none", "no_3nf", "disabled"}
MEASURED_DEUTERON_BINDING_MEV = 2.22457
REQUIRED_LADDERS = ("Np", "Nq", "W1", "J2N", "J3N", "Pade")


PARAMETER_LABELS: dict[str, tuple[str, Any]] = {
    "two_J_3N_max": ("two_j_3n_max", int),
    "Np_WP": ("np_wp", int),
    "Nq_WP": ("nq_wp", int),
    "J_2N_max": ("j_2n_max", int),
    "Nphi": ("nphi", int),
    "Nx": ("nx", int),
    "chebyshev sparseness": ("chebyshev_t", float),
    "chebyshev scale": ("chebyshev_s", float),
    "Np_per_WP": ("np_per_wp", int),
    "Nq_per_WP": ("nq_per_wp", int),
    "Np_per_WP_W1": ("np_per_wp_w1", int),
    "Nq_per_WP_W1": ("nq_per_wp_w1", int),
    "Nangle_3NF": ("nangle_3nf", int),
    "Padé maximum diagonal order": ("pade_max_order", int),
    "Tensor-force on": ("tensor_force", str),
    "Isospin-breaking in 1S0": ("isospin_breaking_1s0", str),
    "Mid-point approximation": ("midpoint_approx", str),
    "Calculate breakup amplitudes": ("include_breakup_channels", str),
    "Potential model": ("potential_model", str),
    "Three-nucleon force": ("three_nucleon_force", str),
    "3NF LEC c_D": ("c_d", float),
    "3NF LEC c_E": ("c_e", float),
    "3NF cutoff Lambda_3NF [MeV]": ("lambda_3nf_mev", float),
    "p-momentum grid type": ("p_grid_type", str),
    "p-momentum grid input file": ("p_grid_file", str),
    "q-momentum grid type": ("q_grid_type", str),
    "q-momentum grid input file": ("q_grid_file", str),
    "Energy input file": ("energy_input_file", str),
}

PAIR_MATCH_KEYS = (
    "two_j_3n_max",
    "np_wp",
    "nq_wp",
    "j_2n_max",
    "nphi",
    "nx",
    "chebyshev_t",
    "chebyshev_s",
    "np_per_wp",
    "nq_per_wp",
    "np_per_wp_w1",
    "nq_per_wp_w1",
    "pade_max_order",
    "tensor_force",
    "isospin_breaking_1s0",
    "midpoint_approx",
    "include_breakup_channels",
    "potential_model",
    "p_grid_type",
    "p_grid_file",
    "q_grid_type",
    "q_grid_file",
    "energy_input_file",
)

LADDER_ALLOWED_KEYS = {
    "Np": {"np_wp"},
    "Nq": {"nq_wp"},
    "W1": {"np_per_wp_w1", "nq_per_wp_w1"},
    "J2N": {"j_2n_max"},
    "J3N": {"two_j_3n_max"},
    "Pade": {"pade_max_order"},
}
LADDER_COMPARISON_KEYS = PAIR_MATCH_KEYS + (
    "three_nucleon_force",
    "c_d",
    "c_e",
    "lambda_3nf_mev",
    "nangle_3nf",
)


@dataclass(frozen=True)
class CurveResult:
    values: np.ndarray
    lower_tlab_mev: float
    upper_tlab_mev: float
    upper_weight: float


@dataclass(frozen=True)
class RunAudit:
    run_dir: str
    parameters: dict[str, Any]
    binding_energy_mev: float | None
    binding_headers_consistent: bool
    expected_blocks: list[str]
    present_blocks: list[str]
    missing_blocks: list[str]
    pade_counts: dict[str, int]


def _convert_value(raw: str, converter: Any) -> Any:
    value = raw.strip()
    if converter is str:
        return value
    return converter(value.split()[0])


def parse_run_parameters(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "run_parameters.txt"
    if not path.is_file():
        raise FileNotFoundError(f"missing solver metadata: {path}")
    parameters: dict[str, Any] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in raw_line:
            continue
        label, raw_value = raw_line.split(":", 1)
        specification = PARAMETER_LABELS.get(label.strip())
        if specification is None:
            continue
        key, converter = specification
        parameters[key] = _convert_value(raw_value, converter)
    return parameters


_BINDING_PATTERN = re.compile(
    r"^# Deuteron BE:\s*([+\-0-9.eE]+)\s+MeV\s*$", re.MULTILINE
)
_BLOCK_PATTERN = re.compile(r"_JP_(\d+)_(-?1)_Jmax_")


def binding_energy_from_outputs(run_dir: Path) -> tuple[float | None, bool]:
    values: list[float] = []
    for path in sorted(run_dir.glob("U_PW_elements_*.txt")):
        match = _BINDING_PATTERN.search(path.read_text(encoding="utf-8"))
        if match:
            values.append(float(match.group(1)))
    if not values:
        return None, False
    consistent = max(values) - min(values) <= 1.0e-7
    return float(sum(values) / len(values)), consistent


def present_blocks(run_dir: Path) -> set[tuple[int, int]]:
    result: set[tuple[int, int]] = set()
    for path in run_dir.glob("U_PW_elements_*.txt"):
        match = _BLOCK_PATTERN.search(path.name)
        if match:
            result.add((int(match.group(1)), int(match.group(2))))
    return result


def expected_blocks(two_j_3n_max: int) -> set[tuple[int, int]]:
    return {
        (two_j, parity)
        for two_j in range(1, two_j_3n_max + 1, 2)
        for parity in (-1, 1)
    }


def format_block(block: tuple[int, int]) -> str:
    two_j, parity = block
    return f"{two_j}/2{'+' if parity == 1 else '-'}"


def pade_counts(run_dir: Path) -> dict[str, int]:
    counts = {"converged": 0, "truncated": 0, "unknown": 0}
    for path in sorted(run_dir.glob("U_PW_convergence_*.txt")):
        for raw_line in path.read_text(encoding="utf-8").splitlines():
            fields = raw_line.split()
            if not fields or fields[0].startswith("#") or len(fields) < 5:
                continue
            status = int(fields[3])
            if status == 1:
                counts["converged"] += 1
            elif status == 2:
                counts["truncated"] += 1
            else:
                counts["unknown"] += 1
    return counts


def inspect_run(run_dir: Path) -> RunAudit:
    run_dir = run_dir.resolve()
    parameters = parse_run_parameters(run_dir)
    maximum_j = parameters.get("two_j_3n_max")
    expected = expected_blocks(maximum_j) if isinstance(maximum_j, int) else set()
    present = present_blocks(run_dir)
    binding, consistent = binding_energy_from_outputs(run_dir)
    return RunAudit(
        run_dir=display_path(run_dir),
        parameters=parameters,
        binding_energy_mev=binding,
        binding_headers_consistent=consistent,
        expected_blocks=[format_block(block) for block in sorted(expected)],
        present_blocks=[format_block(block) for block in sorted(present)],
        missing_blocks=[format_block(block) for block in sorted(expected - present)],
        pade_counts=pade_counts(run_dir),
    )


def paired_parameter_mismatches(
    two_nf: RunAudit, three_nf: RunAudit
) -> dict[str, dict[str, Any]]:
    mismatches: dict[str, dict[str, Any]] = {}
    for key in PAIR_MATCH_KEYS:
        left = two_nf.parameters.get(key)
        right = three_nf.parameters.get(key)
        if left != right:
            mismatches[key] = {"two_nf": left, "two_nf_plus_3nf": right}
    return mismatches


def ladder_parameter_changes(
    main: RunAudit, reference: RunAudit
) -> dict[str, dict[str, Any]]:
    changes: dict[str, dict[str, Any]] = {}
    for key in LADDER_COMPARISON_KEYS:
        main_value = main.parameters.get(key)
        reference_value = reference.parameters.get(key)
        if main_value != reference_value:
            changes[key] = {"main": main_value, "reference": reference_value}
    return changes


def valid_ladder_change(changes: dict[str, Any], dimension: str) -> bool:
    changed_keys = set(changes)
    return bool(changed_keys) and changed_keys.issubset(
        LADDER_ALLOWED_KEYS[dimension]
    )


def energy_to_q_index(
    blocks: Sequence[pw.JPiBlock], tlab_mev: float, ecm_mev: float
) -> int:
    point = min(
        blocks[0].points,
        key=lambda item: abs(item.tlab - tlab_mev) + abs(item.ecm - ecm_mev),
    )
    if abs(point.tlab - tlab_mev) > 1.0e-6 or abs(point.ecm - ecm_mev) > 1.0e-6:
        raise ValueError(
            f"could not map Tlab={tlab_mev}, Ecm={ecm_mev} to a q bin"
        )
    return int(point.q_idx)


def interpolate_curve(
    run_dir: Path, target_tlab_mev: float, angles_deg: Sequence[float]
) -> CurveResult:
    blocks = pw.parse_solver_output(run_dir)
    q_kinematics = pw.parse_q_kinematics(run_dir)
    energies = sorted(pw.list_solver_energies(blocks))
    lower = [entry for entry in energies if entry[0] <= target_tlab_mev]
    upper = [entry for entry in energies if entry[0] >= target_tlab_mev]
    if not lower or not upper:
        available = ", ".join(f"{entry[0]:.6g}" for entry in energies)
        raise ValueError(
            f"Tlab={target_tlab_mev} MeV is not bracketed in {run_dir}; "
            f"available energies: {available}"
        )
    lo = lower[-1]
    hi = upper[0]
    q_lo = energy_to_q_index(blocks, *lo)
    q_hi = energy_to_q_index(blocks, *hi)
    low_values = np.asarray([
        item.Ay_n
        for item in pw.observables_at_angles(
            blocks, q_lo, angles_deg, bin_info=q_kinematics[q_lo]
        )
    ])
    if q_lo == q_hi:
        weight = 0.0
        values = low_values
    else:
        high_values = np.asarray([
            item.Ay_n
            for item in pw.observables_at_angles(
                blocks, q_hi, angles_deg, bin_info=q_kinematics[q_hi]
            )
        ])
        weight = (target_tlab_mev - lo[0]) / (hi[0] - lo[0])
        values = (1.0 - weight) * low_values + weight * high_values
    return CurveResult(values, float(lo[0]), float(hi[0]), float(weight))


def load_experiment(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    angles: list[float] = []
    values: list[float] = []
    uncertainties: list[float] = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if not row or row[0].startswith("#") or row[0] == "theta_cm_deg":
                continue
            angles.append(float(row[0]))
            values.append(float(row[1]))
            uncertainties.append(float(row[2]))
    if not angles:
        raise ValueError(f"no experimental data in {path}")
    return np.asarray(angles), np.asarray(values), np.asarray(uncertainties)


def curve_metrics(
    theta_grid: np.ndarray,
    curve: np.ndarray,
    experiment_theta: np.ndarray,
    experiment_values: np.ndarray,
    experiment_uncertainties: np.ndarray,
) -> dict[str, Any]:
    predictions = np.interp(experiment_theta, theta_grid, curve)
    residuals = predictions - experiment_values
    peak_index = int(np.argmax(curve))
    experiment_peak_index = int(np.argmax(experiment_values))
    valid_uncertainties = experiment_uncertainties > 0.0
    chi2 = float(np.sum(np.square(
        residuals[valid_uncertainties] / experiment_uncertainties[valid_uncertainties]
    )))
    n_chi2 = int(np.count_nonzero(valid_uncertainties))
    return {
        "ay_max": float(curve[peak_index]),
        "theta_at_ay_max_deg": float(theta_grid[peak_index]),
        "experimental_ay_max": float(experiment_values[experiment_peak_index]),
        "experimental_peak_theta_deg": float(experiment_theta[experiment_peak_index]),
        "peak_deficit_model_minus_experiment": float(
            curve[peak_index] - experiment_values[experiment_peak_index]
        ),
        "prediction_at_experimental_peak": float(predictions[experiment_peak_index]),
        "deficit_at_experimental_peak": float(
            predictions[experiment_peak_index] - experiment_values[experiment_peak_index]
        ),
        "mae": float(np.mean(np.abs(residuals))),
        "rmse": float(math.sqrt(np.mean(np.square(residuals)))),
        "maximum_absolute_residual": float(np.max(np.abs(residuals))),
        "chi2": chi2,
        "chi2_points": n_chi2,
        "chi2_per_point": chi2 / n_chi2 if n_chi2 else None,
        "predictions_at_experiment": predictions.tolist(),
        "residuals_at_experiment": residuals.tolist(),
    }


def core_gates(
    two_nf: RunAudit,
    three_nf: RunAudit,
    binding_tolerance_mev: float,
) -> dict[str, dict[str, Any]]:
    mismatches = paired_parameter_mismatches(two_nf, three_nf)
    binding_values = (two_nf.binding_energy_mev, three_nf.binding_energy_mev)
    binding_present = all(value is not None for value in binding_values)
    binding_error = (
        max(abs(abs(float(value)) - MEASURED_DEUTERON_BINDING_MEV)
            for value in binding_values)
        if binding_present else None
    )
    binding_pair_delta = (
        abs(float(binding_values[0]) - float(binding_values[1]))
        if binding_present else None
    )

    def pade_passes(audit: RunAudit) -> bool:
        counts = audit.pade_counts
        return (
            counts["converged"] > 0
            and counts["truncated"] == 0
            and counts["unknown"] == 0
        )

    return {
        "paired_numerical_setup": {
            "passed": not mismatches,
            "mismatches": mismatches,
        },
        "model_identity": {
            "passed": (
                str(two_nf.parameters.get("three_nucleon_force", "")).lower()
                in NO_3NF_MODELS
                and three_nf.parameters.get("three_nucleon_force")
                == COMPLETE_3NF_MODEL
            ),
            "two_nf": two_nf.parameters.get("three_nucleon_force"),
            "two_nf_plus_3nf": three_nf.parameters.get("three_nucleon_force"),
            "required_three_nf": COMPLETE_3NF_MODEL,
        },
        "deuteron_binding": {
            "passed": bool(
                binding_present
                and two_nf.binding_headers_consistent
                and three_nf.binding_headers_consistent
                and binding_error is not None
                and binding_error <= binding_tolerance_mev
                and binding_pair_delta is not None
                and binding_pair_delta <= 1.0e-7
            ),
            "two_nf_mev": binding_values[0],
            "two_nf_plus_3nf_mev": binding_values[1],
            "max_error_from_measured_mev": binding_error,
            "pair_delta_mev": binding_pair_delta,
            "tolerance_mev": binding_tolerance_mev,
        },
        "partial_wave_block_completeness": {
            "passed": not two_nf.missing_blocks and not three_nf.missing_blocks,
            "two_nf_missing": two_nf.missing_blocks,
            "two_nf_plus_3nf_missing": three_nf.missing_blocks,
        },
        "pade_honesty": {
            "passed": pade_passes(two_nf) and pade_passes(three_nf),
            "two_nf": two_nf.pade_counts,
            "two_nf_plus_3nf": three_nf.pade_counts,
        },
    }


def parse_convergence_pair(raw: str) -> tuple[str, Path, Path]:
    if "=" not in raw or "," not in raw:
        raise ValueError(
            "convergence pair must be DIMENSION=TWO_NF_DIR,THREE_NF_DIR"
        )
    dimension, paths = raw.split("=", 1)
    two_nf, three_nf = paths.split(",", 1)
    dimension = dimension.strip()
    if dimension not in REQUIRED_LADDERS:
        raise ValueError(
            f"unknown convergence dimension {dimension!r}; "
            f"choose one of {', '.join(REQUIRED_LADDERS)}"
        )
    return dimension, Path(two_nf.strip()), Path(three_nf.strip())


def path_from_repo(path: Path) -> Path:
    return path if path.is_absolute() else REPO / path


def display_path(path: Path) -> str:
    """Keep repository artifacts portable while preserving external paths."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO))
    except ValueError:
        return str(resolved)


def write_curve_csv(
    path: Path, theta: np.ndarray, two_nf: np.ndarray, three_nf: np.ndarray
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["theta_cm_deg", "Ay_2NF", "Ay_2NF_plus_3NF", "delta_Ay"])
        for angle, left, right in zip(theta, two_nf, three_nf):
            writer.writerow([
                f"{angle:.6f}", f"{left:.12e}", f"{right:.12e}",
                f"{right - left:.12e}",
            ])


def build_report(args: argparse.Namespace) -> tuple[dict[str, Any], np.ndarray, np.ndarray, np.ndarray]:
    two_nf_dir = path_from_repo(args.two_nf_dir).resolve()
    three_nf_dir = path_from_repo(args.three_nf_dir).resolve()
    experiment_path = path_from_repo(args.experiment).resolve()
    theta = np.arange(0.0, 181.0, 1.0)
    experiment_theta, experiment_values, experiment_uncertainties = load_experiment(
        experiment_path
    )

    two_audit = inspect_run(two_nf_dir)
    three_audit = inspect_run(three_nf_dir)
    two_curve = interpolate_curve(two_nf_dir, args.target_tlab_mev, theta)
    three_curve = interpolate_curve(three_nf_dir, args.target_tlab_mev, theta)
    same_bracket = (
        abs(two_curve.lower_tlab_mev - three_curve.lower_tlab_mev) <= 1.0e-9
        and abs(two_curve.upper_tlab_mev - three_curve.upper_tlab_mev) <= 1.0e-9
    )

    gates = core_gates(two_audit, three_audit, args.binding_tolerance_mev)
    gates["matched_energy_treatment"] = {
        "passed": same_bracket,
        "two_nf_bracket_mev": [
            two_curve.lower_tlab_mev, two_curve.upper_tlab_mev,
        ],
        "two_nf_plus_3nf_bracket_mev": [
            three_curve.lower_tlab_mev, three_curve.upper_tlab_mev,
        ],
    }

    convergence: dict[str, Any] = {}
    for raw in args.convergence_pair:
        dimension, reference_two_dir, reference_three_dir = parse_convergence_pair(raw)
        if dimension in convergence:
            raise ValueError(f"duplicate convergence dimension {dimension}")
        reference_two_dir = path_from_repo(reference_two_dir).resolve()
        reference_three_dir = path_from_repo(reference_three_dir).resolve()
        reference_two_audit = inspect_run(reference_two_dir)
        reference_three_audit = inspect_run(reference_three_dir)
        reference_gates = core_gates(
            reference_two_audit, reference_three_audit,
            args.binding_tolerance_mev,
        )
        two_changes = ladder_parameter_changes(two_audit, reference_two_audit)
        three_changes = ladder_parameter_changes(three_audit, reference_three_audit)
        allowed_keys = LADDER_ALLOWED_KEYS[dimension]

        structure_checks = {
            "reference_pair_passes_core_gates": all(
                item["passed"] for item in reference_gates.values()
            ),
            "two_nf_changes_only_requested_dimension": valid_ladder_change(
                two_changes, dimension),
            "two_nf_plus_3nf_changes_only_requested_dimension": (
                valid_ladder_change(three_changes, dimension)
            ),
        }
        reference_two = interpolate_curve(
            reference_two_dir, args.target_tlab_mev, theta
        ).values
        reference_three = interpolate_curve(
            reference_three_dir, args.target_tlab_mev, theta
        ).values
        differences = {
            "two_nf_max_abs_delta_Ay": float(
                np.max(np.abs(two_curve.values - reference_two))
            ),
            "two_nf_plus_3nf_max_abs_delta_Ay": float(
                np.max(np.abs(three_curve.values - reference_three))
            ),
            "three_nf_effect_max_abs_delta_Ay": float(np.max(np.abs(
                (three_curve.values - two_curve.values)
                - (reference_three - reference_two)
            ))),
        }
        convergence[dimension] = {
            "reference_two_nf_dir": str(reference_two_dir),
            "reference_two_nf_plus_3nf_dir": str(reference_three_dir),
            "allowed_parameter_changes": sorted(allowed_keys),
            "two_nf_parameter_changes": two_changes,
            "two_nf_plus_3nf_parameter_changes": three_changes,
            "reference_core_gates": reference_gates,
            "structure_checks": structure_checks,
            **differences,
            "tolerance_Ay": args.ay_convergence_tolerance,
            "passed": (
                all(structure_checks.values())
                and max(differences.values()) <= args.ay_convergence_tolerance
            ),
        }

    missing_ladders = [name for name in REQUIRED_LADDERS if name not in convergence]
    ladder_passed = not missing_ladders and all(
        item["passed"] for item in convergence.values()
    )
    gates["Ay_convergence_ladders"] = {
        "passed": ladder_passed,
        "required": list(REQUIRED_LADDERS),
        "missing": missing_ladders,
        "results": convergence,
    }

    two_metrics = curve_metrics(
        theta, two_curve.values, experiment_theta, experiment_values,
        experiment_uncertainties,
    )
    three_metrics = curve_metrics(
        theta, three_curve.values, experiment_theta, experiment_values,
        experiment_uncertainties,
    )
    effect = three_curve.values - two_curve.values
    effect_max = float(np.max(np.abs(effect)))
    convergence_errors = [
        value
        for item in convergence.values()
        for key, value in item.items()
        if key.endswith("max_abs_delta_Ay")
    ]
    numerical_uncertainty = max(convergence_errors) if convergence_errors else None
    gates["effect_resolved_above_numerical_uncertainty"] = {
        "passed": bool(
            numerical_uncertainty is not None
            and numerical_uncertainty < effect_max
        ),
        "max_abs_three_nf_effect_Ay": effect_max,
        "estimated_numerical_uncertainty_Ay": numerical_uncertainty,
    }

    accepted = all(item["passed"] for item in gates.values())
    report = {
        "schema": "tictac.low_energy_nd_Ay_validation.v1",
        "target_tlab_mev": args.target_tlab_mev,
        "accepted_physics_benchmark": accepted,
        "disposition": (
            "converged complete-N2LO benchmark"
            if accepted else "diagnostic only; failed acceptance gates"
        ),
        "experiment": {
            "path": display_path(experiment_path),
            "theta_deg": experiment_theta.tolist(),
            "Ay": experiment_values.tolist(),
            "uncertainty": experiment_uncertainties.tolist(),
        },
        "runs": {
            "two_nf": asdict(two_audit),
            "two_nf_plus_3nf": asdict(three_audit),
        },
        "gates": gates,
        "energy_interpolation": {
            "two_nf": asdict(two_curve) | {"values": None},
            "two_nf_plus_3nf": asdict(three_curve) | {"values": None},
        },
        "metrics": {
            "two_nf": two_metrics,
            "two_nf_plus_3nf": three_metrics,
            "R_Ay_rmse_ratio": (
                three_metrics["rmse"] / two_metrics["rmse"]
                if two_metrics["rmse"] > 0.0 else None
            ),
            "three_nf_effect": {
                "max_abs_delta_Ay": effect_max,
                "theta_at_max_abs_delta_deg": float(
                    theta[int(np.argmax(np.abs(effect)))]
                ),
                "signed_delta_at_max_abs": float(
                    effect[int(np.argmax(np.abs(effect)))]
                ),
            },
        },
        "interpretation_allowed": accepted,
        "interpretation": (
            "R_Ay < 1 means the complete 3NF reduces the discrepancy; "
            "R_Ay > 1 means it worsens it.  This interpretation is withheld "
            "until every gate passes."
        ),
    }
    return report, theta, two_curve.values, three_curve.values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--two-nf-dir", type=Path, required=True)
    parser.add_argument("--three-nf-dir", type=Path, required=True)
    parser.add_argument("--experiment", type=Path, required=True)
    parser.add_argument("--target-tlab-mev", type=float, default=10.0)
    parser.add_argument("--output", type=Path, required=True,
                        help="JSON output; a sibling CSV curve is also written")
    parser.add_argument("--binding-tolerance-mev", type=float, default=0.02)
    parser.add_argument("--ay-convergence-tolerance", type=float, default=0.002)
    parser.add_argument(
        "--convergence-pair",
        action="append",
        default=[],
        metavar="DIM=TWO_DIR,THREE_DIR",
    )
    parser.add_argument(
        "--require-acceptance",
        action="store_true",
        help="return a failing status when any physics acceptance gate fails",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.binding_tolerance_mev <= 0.0:
        raise SystemExit("binding tolerance must be positive")
    if args.ay_convergence_tolerance <= 0.0:
        raise SystemExit("Ay convergence tolerance must be positive")
    report, theta, two_curve, three_curve = build_report(args)
    output = path_from_repo(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    csv_path = output.with_suffix(".csv")
    write_curve_csv(csv_path, theta, two_curve, three_curve)
    print(json.dumps({
        "accepted_physics_benchmark": report["accepted_physics_benchmark"],
        "failed_gates": [
            name for name, item in report["gates"].items() if not item["passed"]
        ],
        "R_Ay_rmse_ratio": report["metrics"]["R_Ay_rmse_ratio"],
        "json": str(output),
        "csv": str(csv_path),
    }, indent=2))
    if args.require_acceptance and not report["accepted_physics_benchmark"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
