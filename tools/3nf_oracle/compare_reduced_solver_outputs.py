#!/usr/bin/env python3
"""Compare complete-3NF Padé amplitudes with a dense direct solve.

This utility deliberately reports Padé numerical agreement and the solver's
honesty-sidecar status separately.  A small dense cross-check may agree closely
while still being too coarse, or while the configured Padé tail is marked as
maximum-order truncated; neither condition is publication-level convergence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path


U_PREFIX = "U_PW_elements_"
CONVERGENCE_PREFIX = "U_PW_convergence_"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_run_parameters(directory: Path) -> dict[str, str]:
    parameters = {}
    for line in (directory / "run_parameters.txt").read_text(encoding="utf-8").splitlines():
        if ":" not in line or line.startswith("Running program"):
            continue
        key, value = line.split(":", 1)
        parameters[key.strip()] = value.strip()
    return parameters


def parse_u_rows(path: Path) -> tuple[float, list[tuple[float, float, int, list[complex]]]]:
    deuteron_binding = math.nan
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("# Deuteron BE:"):
            deuteron_binding = float(stripped.split(":", 1)[1].split()[0])
            continue
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            tlab = float(parts[0])
            ecm = float(parts[1])
            q_index = int(parts[2])
            amplitudes = [complex(token) for token in parts[3:]]
        except ValueError:
            continue
        rows.append((tlab, ecm, q_index, amplitudes))
    if not rows:
        raise ValueError(f"no U-matrix rows parsed from {path}")
    return deuteron_binding, rows


def convergence_path(u_path: Path) -> Path:
    stem = u_path.stem.replace(U_PREFIX, CONVERGENCE_PREFIX, 1)
    stem = re.sub(r"_PSI_-?\d+$", "", stem)
    return u_path.with_name(stem + u_path.suffix)


def parse_convergence(path: Path) -> list[dict[str, int]]:
    records = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        row, column, q_index, code, best_order = map(int, stripped.split())
        records.append({
            "row": row,
            "column": column,
            "q_index": q_index,
            "code": code,
            "best_pade_order": best_order,
        })
    return records


def compare_file(pade_path: Path, dense_path: Path) -> dict:
    pade_binding, pade_rows = parse_u_rows(pade_path)
    dense_binding, dense_rows = parse_u_rows(dense_path)
    if len(pade_rows) != len(dense_rows):
        raise ValueError(f"row-count mismatch for {pade_path.name}")

    differences = []
    for pade_row, dense_row in zip(pade_rows, dense_rows):
        if pade_row[:3] != dense_row[:3]:
            raise ValueError(f"kinematic-row mismatch for {pade_path.name}")
        if len(pade_row[3]) != len(dense_row[3]):
            raise ValueError(f"amplitude-count mismatch for {pade_path.name}")
        for pade_value, dense_value in zip(pade_row[3], dense_row[3]):
            absolute = abs(pade_value - dense_value)
            relative = absolute / max(abs(dense_value), 1.0e-300)
            differences.append((absolute, relative))

    conv_path = convergence_path(pade_path)
    convergence = parse_convergence(conv_path)
    codes = [item["code"] for item in convergence]
    return {
        "file": pade_path.name,
        "pade_sha256": sha256(pade_path),
        "dense_sha256": sha256(dense_path),
        "convergence_sha256": sha256(conv_path),
        "deuteron_binding_mev": pade_binding,
        "dense_deuteron_binding_mev": dense_binding,
        "num_kinematic_rows": len(pade_rows),
        "num_complex_amplitudes": len(differences),
        "max_absolute_difference_mev": max(value[0] for value in differences),
        "max_relative_difference": max(value[1] for value in differences),
        "rms_absolute_difference_mev": math.sqrt(
            sum(value[0] ** 2 for value in differences) / len(differences)
        ),
        "pade_status": {
            "num_truly_converged": codes.count(1),
            "num_max_order_truncated": codes.count(2),
            "all_truly_converged": bool(codes) and all(code == 1 for code in codes),
            "records": convergence,
        },
    }


def build_report(pade_directory: Path, dense_directory: Path) -> dict:
    pade_parameters = parse_run_parameters(pade_directory)
    dense_parameters = parse_run_parameters(dense_directory)
    parameter_differences = {
        key: {"pade": pade_parameters.get(key), "dense": dense_parameters.get(key)}
        for key in sorted(set(pade_parameters) | set(dense_parameters))
        if pade_parameters.get(key) != dense_parameters.get(key)
    }
    expected_difference = {"Solve Faddeev with LAPACK"}
    unexpected = set(parameter_differences) - expected_difference
    if unexpected:
        raise ValueError(f"solver inputs differ beyond solve method: {sorted(unexpected)}")

    # Cache/output locations are execution details, not physics inputs.  Keep
    # them out of the portable report while preserving hashes of the complete
    # original run-parameter files below.
    for transient_key in ("Cache root", "P123-matrix read/write folder"):
        pade_parameters.pop(transient_key, None)
    energy_path = Path(pade_parameters["Energy input file"])
    try:
        cpp_index = energy_path.parts.index("CPP")
        pade_parameters["Energy input file"] = str(Path(*energy_path.parts[cpp_index:]))
    except ValueError:
        pade_parameters["Energy input file"] = energy_path.name

    pade_files = sorted(pade_directory.glob(f"{U_PREFIX}*.txt"))
    if not pade_files:
        raise ValueError(f"no U-matrix files under {pade_directory}")
    results = []
    for pade_path in pade_files:
        dense_path = dense_directory / pade_path.name
        if not dense_path.is_file():
            raise FileNotFoundError(dense_path)
        results.append(compare_file(pade_path, dense_path))

    return {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat()
                         .replace("+00:00", "Z"),
        "purpose": "algorithm-level complete-N2LO-3NF Pade versus dense cross-check",
        "parameters": pade_parameters,
        "intentional_parameter_differences": parameter_differences,
        "source_hashes": {
            "pade_run_parameters_sha256": sha256(pade_directory / "run_parameters.txt"),
            "dense_run_parameters_sha256": sha256(dense_directory / "run_parameters.txt"),
        },
        "files": results,
        "aggregate": {
            "num_files": len(results),
            "num_complex_amplitudes": sum(item["num_complex_amplitudes"] for item in results),
            "max_absolute_difference_mev": max(
                item["max_absolute_difference_mev"] for item in results
            ),
            "max_relative_difference": max(
                item["max_relative_difference"] for item in results
            ),
            "all_pade_elements_truly_converged": all(
                item["pade_status"]["all_truly_converged"] for item in results
            ),
        },
        "limitations": [
            "Np_WP=4 gives an unphysical deuteron binding energy near -0.0114 MeV.",
            "Np_per_WP_W1=Nq_per_WP_W1=1 is the unconverged midpoint diagnostic.",
            "Nangle_3NF=2 is not an angular-convergence certificate.",
            "Only one J^pi block is compared.",
            "Numerical agreement with dense inversion does not override Pade honesty codes.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pade-dir", type=Path, required=True)
    parser.add_argument("--dense-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-abs-tol", type=float, default=1.0e-7)
    args = parser.parse_args()

    report = build_report(args.pade_dir, args.dense_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")

    maximum = report["aggregate"]["max_absolute_difference_mev"]
    return 0 if maximum <= args.max_abs_tol else 1


if __name__ == "__main__":
    raise SystemExit(main())
