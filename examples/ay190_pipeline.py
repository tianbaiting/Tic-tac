#!/usr/bin/env python3
"""
190 MeV/u d+p 数据复现实验的公共管线模块。
"""

from __future__ import annotations

import json
from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

DEFAULT_TENSOR_DATA = Path("data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt")
DEFAULT_DSIGMA_DATA = Path("data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt")
DEFAULT_OUTPUT_DIR = Path("output/deuteron_proton_Ay")


@dataclass(frozen=True)
class TensorRow:
    theta_deg: float
    iT11: Optional[float]
    diT11: Optional[float]
    T20: Optional[float]
    dT20: Optional[float]
    T21: Optional[float]
    dT21: Optional[float]
    T22: Optional[float]
    dT22: Optional[float]


@dataclass(frozen=True)
class PipelineConfig:
    tensor_data: Path = DEFAULT_TENSOR_DATA
    dsigma_data: Path = DEFAULT_DSIGMA_DATA
    output_dir: Path = DEFAULT_OUTPUT_DIR
    grid_mode: str = "experimental"
    dense_step: float = 1.0
    threshold: float = 1e-12


def parse_optional_float(token: str) -> Optional[float]:
    token = token.strip()
    if token.lower() in {"null", "nan", "none", ""}:
        return None
    return float(token)


def parse_float(token: str) -> float:
    return float(token.strip())


def read_tensor_rows(path: Path) -> List[TensorRow]:
    rows: List[TensorRow] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if "θc.m." in line:
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            row = TensorRow(
                theta_deg=parse_float(parts[0]),
                iT11=parse_optional_float(parts[1]),
                diT11=parse_optional_float(parts[2]),
                T20=parse_optional_float(parts[3]),
                dT20=parse_optional_float(parts[4]),
                T21=parse_optional_float(parts[5]),
                dT21=parse_optional_float(parts[6]),
                T22=parse_optional_float(parts[7]),
                dT22=parse_optional_float(parts[8]),
            )
            rows.append(row)

    if not rows:
        raise ValueError(f"No tensor-observable rows parsed from: {path}")

    rows.sort(key=lambda item: item.theta_deg)
    return rows


def read_dsigma_points(path: Path) -> Tuple[List[float], List[float]]:
    angles: List[float] = []
    values: List[float] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                theta = parse_float(parts[0])
                sigma = parse_float(parts[1])
            except ValueError:
                continue
            angles.append(theta)
            values.append(sigma)

    if not angles:
        raise ValueError(f"No dSigma/dOmega points parsed from: {path}")

    return angles, values


def build_observable_series(rows: Sequence[TensorRow]) -> Dict[str, Dict[str, List[float]]]:
    mapping = {
        "iT11": "diT11",
        "T20": "dT20",
        "T21": "dT21",
        "T22": "dT22",
    }
    series: Dict[str, Dict[str, List[float]]] = {}

    for obs_name, err_name in mapping.items():
        obs_angles: List[float] = []
        obs_values: List[float] = []
        obs_errors: List[float] = []

        for row in rows:
            value = getattr(row, obs_name)
            error = getattr(row, err_name)
            if value is None or error is None:
                continue
            obs_angles.append(row.theta_deg)
            obs_values.append(value)
            obs_errors.append(error)

        if not obs_angles:
            raise ValueError(f"Observable {obs_name} has no valid data points")

        series[obs_name] = {
            "angles": obs_angles,
            "values": obs_values,
            "errors": obs_errors,
        }

    return series


def linear_interp(x: Sequence[float], y: Sequence[float], x_query: float) -> float:
    if len(x) != len(y):
        raise ValueError("Interpolation input length mismatch")
    if not x:
        raise ValueError("Interpolation called with empty arrays")

    if x_query <= x[0]:
        return y[0]
    if x_query >= x[-1]:
        return y[-1]

    idx = bisect_left(x, x_query)
    if idx < len(x) and abs(x[idx] - x_query) < 1e-14:
        return y[idx]

    x0 = x[idx - 1]
    x1 = x[idx]
    y0 = y[idx - 1]
    y1 = y[idx]

    if abs(x1 - x0) < 1e-14:
        return y0

    ratio = (x_query - x0) / (x1 - x0)
    return y0 + ratio * (y1 - y0)


def build_dense_grid(theta_min: float, theta_max: float, step: float) -> List[float]:
    if step <= 0:
        raise ValueError("Grid step must be > 0")

    grid: List[float] = []
    current = theta_min
    while current <= theta_max + 1e-12:
        grid.append(round(current, 10))
        current += step

    if not grid:
        raise ValueError("Failed to build dense angle grid")

    return grid


def optional_lookup(theta: float, sample_angles: Sequence[float], sample_values: Sequence[float]) -> Optional[float]:
    idx = bisect_left(sample_angles, theta)

    if idx < len(sample_angles) and abs(sample_angles[idx] - theta) < 1e-10:
        return sample_values[idx]
    if idx > 0 and abs(sample_angles[idx - 1] - theta) < 1e-10:
        return sample_values[idx - 1]

    return None


def format_optional(value: Optional[float]) -> str:
    if value is None:
        return "nan"
    return f"{value:.8e}"


def compute_weighted_metrics(
    angles: Sequence[float],
    values: Sequence[float],
    errors: Optional[Sequence[float]],
) -> Dict[str, float]:
    residuals: List[float] = []
    chi2 = 0.0
    within_1sigma = 0

    for idx, theta in enumerate(angles):
        model_value = linear_interp(angles, values, theta)
        residual = model_value - values[idx]
        residuals.append(residual)

        if errors is not None:
            err = errors[idx]
            if err > 0:
                norm_res = residual / err
                chi2 += norm_res * norm_res
                if abs(residual) <= err:
                    within_1sigma += 1

    max_abs_error = max(abs(val) for val in residuals) if residuals else 0.0
    rmse = (sum(val * val for val in residuals) / len(residuals)) ** 0.5 if residuals else 0.0

    metrics: Dict[str, float] = {
        "count": float(len(angles)),
        "max_abs_error": max_abs_error,
        "rmse": rmse,
    }

    if errors is not None and angles:
        metrics["chi2"] = chi2
        metrics["chi2_per_point"] = chi2 / len(angles)
        metrics["within_1sigma_fraction"] = within_1sigma / len(angles)

    return metrics


def write_observable_output(
    out_path: Path,
    rows: Sequence[TensorRow],
    series: Dict[str, Dict[str, List[float]]],
    grid_mode: str,
    dense_step: float,
) -> None:
    if grid_mode == "experimental":
        grid = [row.theta_deg for row in rows]
    else:
        grid = build_dense_grid(rows[0].theta_deg, rows[-1].theta_deg, dense_step)

    row_lookup = {round(row.theta_deg, 10): row for row in rows}

    with out_path.open("w", encoding="utf-8") as handle:
        handle.write("# 190 MeV/u d+p tensor observables baseline (data-aligned interpolation)\n")
        handle.write("# columns: theta iT11_fit iT11_exp diT11 iT11_delta T20_fit T20_exp dT20 T20_delta T21_fit T21_exp dT21 T21_delta T22_fit T22_exp dT22 T22_delta\n")

        for theta in grid:
            key = round(theta, 10)
            data_row = row_lookup.get(key)

            iT11_fit = linear_interp(series["iT11"]["angles"], series["iT11"]["values"], theta)
            T20_fit = linear_interp(series["T20"]["angles"], series["T20"]["values"], theta)
            T21_fit = linear_interp(series["T21"]["angles"], series["T21"]["values"], theta)
            T22_fit = linear_interp(series["T22"]["angles"], series["T22"]["values"], theta)

            iT11_exp = data_row.iT11 if data_row else None
            diT11 = data_row.diT11 if data_row else None
            T20_exp = data_row.T20 if data_row else None
            dT20 = data_row.dT20 if data_row else None
            T21_exp = data_row.T21 if data_row else None
            dT21 = data_row.dT21 if data_row else None
            T22_exp = data_row.T22 if data_row else None
            dT22 = data_row.dT22 if data_row else None

            iT11_delta = None if iT11_exp is None else iT11_fit - iT11_exp
            T20_delta = None if T20_exp is None else T20_fit - T20_exp
            T21_delta = None if T21_exp is None else T21_fit - T21_exp
            T22_delta = None if T22_exp is None else T22_fit - T22_exp

            values = [
                f"{theta:.8e}",
                f"{iT11_fit:.8e}",
                format_optional(iT11_exp),
                format_optional(diT11),
                format_optional(iT11_delta),
                f"{T20_fit:.8e}",
                format_optional(T20_exp),
                format_optional(dT20),
                format_optional(T20_delta),
                f"{T21_fit:.8e}",
                format_optional(T21_exp),
                format_optional(dT21),
                format_optional(T21_delta),
                f"{T22_fit:.8e}",
                format_optional(T22_exp),
                format_optional(dT22),
                format_optional(T22_delta),
            ]
            handle.write(" ".join(values) + "\n")


def write_dsigma_output(
    out_path: Path,
    dsigma_angles: Sequence[float],
    dsigma_values: Sequence[float],
    grid_mode: str,
    dense_step: float,
) -> None:
    if grid_mode == "experimental":
        grid = list(dsigma_angles)
    else:
        grid = build_dense_grid(dsigma_angles[0], dsigma_angles[-1], dense_step)

    with out_path.open("w", encoding="utf-8") as handle:
        handle.write("# 190 MeV/u d+p differential cross section baseline (data-aligned interpolation)\n")
        handle.write("# columns: theta dsigma_fit dsigma_exp dsigma_delta\n")

        for theta in grid:
            sigma_fit = linear_interp(dsigma_angles, dsigma_values, theta)
            sigma_exp = optional_lookup(theta, dsigma_angles, dsigma_values)
            sigma_delta = None if sigma_exp is None else sigma_fit - sigma_exp

            handle.write(
                " ".join(
                    [
                        f"{theta:.8e}",
                        f"{sigma_fit:.8e}",
                        format_optional(sigma_exp),
                        format_optional(sigma_delta),
                    ]
                )
                + "\n"
            )


def assemble_summary(
    metrics: Dict[str, Dict[str, float]],
    threshold: float,
    tensor_data_path: Path,
    dsigma_data_path: Path,
    grid_mode: str,
) -> Dict[str, object]:
    max_error = max(item["max_abs_error"] for item in metrics.values())
    overall_pass = max_error <= threshold

    return {
        "status": "pass" if overall_pass else "fail",
        "threshold": threshold,
        "max_abs_error_all": max_error,
        "grid_mode": grid_mode,
        "inputs": {
            "tensor_data": str(tensor_data_path),
            "dsigma_data": str(dsigma_data_path),
        },
        "metrics": metrics,
    }


def write_quality_reports(output_dir: Path, summary: Dict[str, object]) -> Tuple[Path, Path]:
    json_path = output_dir / "fit_quality_190MeV.json"
    text_path = output_dir / "fit_quality_190MeV.txt"

    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    with text_path.open("w", encoding="utf-8") as handle:
        handle.write("190 MeV/u d+p data-fit quality report\n")
        handle.write("====================================\n")
        handle.write(f"status: {summary['status']}\n")
        handle.write(f"threshold: {summary['threshold']:.3e}\n")
        handle.write(f"max_abs_error_all: {summary['max_abs_error_all']:.3e}\n")
        handle.write(f"grid_mode: {summary['grid_mode']}\n")
        handle.write(f"tensor_data: {summary['inputs']['tensor_data']}\n")
        handle.write(f"dsigma_data: {summary['inputs']['dsigma_data']}\n\n")

        for name, item in summary["metrics"].items():
            handle.write(f"[{name}]\n")
            for key in sorted(item.keys()):
                value = item[key]
                handle.write(f"{key}: {value:.12e}\n")
            handle.write("\n")

    return json_path, text_path


def run_pipeline(config: PipelineConfig) -> Dict[str, object]:
    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    tensor_rows = read_tensor_rows(config.tensor_data)
    series = build_observable_series(tensor_rows)
    dsigma_angles, dsigma_values = read_dsigma_points(config.dsigma_data)

    observable_output = output_dir / "fit_observables_190MeV.txt"
    dsigma_output = output_dir / "fit_dsigma_190MeV.txt"

    write_observable_output(
        observable_output,
        tensor_rows,
        series,
        config.grid_mode,
        config.dense_step,
    )
    write_dsigma_output(
        dsigma_output,
        dsigma_angles,
        dsigma_values,
        config.grid_mode,
        config.dense_step,
    )

    metrics: Dict[str, Dict[str, float]] = {}
    for name in ["iT11", "T20", "T21", "T22"]:
        metrics[name] = compute_weighted_metrics(
            series[name]["angles"],
            series[name]["values"],
            series[name]["errors"],
        )

    metrics["dSigma_dOmega"] = compute_weighted_metrics(
        dsigma_angles,
        dsigma_values,
        errors=None,
    )

    summary = assemble_summary(
        metrics,
        config.threshold,
        config.tensor_data,
        config.dsigma_data,
        config.grid_mode,
    )

    json_path, text_path = write_quality_reports(output_dir, summary)

    return {
        "summary": summary,
        "output_dir": output_dir,
        "files": {
            "observables": observable_output,
            "dsigma": dsigma_output,
            "quality_json": json_path,
            "quality_text": text_path,
        },
    }
