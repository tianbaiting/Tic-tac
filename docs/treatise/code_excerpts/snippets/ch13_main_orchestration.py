# ===============================================================
# 抽取自仓库 [current]: examples/compare_Ay_experiment.py
# 行号区段：622..720
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
    parser.add_argument("--energy-delta-pass", type=float, default=40.0)
    return parser


def _serialize_complex(z: complex) -> Dict[str, float]:
    return {"re": z.real, "im": z.imag}


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    dsigma_output_unit = normalize_dsigma_unit(args.dsigma_unit)

    # Resolve all IO roots first so the rest of the script only handles absolute paths.
    work_dir = (root / args.work_dir).resolve()
    solver_out_dir = (root / args.solver_out_dir).resolve()
    tensor_data = (root / args.tensor_data).resolve()
    dsigma_data = (root / args.dsigma_data).resolve()

    work_dir.mkdir(parents=True, exist_ok=True)
    run_solver_if_requested(work_dir, args.regenerate, args.target_tlab_mev)

    # Solver branch: U files -> parsed channel points.
    u_files = find_latest_solver_u_files(solver_out_dir)
    if not u_files:
        raise RuntimeError(f"No U_PW_elements files found in {solver_out_dir}")

    points: List[SolverChannelPoint] = []
    for path in u_files:
        points.extend(parse_u_file(path))

    if not points:
        raise RuntimeError("No valid U-matrix rows parsed from solver output")

    # Experiment branch: raw measurement tables -> typed arrays.
    exp_it11 = read_experimental_iT11(tensor_data)
    exp_dsigma = read_experimental_dsigma(dsigma_data, target_unit=dsigma_output_unit)

    # Join branch: combine solver channels by energy and evaluate every available energy point.
    combined = combine_channels_by_energy(points)
    if not combined:
        raise RuntimeError("Failed to combine solver channels into energy points")

    enriched: List[Dict[str, object]] = []
    for item in combined:
        fit = _predict_observables_from_u(
            item["u_eff"],
            exp_it11,
            exp_dsigma,
            dsigma_output_unit=dsigma_output_unit,
        )

        it11_metrics = fit["it11_metrics"]
        dsigma_metrics = fit["dsigma_metrics"]

        enriched_item: Dict[str, object] = {
            "tlab_mev": item["tlab"],
            "ecm_weighted_mev": item["ecm_weighted"],
            "abs_delta_tlab_mev": abs(float(item["tlab"]) - args.target_tlab_mev),
            "it11_mae": it11_metrics["mae"],
            "it11_rmse": it11_metrics["rmse"],
            "it11_max_abs_error": it11_metrics["max_abs_error"],
            "dsigma_mae": dsigma_metrics["mae"],
            "dsigma_rmse": dsigma_metrics["rmse"],
            "dsigma_max_abs_error": dsigma_metrics["max_abs_error"],
            "dsigma_rel_rmse": fit["dsigma_rel_rmse"],
            "phase_sign": fit["phase_sign"],
            "phase_indicator": fit["phase_indicator"],
            "u_norm": fit["u_norm"],
            "inv_it11": fit["inv_it11"],
            "inv_t20": fit["inv_t20"],
            "inv_t21": fit["inv_t21"],
            "inv_t22": fit["inv_t22"],
            "u_eff": {
                "u00": _serialize_complex(item["u_eff"]["u00"]),
                "u01": _serialize_complex(item["u_eff"]["u01"]),
                "u10": _serialize_complex(item["u_eff"]["u10"]),
                "u11": _serialize_complex(item["u_eff"]["u11"]),
            },
            "it11_curve": {
                "angles_deg": exp_it11["angles"],
                "exp": exp_it11["values"],
                "pred": fit["it11_pred"],
            },
            "dsigma_curve": {
                "angles_deg": exp_dsigma["angles"],
                "exp": exp_dsigma["values"],
                "pred": fit["dsigma_pred"],
            },
        }
        enriched.append(enriched_item)

    # Select closest available solver energy to the requested target and emit artifacts.
    best = min(enriched, key=lambda item: float(item["abs_delta_tlab_mev"]))

    pass_flag = (
        float(best["it11_rmse"]) <= args.ay_rmse_pass
        and float(best["dsigma_rel_rmse"]) <= args.dsigma_rel_rmse_pass
        and float(best["abs_delta_tlab_mev"]) <= args.energy_delta_pass
