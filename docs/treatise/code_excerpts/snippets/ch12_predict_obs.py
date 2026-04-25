# ===============================================================
# 抽取自仓库 [current]: examples/compare_Ay_experiment.py
# 行号区段：310..356
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
def _predict_observables_from_u(
    u_eff: Dict[str, complex],
    exp_it11: Dict[str, List[float]],
    exp_dsigma: Dict[str, object],
    dsigma_output_unit: str = UNIT_MB_PER_SR,
) -> Dict[str, object]:
    u00 = u_eff["u00"]
    u01 = u_eff["u01"]
    u10 = u_eff["u10"]
    u11 = u_eff["u11"]

    u_norm = abs(u00) ** 2 + abs(u01) ** 2 + abs(u10) ** 2 + abs(u11) ** 2
    u_norm = max(u_norm, 1e-18)

    phase_indicator = ((u00 + u11).conjugate() * (u01 - u10)).imag
    phase_sign = 1.0 if phase_indicator >= 0.0 else -1.0

    inv_it11 = phase_indicator / u_norm
    inv_t20 = (abs(u00) ** 2 + abs(u11) ** 2 - abs(u01) ** 2 - abs(u10) ** 2) / u_norm
    inv_t21 = ((u00 - u11).conjugate() * (u01 + u10)).real / u_norm
    inv_t22 = (u00.conjugate() * u11 - u01.conjugate() * u10).real / u_norm

    ay_pred = [_predict_it11(angle, inv_it11) for angle in exp_it11["angles"]]
    sigma_mb = [_predict_dsigma(angle, u_norm, inv_t20, inv_t22) for angle in exp_dsigma["angles"]]
    sigma_pred = [
        convert_dsigma_value(value, UNIT_MB_PER_SR, dsigma_output_unit) for value in sigma_mb
    ]

    ay_metrics = _residual_metrics(ay_pred, exp_it11["values"])
    dsigma_metrics = _residual_metrics(sigma_pred, exp_dsigma["values"])
    dsigma_rel_rmse = _relative_rmse(sigma_pred, exp_dsigma["values"])

    return {
        "u_norm": u_norm,
        "phase_indicator": phase_indicator,
        "phase_sign": phase_sign,
        "inv_it11": inv_it11,
        "inv_t20": inv_t20,
        "inv_t21": inv_t21,
        "inv_t22": inv_t22,
        "it11_pred": ay_pred,
        "dsigma_pred": sigma_pred,
        "it11_metrics": ay_metrics,
        "dsigma_metrics": dsigma_metrics,
        "dsigma_rel_rmse": dsigma_rel_rmse,
    }

