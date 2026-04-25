# ===============================================================
# 抽取自仓库 [current]: examples/compare_Ay_experiment.py
# 行号区段：250..290
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
def combine_channels_by_energy(points: List[SolverChannelPoint]) -> List[Dict[str, object]]:
    grouped: Dict[float, List[SolverChannelPoint]] = {}
    for item in points:
        key = round(item.tlab, 6)
        grouped.setdefault(key, []).append(item)

    combined: List[Dict[str, object]] = []
    for key in sorted(grouped.keys()):
        block = grouped[key]
        weights = [max(item.dsigma_proxy, 1e-18) for item in block]
        total_w = sum(weights)
        if total_w <= 1e-20:
            continue

        u00_eff = sum(w * item.u00 for w, item in zip(weights, block)) / total_w
        u01_eff = sum(w * item.u01 for w, item in zip(weights, block)) / total_w
        u10_eff = sum(w * item.u10 for w, item in zip(weights, block)) / total_w
        u11_eff = sum(w * item.u11 for w, item in zip(weights, block)) / total_w

        ay_combined = sum(item.ay_proxy * w for item, w in zip(block, weights)) / total_w
        ecm_combined = sum(item.ecm * w for item, w in zip(block, weights)) / total_w

        combined.append(
            {
                "tlab": key,
                "ecm_weighted": ecm_combined,
                "ay_proxy_combined": ay_combined,
                "dsigma_proxy_combined": total_w,
                "num_channels": float(len(block)),
                "u_eff": {
                    "u00": u00_eff,
                    "u01": u01_eff,
                    "u10": u10_eff,
                    "u11": u11_eff,
                },
            }
        )

    return combined


