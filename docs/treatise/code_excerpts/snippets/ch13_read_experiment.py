# ===============================================================
# 抽取自仓库 [current]: examples/compare_Ay_experiment.py
# 行号区段：105..172
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
def _legendre_p(n: int, x: float) -> float:
    if n == 0:
        return 1.0
    if n == 1:
        return x
    p0 = 1.0
    p1 = x
    for k in range(2, n + 1):
        pk = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k
        p0, p1 = p1, pk
    return p1


def read_experimental_iT11(path: Path) -> Dict[str, List[float]]:
    angles: List[float] = []
    values: List[float] = []
    errors: List[float] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "θc.m." in line:
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        if parts[1].lower() == "null" or parts[2].lower() == "null":
            continue

        angles.append(float(parts[0]))
        values.append(float(parts[1]))
        errors.append(float(parts[2]))

    if not values:
        raise ValueError(f"No valid iT11 rows parsed from {path}")

    return {"angles": angles, "values": values, "errors": errors}


def read_experimental_dsigma(path: Path, target_unit: str = UNIT_MB_PER_SR) -> Dict[str, object]:
    lines = path.read_text(encoding="utf-8").splitlines()
    source_unit = infer_dsigma_unit_from_lines(lines, fallback_unit=UNIT_MB_PER_SR)
    output_unit = normalize_dsigma_unit(target_unit)

    angles: List[float] = []
    values: List[float] = []

    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            angles.append(float(parts[0]))
            values.append(float(parts[1]))
        except ValueError:
            continue

    if not values:
        raise ValueError(f"No valid dSigma/dOmega rows parsed from {path}")

    return {
        "angles": angles,
        "values": convert_dsigma_series(values, source_unit, output_unit),
        "source_unit": source_unit,
        "output_unit": output_unit,
    }
