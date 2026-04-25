# ===============================================================
# 抽取自仓库 [current]: examples/compare_Ay_experiment.py
# 行号区段：195..245
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
def parse_u_file(path: Path) -> List[SolverChannelPoint]:
    parity = detect_parity_from_filename(path)
    two_j = detect_two_j_from_filename(path)
    points: List[SolverChannelPoint] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < 7:
            continue

        try:
            tlab = float(parts[0])
            ecm = float(parts[1])
            q_idx = int(parts[2])
            u00 = complex(parts[3])
            u01 = complex(parts[4])
            u10 = complex(parts[5])
            u11 = complex(parts[6])
        except Exception:
            continue

        f_no_flip = u00 + u11
        f_flip = u01 + u10
        dsigma_proxy = abs(f_no_flip) ** 2 + abs(f_flip) ** 2
        ay_proxy = 0.0
        if dsigma_proxy > 1e-15:
            ay_proxy = (f_no_flip.conjugate() * f_flip).imag / dsigma_proxy

        points.append(
            SolverChannelPoint(
                tlab=tlab,
                ecm=ecm,
                q_idx=q_idx,
                two_j=two_j,
                parity=parity,
                u00=u00,
                u01=u01,
                u10=u10,
                u11=u11,
                ay_proxy=ay_proxy,
                dsigma_proxy=dsigma_proxy,
            )
        )

    return points


