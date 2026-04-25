# ===============================================================
# 抽取自仓库 [current]: examples/extract_phase_shifts.py
# 行号区段：59..135
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
def parse_u_pw_file(path: Path):
    """Return (channels, rows) where:
        channels = list[dict(idx, l, two_j)]
        rows     = list[dict(tlab, ecm, q_idx, U_matrix (ndim x ndim complex))]
    """
    channels = []
    rows = []
    in_table1 = False
    in_table2 = False
    header_cols = None

    with path.open() as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if "Name" in stripped and "row-idx" in stripped and "l'" in stripped:
                    in_table1 = True
                    in_table2 = False
                    continue
                if "Tlab" in stripped and "U00" in stripped:
                    in_table1 = False
                    in_table2 = True
                    # Determine how many complex U columns are present.
                    header_cols = [c for c in stripped.strip("#").split() if c.startswith("U")]
                    continue
                if "####" in stripped:
                    in_table1 = False
                    # keep in_table2 as is; the footer "####" after data closes it
                    continue
                continue

            if in_table1:
                # label row-idx col-idx l' 2j' l 2j
                toks = stripped.split()
                if len(toks) != 7:
                    continue
                label, r, c, lp, jp, l, j = toks
                channels.append(
                    dict(
                        label=label,
                        row=int(r), col=int(c),
                        lp=int(lp), two_jp=int(jp),
                        l=int(l), two_j=int(j),
                    )
                )
                continue

            if in_table2:
                toks = stripped.split()
                if len(toks) < 4:
                    continue
                try:
                    tlab = float(toks[0])
                    ecm = float(toks[1])
                    q_idx = int(toks[2])
                except ValueError:
                    continue
                u_tokens = toks[3:]  # header slot count is a 3x3 template; data is n²
                n = int(round(math.sqrt(len(u_tokens))))
                if n * n != len(u_tokens):
                    raise ValueError(
                        f"U token count {len(u_tokens)} is not a perfect square in {path}"
                    )
                U = np.zeros((n, n), dtype=complex)
                for k, tok in enumerate(u_tokens):
                    U[k // n, k % n] = parse_complex_token(tok)
                rows.append(dict(tlab=tlab, ecm=ecm, q_idx=q_idx, U=U))

    # channels table lists every (row, col) pair; dedupe by row-index.
    by_row = {}
    for c in channels:
        if c["row"] == c["col"]:
            by_row[c["row"]] = dict(row=c["row"], l=c["lp"], two_j=c["two_jp"])
    channels = [by_row[k] for k in sorted(by_row)]
    return channels, rows
