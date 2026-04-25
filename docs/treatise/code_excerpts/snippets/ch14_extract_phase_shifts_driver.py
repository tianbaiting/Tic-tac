# ===============================================================
# 抽取自仓库 [current]: examples/extract_phase_shifts.py
# 行号区段：222..304
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================

# --------------------------- driver --------------------------- #

def extract_for_file(u_path: Path, q_kin_path: Path, tlab_target: float | None,
                     verbose: bool = True):
    channels, rows = parse_u_pw_file(u_path)
    boundaries, midpoints = parse_q_kinematics(q_kin_path)

    if not rows:
        print(f"[warn] no U rows in {u_path.name}")
        return

    # header metadata
    header_text = u_path.read_text().splitlines()[:12]
    meta = {}
    for ln in header_text:
        if "JP:" in ln:
            meta["JP"] = ln.split(":", 1)[1].strip()
        if "Potential:" in ln:
            meta["Potential"] = ln.split(":", 1)[1].strip()

    print(f"\n=== {u_path.name} ===")
    for k, v in meta.items():
        print(f"  {k}: {v}")
    print(f"  channels (row-idx : l,2j): "
          + ", ".join(f"{c['row']}:{c['l']},{c['two_j']}" for c in channels))

    # For each Tlab row, extract phase shifts
    results = []
    for r in rows:
        q_idx = r["q_idx"]
        q0 = midpoints[q_idx]["q"]                       # MeV
        # Bin width in energy (ΔE_j = E_{j+1}^bnd − E_j^bnd for energy WP)
        # boundaries has Nq+1 entries indexed 0..Nq; bin j spans [j, j+1]
        dE = boundaries[q_idx + 1]["ecm"] - boundaries[q_idx]["ecm"]

        U_wp = r["U"]
        # Be defensive: keep only the N=len(channels) sub-block
        n_ch = len(channels)
        U_wp_chn = U_wp[:n_ch, :n_ch]

        U_plane = wp_to_plane_wave(U_wp_chn, q0, dE)
        S = plane_wave_to_S(U_plane, q0)
        eigs, deltas = eigenphases(S)

        # Also directly extract δ from diagonal S (no mixing)
        delta_diag = 0.5 * np.angle(np.diag(S))
        inelast_diag = np.abs(np.diag(S))

        results.append(dict(
            tlab=r["tlab"], ecm=r["ecm"], q0=q0, dE=dE,
            S=S, eigvals=eigs, delta_eig=deltas,
            delta_diag=delta_diag, inelast_diag=inelast_diag,
        ))

    # Pretty-print selected Tlab
    for res in results:
        if tlab_target is not None and abs(res["tlab"] - tlab_target) > 5.0:
            continue
        print(f"\n  Tlab = {res['tlab']:.3f} MeV  |  Ecm = {res['ecm']:.3f} MeV  "
              f"|  q0 = {res['q0']:.3f} MeV  |  ΔE_j = {res['dE']:.3f} MeV")
        # Debug: print S and asymmetry / unitarity diagnostics
        S = res["S"]
        asym = np.linalg.norm(S - S.T)
        unit_def = np.linalg.norm(S @ S.conj().T - np.eye(S.shape[0]))
        print(f"    ||S - S^T|| = {asym:.4e}   ||SS†-1|| = {unit_def:.4e}")
        # Diagonal readout (no recoupling / no eigenphase mix)
        for k in range(len(res["delta_diag"])):
            d = res["delta_diag"][k]
            lab = channels[k]
            re, im = canonical_delta_deg(d)
            print(f"    δ_diag [row={k}, (l,2j)=({lab['l']},{lab['two_j']})]  "
                  f"= {re:+8.3f}° + i·{im:+7.3f}°   "
                  f"|S_kk| = {res['inelast_diag'][k]:.4f}")
        # Eigenphase (2x2 Blatt–Biedenharn — quick)
        for k, (ev, de) in enumerate(zip(res["eigvals"], res["delta_eig"])):
            re, im = canonical_delta_deg(de)
            print(f"    δ_eig  [k={k}]          "
                  f"= {re:+8.3f}° + i·{im:+7.3f}°   "
                  f"|λ_k|   = {abs(ev):.4f}")

    return results

