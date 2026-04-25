# ===============================================================
# 抽取自仓库 [current]: examples/extract_phase_shifts.py
# 行号区段：177..220
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
# ----------------------- physics conversion ----------------------- #

def wp_to_plane_wave(U_wp: np.ndarray, q0_mev: float, dE_mev: float) -> np.ndarray:
    """Miller Eq. (D5) for energy WPs:
        U^plane(q0) = U^WP * f̄²(q0) / (q0² ΔE_j)
                    = U^WP / (q0 · μ1 · ΔE_j)
    Units: U^plane has [MeV]⁻².
    """
    return U_wp / (q0_mev * MU1_ND * dE_mev)


def plane_wave_to_S(U_plane: np.ndarray, q0_mev: float) -> np.ndarray:
    """Miller Eq. (D3): S = I − 2πi q0 m_N U^plane (dimensionless)."""
    n = U_plane.shape[0]
    return np.eye(n, dtype=complex) - 2.0j * math.pi * q0_mev * MN_MILLER * U_plane


def eigenphases(S: np.ndarray):
    """Blatt–Biedenharn style eigenphases: S = U exp(2iδ) U^T.
    For complex S above threshold we use symmetric decomposition: Takagi's
    factorization of a complex symmetric matrix. For a quick readout we just
    diagonalize S and take δ_k = ½ arg(λ_k). Mixing-angle reconstruction is
    out of scope for the first sanity check.
    """
    eigvals = np.linalg.eigvals(S)
    deltas = 0.5 * np.angle(eigvals)  # radians, principal branch
    # Also return |λ| as "inelasticity" diagnostic (|λ|=1 means pure elastic).
    return eigvals, deltas


def rad_to_deg(x: complex | float) -> float:
    return float(x) * 180.0 / math.pi


def canonical_delta_deg(delta_complex: complex) -> tuple[float, float]:
    """Map (Re δ, Im δ) into the Miller plot convention:
    Re δ ∈ (0°, 180°] by using δ ≡ δ+180° (same S), pick branch with positive Re.
    Returns (Re_deg, Im_deg).
    """
    re = delta_complex.real * 180.0 / math.pi
    im = delta_complex.imag * 180.0 / math.pi
    # shift Re into [0, 180)
    re = re - 180.0 * math.floor(re / 180.0)
    return re, im
