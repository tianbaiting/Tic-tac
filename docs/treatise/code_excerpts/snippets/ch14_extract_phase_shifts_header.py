# ===============================================================
# 抽取自仓库 [current]: examples/extract_phase_shifts.py
# 行号区段：1..50
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
#!/usr/bin/env python3
"""
Extract nd elastic-scattering partial-wave phase shifts from Tic-tac
U_PW_elements_*.txt output. Implements Miller, Ekström, Hebeler
PRC 106, 024001 (2022), Appendix D/E.

Conversion path (Miller eqs. D3, D5):

    U^plane(q0)  =  U^WP  *  f̄²(q0) / (q0² * ΔE_j)
                 =  U^WP  / (q0 * μ1 * ΔE_j)           (energy WPs: f̄² = q/μ1)

    S           =  I  −  2πi · q0 · m_N · U^plane(q0)

    δ_diag      =  ½ arg(S_diag)              (quick diagonal readout)

    S = U exp(2iδ) U^T     (full eigenphase analysis, Blatt–Biedenharn/E2)

For J^Π = 1/2+, the 2-channel doublet contains (l=0, j=1/2) ²S₁/₂ and
(l=2, j=3/2). Miller Fig. 1 at E_lab = 13 MeV (Nijmegen-I) gives the dominant
doublet phase: Re δ ≈ 105°, Im δ ≈ 15°.

Usage:
  python extract_phase_shifts.py --solver-out-dir <dir> [--tlab 13.0]
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np


HBARC_MEV_FM = 197.327
M_PROTON_MEV = 938.272
M_NEUTRON_MEV = 939.565
M_DEUTERON_MEV = 1875.61

# Miller's m_N for the S-matrix prefactor (PRC 106, text after Eq. D3).
MN_MILLER = 2.0 * M_PROTON_MEV * M_NEUTRON_MEV / (M_PROTON_MEV + M_NEUTRON_MEV)

# Reduced mass of the N-d system (spectator-pair system) for f̄²(q) = q/μ1 and
# ΔE_j = (q_{j+1}² - q_j²)/(2 μ1).
MU1_ND = M_NEUTRON_MEV * M_DEUTERON_MEV / (M_NEUTRON_MEV + M_DEUTERON_MEV)


# --------------------------- parsing helpers --------------------------- #

_COMPLEX_RE = re.compile(r"([+\-]\d\.\d+e[+\-]\d+)([+\-]\d\.\d+e[+\-]\d+)j")

