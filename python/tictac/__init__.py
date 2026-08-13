"""Tic-tac few-body observable library.

The partial-wave reconstruction of the d+p elastic scattering amplitude and the
spin-1 polarization observables live here. Example scripts under ``examples/``
are thin drivers that import this package; they no longer contain the physics.

Math flow (see docs/architecture.md §6):
    U^{Jpi}  ->  M(theta)            (tictac.amplitudes.assemble_m_matrix)
    M(theta) ->  observables         (tictac.observables.observables_from_M)
        dSigma/dOmega, iT11, T20, T21, T22
        = Tr[M tau_kq M^dagger] / Tr[M M^dagger]
"""

from .amplitudes import (           # noqa: F401  (public API re-export)
    JPiBlock,
    WPGridBin,
    assemble_m_matrix,
    calibrate_dsigma_scale,
    clebsch_gordan,
    list_solver_energies,
    observables_from_M,
    parse_q_kinematics,
    parse_solver_output,
    wigner_6j,
)
from .io import (                   # noqa: F401
    UFileMeta,
    parse_u_file_meta,
    select_latest_u_file_family,
)

__all__ = [
    "JPiBlock",
    "WPGridBin",
    "assemble_m_matrix",
    "calibrate_dsigma_scale",
    "clebsch_gordan",
    "list_solver_energies",
    "observables_from_M",
    "parse_q_kinematics",
    "parse_solver_output",
    "wigner_6j",
    "UFileMeta",
    "parse_u_file_meta",
    "select_latest_u_file_family",
]
