"""Spin-1 polarization observables.

Public entry point for the observable layer. The implementation of
``M(theta) -> observables`` (the spin-density trace formulas
``T_{kq} = Tr[M tau_{kq} M^dagger] / Tr[M M^dagger]`` giving dSigma/dOmega,
iT11, T20, T21, T22) currently lives in :mod:`tictac.amplitudes` because it is
tightly coupled to the assembled scattering amplitude ``M``; this module
re-exports the canonical names so callers can depend on the stable
``tictac.observables`` interface.
"""

from .amplitudes import observables_from_M  # noqa: F401

__all__ = ["observables_from_M"]
