#!/usr/bin/env python3
"""Compatibility shim.

The canonical home for the d+p partial-wave amplitude / spin-1 observable
physics is now the importable package ``tictac``
(``python/tictac/amplitudes.py``). This thin re-export keeps every existing
``from pw_amplitudes import ...`` working unchanged when example scripts are run
as ``python3 examples/<script>.py`` (which places ``examples/`` on
``sys.path``).

See docs/architecture.md §6 for the U^{Jpi} -> M(theta) -> observables chain.
"""

import os as _os
import sys as _sys

_PKG = _os.path.normpath(_os.path.join(_os.path.dirname(__file__), "..", "python"))
if _PKG not in _sys.path:
    _sys.path.insert(0, _PKG)

from tictac.amplitudes import *  # noqa: F401,F403  (re-export the public API)
