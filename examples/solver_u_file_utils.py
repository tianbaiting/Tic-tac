#!/usr/bin/env python3
"""Compatibility shim.

The canonical home for U-file parsing helpers is now the importable package
``tictac`` (``python/tictac/io.py``). This thin re-export keeps every existing
``from solver_u_file_utils import ...`` working unchanged when example scripts
are run as ``python3 examples/<script>.py`` (which places ``examples/`` on
``sys.path``).
"""

import os as _os
import sys as _sys

_PKG = _os.path.normpath(_os.path.join(_os.path.dirname(__file__), "..", "python"))
if _PKG not in _sys.path:
    _sys.path.insert(0, _PKG)

from tictac.io import *  # noqa: F401,F403  (re-export the public helpers)
from tictac.io import (  # noqa: F401  (explicit, for IDE / starred-import safety)
    UFileMeta,
    detect_parity_symbol,
    detect_two_j,
    parse_u_file_meta,
    required_p123_sparse_names,
    select_latest_u_file_family,
)
