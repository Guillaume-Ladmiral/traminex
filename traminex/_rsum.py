"""R's ``sum()`` of doubles, for the modules that reproduce R's rounding.

Private helper shared by :mod:`traminex.disstree` and :mod:`traminex.events._mining`.
"""

from __future__ import annotations

import numpy as np


def r_sum(values) -> float:
    """R's ``sum()`` of doubles, in R's order and precision.

    R adds the values one at a time into a C ``long double`` and rounds the
    total to double once: in R 4.6.1 ``sum(c(1, 1e-16, 1e-16))`` (and
    ``rowSums``/``colSums`` of the same values) is 1.0000000000000002, where
    double arithmetic gives 1. ``np.cumsum`` with ``dtype=np.longdouble`` (the
    platform's C long double) adds in the same order and precision, so the
    result is R's where NumPy's ``longdouble`` is the 80-bit type R uses on
    x86-64; where it is plain double (e.g. Windows builds of NumPy) this is the
    sequential double sum. A matrix is read column by column, as R stores it.
    """
    flat = np.asarray(values, dtype=float).ravel(order="F")
    if flat.size == 0:
        return 0.0
    return float(np.cumsum(flat, dtype=np.longdouble)[-1])
