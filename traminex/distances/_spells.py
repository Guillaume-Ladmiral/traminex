"""Spell (run-length) encoding helpers shared by the OMspell/OMslen distance family.

A *spell* is a maximal run of an identical state. The DSS (distinct successive states)
representation collapses each spell to a single token; durations record the run lengths.
"""

from __future__ import annotations

import math

import numpy as np


def spell_encode(seq: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the DSS ``(states, durations)`` of a single sequence.

    Example: ``[0, 0, 1, 2, 2]`` -> states ``[0, 1, 2]``, durations ``[2.0, 1.0, 2.0]``.
    """
    seq = np.asarray(seq)
    n = seq.shape[0]
    if n == 0:
        return np.array([], dtype=seq.dtype), np.array([], dtype=float)
    states: list = []
    durs: list[float] = []
    cur = seq[0]
    run = 1
    for t in range(1, n):
        if seq[t] == cur:
            run += 1
        else:
            states.append(cur)
            durs.append(float(run))
            cur = seq[t]
            run = 1
    states.append(cur)
    durs.append(float(run))
    return np.asarray(states, dtype=seq.dtype), np.asarray(durs, dtype=float)


def spell_length_per_position(seq: np.ndarray) -> np.ndarray:
    """Return, for each position, the length of the spell that position belongs to.

    Example: ``[0, 0, 1, 2, 2]`` -> ``[2.0, 2.0, 1.0, 2.0, 2.0]``.
    """
    seq = np.asarray(seq)
    n = seq.shape[0]
    out = np.empty(n, dtype=float)
    _, durs = spell_encode(seq)
    pos = 0
    for d in durs:
        di = int(d)
        out[pos : pos + di] = d
        pos += di
    return out


def spell_length_power(durations: np.ndarray, exponent: float) -> np.ndarray:
    """Spell lengths raised to ``exponent`` as R's ``^`` raises them.

    seqdist uses R's ``^`` for OMslen's ``dur.mat ^ (-1 * h)`` and for the
    ``seqdur(...) ^ tpow`` of OMspell, NMSMST and SVRspell (seqdist.R:619, 626).
    It equals the C library's ``pow``, which :func:`math.pow` calls (7,000 of
    7,000 values: lengths 1-1000, seven exponents), while ``np.power`` differs in
    the last bit for some lengths (TX-AUD-20261019-030). Where ``math.pow``
    raises, R returns a value: an overflow and 0 to a negative power give Inf.
    """
    exponent = float(exponent)
    out = np.empty(np.shape(durations), dtype=float)
    flat = out.reshape(-1)
    for index, value in enumerate(np.asarray(durations, dtype=float).reshape(-1)):
        try:
            flat[index] = math.pow(value, exponent)
        except OverflowError:
            flat[index] = math.inf
        except ValueError:
            flat[index] = math.inf if value == 0.0 and exponent < 0.0 else math.nan
    return out
