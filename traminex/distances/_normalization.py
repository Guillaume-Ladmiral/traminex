"""OM-family distance normalization, ported from TraMineR 2.2-14.

TraMineR applies post-normalization inside each method's C++ ``distance()`` via
``DistanceCalculator::normalizeDistance`` (``SOURCES/TraMineR/src/distancecalculator.h``
lines 45-70). The transform is selected by an integer ``norm`` code and depends on the
raw distance plus a per-method *maximum possible distance* and a pair of *lengths*::

    code 0 none      -> rawdist
    code 1 maxlength -> rawdist / max(l1, l2)
    code 2 gmean     -> 1 - (maxdist - rawdist) / (2 * sqrt(l1) * sqrt(l2))
    code 3 maxdist   -> rawdist / maxdist
    code 4 YujianBo  -> 2 * rawdist / (rawdist + maxdist)

``rawdist == 0`` always returns 0 (identical sequences -> 0). The integer codes follow
``norms[-1]`` order in ``seqdist.R`` (``maxlength, gmean, maxdist, YujianBo``).

The per-method ``maxdist`` / ``l1`` / ``l2`` are transcribed from the corresponding C++
``distance()``; see the per-method helpers below. ``norm="auto"`` resolves per
method (seqdist.R lines 317-326): ``maxlength`` for OM/HAM/DHD, ``gmean`` for
LCS/LCP/RLCP, and ``YujianBo`` for the spell/warp methods (OMloc, OMstran,
OMspell, OMslen, TWED) -- the ``_AUTO_TARGETS`` table below is authoritative.

The shared resolver covers every R-facing method that accepts this vocabulary:
OM, LCS/LCP/RLCP, HAM/DHD, OMloc, TWED, OMslen, OMspell, and OMstran.  Each
dispatcher supplies its method-specific bound and length terms.  TraMineR
*refuses* ``norm`` for NMS/NMSMST/SVRspell and all non-auto CHI2/EUCLID forms;
those boundaries are enforced by the distance validation/row adapters.  The
lower-level ``pairwise`` API retains its historical ``length``/``max`` aliases.
"""

from __future__ import annotations

import math

import numpy as np

from ._spells import spell_length_power

# TraMineR norm names -> integer code (distancecalculator.h switch).
_NORM_CODES = {"maxlength": 1, "gmean": 2, "maxdist": 3, "yujianbo": 4}
# Friendly aliases matching the classic-OM spelling used elsewhere in traminex.
_NORM_ALIASES = {"length": "maxlength", "max": "maxdist"}
# TraMineR's method-specific ``norm="auto"`` choices (seqdist.R lines 317-326).
_AUTO_TARGETS = {
    "OM": "maxlength",
    "HAM": "maxlength",
    "DHD": "maxlength",
    "LCS": "gmean",
    "LCP": "gmean",
    "RLCP": "gmean",
}
_AUTO_TARGET = "yujianbo"


def resolve_norm_code(method: str, norm) -> int:
    """Map a ``norm`` argument to its TraMineR integer code for an OM-family method.

    Accepts ``None``/``"none"`` (-> 0), ``"auto"``, the TraMineR names
    (``maxlength``/``gmean``/``maxdist``/``YujianBo``), and the traminex aliases
    ``length`` (-> maxlength) and ``max`` (-> maxdist). Raises ``ValueError`` otherwise.
    """
    if norm is None:
        return 0
    n = str(norm).strip().lower()
    if n in ("none", ""):
        return 0
    if n == "auto":
        n = _AUTO_TARGETS.get(method.upper(), _AUTO_TARGET)
    n = _NORM_ALIASES.get(n, n)
    code = _NORM_CODES.get(n)
    if code is None:
        raise ValueError(
            f"norm={norm!r} is not supported for method {method}; use one of "
            "None, 'maxlength' ('length'), 'gmean', 'maxdist' ('max'), 'YujianBo', 'auto'"
        )
    return code


def normalize_value(rawdist: float, maxdist: float, l1: float, l2: float, code: int) -> float:
    """Apply ``DistanceCalculator::normalizeDistance`` (distancecalculator.h:45-70)."""
    if rawdist == 0:
        return 0.0
    if code == 0:
        return rawdist
    if code == 1:  # maxlength
        if l1 > l2:
            return rawdist / l1
        if l2 > 0:
            return rawdist / l2
        return 0.0
    if code == 2:  # gmean
        if l1 * l2 == 0:
            return 1.0 if l1 != l2 else 0.0
        return 1.0 - ((maxdist - rawdist) / (2.0 * math.sqrt(l1) * math.sqrt(l2)))
    if code == 3:  # maxdist
        if maxdist == 0:
            return 1.0
        return rawdist / maxdist
    if code == 4:  # YujianBo
        if maxdist == 0:
            return 1.0
        return (2.0 * rawdist) / (rawdist + maxdist)
    return rawdist


def _maxscost(max_offdiag: float, indel: float, code: int) -> float:
    """OMdistance::setParameters maxscost (OMdistance.cpp:23-34).

    For YujianBo (code 4) the maximum substitution cost is fixed at ``2*indel``
    (Yujian & Bo 2007); otherwise it is ``min(max offdiagonal subcost, 2*indel)``.
    """
    if code == 4:
        return 2.0 * indel
    return min(max_offdiag, 2.0 * indel)


def _spell_term_maxima(length: int, tpow: float) -> np.ndarray:
    """For p = 1..length spells, the largest sum of ``t**tpow - 1`` over the spell lengths.

    ``t**tpow - 1`` is convex for tpow >= 1 (and 0 for tpow == 0), so one long spell and
    p - 1 of length 1 maximize the sum; for 0 < tpow < 1 it is concave and the most even
    lengths do (majorization). The powers are :func:`spell_length_power`'s, as in the
    distance itself.
    """
    p = np.arange(1, length + 1)
    if tpow >= 1.0 or tpow == 0.0:
        return spell_length_power(length - p + 1.0, tpow) - 1.0
    base, rest = np.divmod(length, p)
    return rest * (spell_length_power(base + 1.0, tpow) - 1.0) + (p - rest) * (
        spell_length_power(base.astype(float), tpow) - 1.0
    )


def omspell_bound(
    mm: int,
    nn: int,
    *,
    indel: float,
    substitution: float,
    timecost: float,
    tpow: float,
    yujianbo: bool,
) -> float:
    """An upper bound of the OMspell distance between sequences of ``mm`` and ``nn`` positions.

    R bounds it by ``|nn - mm| * indel + maxscost * min(mm, nn)`` (OMPerdistance.cpp:108),
    which leaves out the spell-duration costs (a spell of length t costs
    ``indel + timecost * (t**tpow - 1)`` to insert or delete, OMPerdistance.h) and so is
    exceeded, even at ``timecost`` 0 when the substitution cost is below the indel
    (TX-REP-20261021-004). This is the largest cost, over the spell structures of the two
    lengths, of pairing spells in order (substituting, or deleting and inserting) and
    indelling the rest: with p and q spells,
    ``timecost * (g_p(mm) + g_q(nn)) + min(p, q) * pair + |p - q| * indel``, ``g`` the largest
    duration term and ``pair = min(substitution, 2 * indel)`` (``2 * indel`` under YujianBo, as
    R's convention takes two indels). For tpow < 0 the duration terms are at most 0 and a
    pair of spells of one state costs below ``timecost``, so ``pair`` becomes
    ``min(max(substitution, timecost), 2 * indel)`` and the terms are left out.
    """
    if tpow >= 0.0:
        pair = 2.0 * indel if yujianbo else min(substitution, 2.0 * indel)
        gx = timecost * _spell_term_maxima(mm, tpow) if mm else np.zeros(1)
        gy = timecost * _spell_term_maxima(nn, tpow) if nn else np.zeros(1)
    else:
        pair = 2.0 * indel if yujianbo else min(max(substitution, timecost), 2.0 * indel)
        gx, gy = np.zeros(max(mm, 1)), np.zeros(max(nn, 1))
    p = np.arange(1, mm + 1, dtype=float) if mm else np.zeros(1)
    q = np.arange(1, nn + 1, dtype=float) if nn else np.zeros(1)
    # min(p, q) * pair + |p - q| * indel is p * indel + q * (pair - indel) for q <= p and
    # p * (pair - indel) + q * indel for q >= p: prefix and suffix maxima over q.
    below = np.maximum.accumulate(gy + q * (pair - indel))
    above = np.maximum.accumulate((gy + q * indel)[::-1])[::-1]
    upto = np.searchsorted(q, p, side="right")
    first = np.searchsorted(q, p, side="left")
    best = np.full(p.shape, -np.inf)
    has_below = upto > 0
    best[has_below] = p[has_below] * indel + below[upto[has_below] - 1]
    has_above = first < q.size
    best[has_above] = np.maximum(
        best[has_above], p[has_above] * (pair - indel) + above[first[has_above]]
    )
    return float(np.max(gx + best))
