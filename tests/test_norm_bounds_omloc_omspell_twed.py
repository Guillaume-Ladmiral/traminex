"""Normalization bounds of OMloc, OMspell and TWED (TX-REP-20261021-003 to -005).

seqdist.Rd: "m the maximum possible of the distance given the lengths p and q ... For some distances, m is only a
possibly non-reachable upper bound"; maxdist is d/m, gmean 1 - (m - d)/(2 sqrt(l1 l2)) in R's normalizer and
YujianBo 2d/(d + m), each at most 1 exactly when d <= m. R's m is too small for these three methods, so their
normalized values exceed 1 in R and traminex alike: an R defect corrected. traminex keeps R's raw distances and
lengths l1, l2 and takes as m the cost of an alignment that always exists, priced with the largest cost each of its
operations can have (the construction of R's own OM bound). ``scripts/r_goldens/make_goldens_norm_bounds.R`` records
R's raw and normalized values.

The bounds below are written from that derivation, independently of the code under test; the exhaustive checks
show d <= m for every pair of short sequences over a grid of costs and parameters.

gmean keeps R's own bound (owner decision of 2026-10-22, TX-REP-20261022-002, option A): gmean's lengths keep it at
least 0 only for a bound close to the true maximum, so it is R's value, given R's raw distance, and can exceed 1 where
R's does. maxdist and YujianBo (and so auto) take the corrected bounds.
"""

from __future__ import annotations

import itertools
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "norm_bounds.csv", dtype=str
).set_index(["case", "norm"])
AB = ["A", "B"]
ASYM = np.array([[0.0, 0.1], [3.0, 0.0]])
SYM3 = np.array([[0.0, 3.0], [3.0, 0.0]])
ASYM3 = np.array([[0.0, 1.0, 0.2], [2.5, 0.0, 1.5], [0.4, 3.0, 0.0]])
CONSTANT2 = np.array([[0.0, 2.0], [2.0, 0.0]])
DATA = {
    "a_bbb": (["A", "BBB"], AB),
    "four3": (["BBB", "AAA", "BAB", "ABB"], AB),
    "ragged2": (["BBB", "A", "BAB", "AB"], AB),
    "three": (["ACB", "CC", "BABC", "C", "ABCA"], ["A", "B", "C"]),
}
NORMS = ("maxdist", "gmean", "YujianBo")
BOUND_NORMS = ("maxdist", "YujianBo")  # the norms that take the corrected bound


def _seqs(rows: list[str], alphabet: list[str]) -> tx.STSList:
    width = max(len(row) for row in rows)
    cells = [list(row) + [None] * (width - len(row)) for row in rows]
    return tx.seqdef(np.array(cells, dtype=object), alphabet=alphabet, right="DEL")


def _golden(case: str, norm: str, size: int) -> np.ndarray:
    values = [float(v) for v in GOLDEN.loc[(case, norm), "value"].split()]
    return np.array(values).reshape(size, size, order="F")


def _normalize(d: float, m: float, l1: float, l2: float, norm: str) -> float:
    """R's DistanceCalculator::normalizeDistance (distancecalculator.h)."""
    if d == 0:
        return 0.0
    if norm == "maxdist":
        return 1.0 if m == 0 else d / m
    if norm == "gmean":
        if l1 * l2 == 0:
            return 1.0 if l1 != l2 else 0.0
        return 1.0 - (m - d) / (2.0 * math.sqrt(l1) * math.sqrt(l2))
    return 1.0 if m == 0 else 2.0 * d / (d + m)


def _upper_max(sm: np.ndarray) -> float:
    upper = sm[np.triu_indices(sm.shape[0], k=1)]
    return float(upper.max()) if upper.size else 0.0


def _offdiag_max(sm: np.ndarray) -> float:
    return float(sm[~np.eye(sm.shape[0], dtype=bool)].max())


def _omloc_parts(sm: np.ndarray, expcost: float, context: float, norm: str):
    """R's normalization indel and maxscost, the largest local indel I and the bound's substitution term."""
    indel_r = float(sm.max()) * expcost + context  # seqdist.R:723
    maxscost_r = (
        2 * indel_r if norm == "YujianBo" else min(_upper_max(sm), 2 * indel_r)
    )  # OMdistance.cpp:23-34
    # An indel costs timecost + context * (sm[prev, state] + sm[next, state]) / 2 (OMVIdistance.h).
    local = expcost * maxscost_r + context * float(sm.max())
    sub = 2 * local if norm == "YujianBo" else min(_offdiag_max(sm), 2 * local)
    return indel_r, maxscost_r, local, sub


def _omloc_bound(m: int, n: int, local: float, sub: float) -> float:
    """Substitute along the shorter length (or indel both), indel the rest."""
    return abs(n - m) * local + sub * min(m, n)


def _assert_r_gmean(
    actual: np.ndarray, raw: np.ndarray, r_bound: np.ndarray, lengths: np.ndarray
) -> None:
    """gmean is R's normalizer with R's bound and lengths (owner decision, TX-REP-20261022-002)."""
    for i, j in itertools.product(range(raw.shape[0]), repeat=2):
        expected = _normalize(raw[i, j], r_bound[i, j], lengths[i], lengths[j], "gmean")
        assert actual[i, j] == pytest.approx(expected, rel=1e-12, abs=1e-12), (i, j)


OMLOC_CASES = {
    "omloc_constant_e025": ("a_bbb", CONSTANT2, 0.25, None),
    "omloc_asym_e05": ("four3", ASYM, 0.5, None),
    "omloc_asym_e025": ("four3", ASYM, 0.25, None),
    "omloc_asym_e0": ("four3", ASYM, 0.0, None),
    "omloc_sym_ragged_e01": ("ragged2", SYM3, 0.1, None),
    "omloc_asym3_context15": ("three", ASYM3, 0.3, 1.5),
}


def _omloc_options(sm: np.ndarray, expcost: float, context, norm: str, constant: bool) -> dict:
    options = dict(method="OMloc", subm="CONSTANT" if constant else sm, expcost=expcost, norm=norm)
    if context is not None:
        options["context"] = context
    return options


@pytest.mark.parametrize("case", sorted(OMLOC_CASES))
def test_omloc_raw_distances_are_r_s(case: str) -> None:
    data, sm, expcost, context = OMLOC_CASES[case]
    rows, alphabet = DATA[data]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = np.asarray(
            tx.seqdist(
                _seqs(rows, alphabet),
                **_omloc_options(sm, expcost, context, "none", case == "omloc_constant_e025"),
            ),
            dtype=float,
        )
    np.testing.assert_allclose(raw, _golden(case, "none", len(rows)), rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("norm", BOUND_NORMS)
@pytest.mark.parametrize("case", sorted(OMLOC_CASES))
def test_omloc_normalized_values_use_a_true_bound(case: str, norm: str) -> None:
    data, sm, expcost, context = OMLOC_CASES[case]
    rows, alphabet = DATA[data]
    ctx = 1 - 2 * expcost if context is None else context
    indel_r, maxscost_r, local, sub = _omloc_parts(sm, expcost, ctx, norm)
    r_raw = _golden(case, "none", len(rows))
    r_norm = _golden(case, norm, len(rows))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(
            tx.seqdist(
                _seqs(rows, alphabet),
                **_omloc_options(sm, expcost, context, norm, case == "omloc_constant_e025"),
            ),
            dtype=float,
        )
    for i, j in itertools.product(range(len(rows)), repeat=2):
        m, n = len(rows[i]), len(rows[j])
        if norm == "YujianBo":
            # R's raw distance under YujianBo (its maxscost feeds the local indels), from R's value and R's bound.
            r_bound = abs(n - m) * indel_r + maxscost_r * min(m, n)
            d = 0.0 if r_norm[i, j] == 0 else r_norm[i, j] * r_bound / (2 - r_norm[i, j])
        else:
            d = r_raw[i, j]
        expected = _normalize(d, _omloc_bound(m, n, local, sub), m * indel_r, n * indel_r, norm)
        assert actual[i, j] == pytest.approx(expected, rel=1e-12, abs=1e-14)
        assert actual[i, j] <= 1 + 1e-12
        if r_norm[i, j] > 1 + 1e-9:
            assert actual[i, j] < r_norm[i, j]


def _all_sequences(alphabet: list[str], longest: int) -> list[str]:
    return [
        "".join(p) for k in range(1, longest + 1) for p in itertools.product(alphabet, repeat=k)
    ]


@pytest.mark.parametrize(
    "alphabet, longest", [(AB, 4), (["A", "B", "C"], 3)], ids=["2-states", "3-states"]
)
def test_omloc_bound_holds_for_every_pair(alphabet: list[str], longest: int) -> None:
    rows = _all_sequences(alphabet, longest)
    seqs = _seqs(rows, alphabet)
    rng = np.random.default_rng(20261022)
    k = len(alphabet)
    matrices = [np.where(np.eye(k, dtype=bool), 0.0, 2.0)]
    for symmetric in (True, False, False):
        sm = rng.uniform(0.05, 4.0, (k, k))
        sm = (sm + sm.T) / 2 if symmetric else sm
        np.fill_diagonal(sm, 0.0)
        matrices.append(sm)
    lengths = np.array([len(row) for row in rows])
    checked = 0
    for sm, expcost, context in itertools.product(matrices, (0.0, 0.1, 0.25, 0.5), (None, 1.5)):
        ctx = 1 - 2 * expcost if context is None else context
        options = dict(method="OMloc", subm=sm, expcost=expcost)
        if context is not None:
            options["context"] = context
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            raw = np.asarray(tx.seqdist(seqs, norm="none", **options), dtype=float)
            normalized = {
                norm: np.asarray(tx.seqdist(seqs, norm=norm, **options), dtype=float)
                for norm in NORMS
            }
        _, _, local, sub = _omloc_parts(sm, expcost, ctx, "maxdist")
        bound = np.abs(lengths[:, None] - lengths[None, :]) * local + sub * np.minimum.outer(
            lengths, lengths
        )
        assert np.all(raw <= bound * (1 + 1e-12) + 1e-12)
        for norm in BOUND_NORMS:
            assert np.all(normalized[norm] <= 1 + 1e-12), (norm, expcost, context)
        indel_r, maxscost_r, _, _ = _omloc_parts(sm, expcost, ctx, "gmean")
        r_bound = np.abs(
            lengths[:, None] - lengths[None, :]
        ) * indel_r + maxscost_r * np.minimum.outer(lengths, lengths)
        _assert_r_gmean(normalized["gmean"], raw, r_bound, lengths * indel_r)
        checked += raw.size
    assert checked == len(matrices) * 8 * len(rows) ** 2


# --- OMspell (TX-REP-20261021-004) ---------------------------------------------------------------
S1 = np.array(
    [
        [0, 1.871, 1.957, 2.047],
        [1.871, 0, 3.639, 2.424],
        [1.957, 3.639, 0, 1.551],
        [2.047, 2.424, 1.551, 0],
    ]
)
SMALL = np.array([[0.0, 0.1], [0.1, 0.0]])
ABCD = ["A", "B", "C", "D"]
DATA.update(
    {
        "alternating": (["AAAAAAAAAA", "ABABABABAB"], AB),
        "x1": (["BDDDDC", "CBCBC", "BADBCA", "DCB", "DB", "A"], ABCD),
    }
)
CONSTANT4 = np.where(np.eye(4, dtype=bool), 0.0, 2.0)
OMSPELL_CASES = {
    "omspell_small_e0": ("alternating", dict(subm=SMALL, indel=1.0, expcost=0.0), SMALL),
    "omspell_x1_tpow2": ("x1", dict(subm=S1, indel=0.5, expcost=1.0, tpow=2.0), S1),
    "omspell_x1_tpow05": ("x1", dict(subm=S1, indel=1.0, expcost=2.0, tpow=0.5), S1),
    "omspell_x1_tpowm05": ("x1", dict(subm=S1, indel=1.0, expcost=2.0, tpow=-0.5), S1),
    "omspell_x1_default": ("x1", dict(subm="CONSTANT"), CONSTANT4),
}


def _compositions(length: int):
    """Every way to cut ``length`` positions into spells (ordered, each at least 1)."""
    if length == 0:
        yield ()
        return
    for cuts in itertools.product((False, True), repeat=length - 1):
        parts, run = [], 1
        for cut in cuts:
            if cut:
                parts.append(run)
                run = 1
            else:
                run += 1
        yield tuple(parts + [run])


def _omspell_bound(
    mm: int, nn: int, indel: float, sub: float, expcost: float, tpow: float, yujianbo: bool
) -> float:
    """The largest cost, over every spell structure of the two lengths, of pairing spells in order (substituting,
    or deleting and inserting) and indelling the rest; a spell of length t costs indel + expcost * (t**tpow - 1)
    (OMPerdistance.h). For tpow < 0 the duration terms are at most 0 and the duration difference of two spells of
    one state is below 1, so a pair costs at most min(max(S, expcost), 2 * indel)."""

    def term(t: int) -> float:
        return math.pow(t, tpow) - 1.0

    best = -math.inf
    for x in _compositions(mm):
        for y in _compositions(nn):
            paired, rest = min(len(x), len(y)), abs(len(x) - len(y))
            if tpow >= 0:
                pair = 2 * indel if yujianbo else min(sub, 2 * indel)
                cost = (
                    expcost * (sum(map(term, x)) + sum(map(term, y))) + paired * pair + rest * indel
                )
            else:
                pair = 2 * indel if yujianbo else min(max(sub, expcost), 2 * indel)
                cost = paired * pair + rest * indel
            best = max(best, cost)
    return best


@pytest.mark.parametrize("case", sorted(OMSPELL_CASES))
def test_omspell_raw_distances_are_r_s(case: str) -> None:
    data, options, _ = OMSPELL_CASES[case]
    rows, alphabet = DATA[data]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = np.asarray(
            tx.seqdist(_seqs(rows, alphabet), method="OMspell", norm="none", **options), dtype=float
        )
    np.testing.assert_allclose(raw, _golden(case, "none", len(rows)), rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("norm", BOUND_NORMS)
@pytest.mark.parametrize("case", sorted(OMSPELL_CASES))
def test_omspell_normalized_values_use_a_true_bound(case: str, norm: str) -> None:
    data, options, sm = OMSPELL_CASES[case]
    rows, alphabet = DATA[data]
    indel = float(
        options.get("indel", 1.0)
    )  # seqdist's "auto" indel for CONSTANT costs: cval / 2 = 1
    expcost, tpow = float(options.get("expcost", 1.0)), float(options.get("tpow", 1.0))
    r_raw, r_norm = _golden(case, "none", len(rows)), _golden(case, norm, len(rows))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(
            tx.seqdist(_seqs(rows, alphabet), method="OMspell", norm=norm, **options), dtype=float
        )
    for i, j in itertools.product(range(len(rows)), repeat=2):
        m, n = len(rows[i]), len(rows[j])
        bound = _omspell_bound(m, n, indel, _offdiag_max(sm), expcost, tpow, norm == "YujianBo")
        expected = _normalize(r_raw[i, j], bound, m * indel, n * indel, norm)
        assert actual[i, j] == pytest.approx(expected, rel=1e-12, abs=1e-14)
        assert actual[i, j] <= 1 + 1e-12
        if r_norm[i, j] > 1 + 1e-9:
            assert actual[i, j] < r_norm[i, j]


def test_omspell_defaults_keep_r_s_values() -> None:
    """With R's defaults (CONSTANT costs, indel 1, expcost 1, tpow 1) the bound is R's own."""
    rows, alphabet = DATA["x1"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for norm in NORMS:
            actual = np.asarray(
                tx.seqdist(_seqs(rows, alphabet), method="OMspell", subm="CONSTANT", norm=norm),
                dtype=float,
            )
            np.testing.assert_allclose(
                actual, _golden("omspell_x1_default", norm, len(rows)), rtol=1e-12, atol=1e-14
            )


def test_omspell_bound_holds_for_every_pair() -> None:
    rows = _all_sequences(AB, 5)
    seqs = _seqs(rows, AB)
    lengths = [len(row) for row in rows]
    rng = np.random.default_rng(20261022)
    matrices = [SMALL, np.array([[0.0, 3.0], [3.0, 0.0]]), np.array([[0.0, 0.4], [2.6, 0.0]])]
    grid = list(
        itertools.product(range(3), (-0.5, 0.0, 0.5, 1.0, 2.0), (0.0, 0.5, 2.0), (0.5, 1.0, 2.0))
    )
    checked = 0
    for index in rng.choice(len(grid), size=30, replace=False):
        which, tpow, expcost, indel = grid[index]
        sm = matrices[which]
        options = dict(method="OMspell", subm=sm, indel=indel, expcost=expcost, tpow=tpow)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            raw = np.asarray(tx.seqdist(seqs, norm="none", **options), dtype=float)
            normalized = {
                norm: np.asarray(tx.seqdist(seqs, norm=norm, **options), dtype=float)
                for norm in NORMS
            }
        bounds = {
            (m, n): _omspell_bound(m, n, indel, _offdiag_max(sm), expcost, tpow, False)
            for m in range(1, 6)
            for n in range(1, 6)
        }
        bound = np.array([[bounds[(m, n)] for n in lengths] for m in lengths])
        assert np.all(raw <= bound * (1 + 1e-12) + 1e-12), (tpow, expcost, indel, which)
        for norm in BOUND_NORMS:
            assert np.all(normalized[norm] <= 1 + 1e-12), (norm, tpow, expcost, indel, which)
        size = np.array(lengths)
        r_bound = np.abs(size[:, None] - size[None, :]) * indel + min(
            _upper_max(sm), 2 * indel
        ) * np.minimum.outer(size, size)
        _assert_r_gmean(normalized["gmean"], raw, r_bound, size * indel)
        checked += raw.size
    assert checked == 30 * len(rows) ** 2


def test_omspell_bound_closed_forms_equal_the_enumeration() -> None:
    """The production bound (closed-form spell-term maxima, prefix and suffix maxima over the spell counts) equals
    the enumeration over every spell structure, for lengths 0-6 and a grid of parameters."""
    from traminex.distances._normalization import omspell_bound

    for tpow, expcost, indel, sub, yujianbo in itertools.product(
        (-0.5, 0.0, 0.3, 1.0, 2.5), (0.0, 0.7, 2.0), (0.5, 1.0), (0.1, 1.5, 4.0), (False, True)
    ):
        for m, n in itertools.product(range(7), repeat=2):
            expected = _omspell_bound(m, n, indel, sub, expcost, tpow, yujianbo) if m or n else 0.0
            actual = omspell_bound(
                m, n, indel=indel, substitution=sub, timecost=expcost, tpow=tpow, yujianbo=yujianbo
            )
            assert actual == pytest.approx(expected, rel=1e-12, abs=1e-12), (
                m,
                n,
                tpow,
                expcost,
                indel,
                sub,
            )


# --- TWED (TX-REP-20261021-005) ------------------------------------------------------------------
S2 = np.array(
    [
        [0, 2.631, 1.99, 1.159],
        [2.631, 0, 3.697, 2.168],
        [1.99, 3.697, 0, 1.425],
        [1.159, 2.168, 1.425, 0],
    ]
)
DATA.update(
    {
        "x2": (["BDBAB", "DAA", "A", "DABDC", "CCDADDC", "ABBBDA"], ABCD),
        "empty_aba": (["", "ABA"], AB),
    }
)
TWED_CASES = {
    "twed_x1_h2_nu2": ("x1", dict(subm=S1, indel=0.5, h=2.0, nu=2.0), S1),
    "twed_x2_h0_nu05": ("x2", dict(subm=S2, indel=0.5, h=0.0, nu=0.5), S2),
    "twed_empty": ("empty_aba", dict(subm="CONSTANT", nu=0.5, h=0.5), CONSTANT2),
    "twed_x1_auto_indel": ("x1", dict(subm=S1, h=0.5, nu=0.5), S1),
}


def _twed_bound(
    m: int, n: int, nu: float, h: float, sub: float, indel: float, yujianbo: bool
) -> float:
    """Match along the shorter length (each match costs at most twice the largest cost, with no time penalty on the
    diagonal) and warp the rest (each warp costs at most the largest cost + nu + h), priced as R's own bound
    (TWEDdistance.cpp:153) with the uncapped largest cost; an empty sequence is at R's fixed |n - m| * indel.
    """
    if min(m, n) == 0:
        return abs(n - m) * indel
    c = max(sub, 2 * indel) if yujianbo else sub
    return abs(n - m) * (nu + h + c) + 2 * (c + nu) * min(m, n)


def _twed_indel(options: dict, sm: np.ndarray) -> float:
    if "indel" in options:
        return float(options["indel"])
    return (
        2 * float(sm.max()) + float(options["nu"]) + float(options["h"])
    )  # seqdist.R:359-365, 436-439


def _rows_with_empty(rows: list[str], alphabet: list[str]) -> tx.STSList:
    width = max(1, max(len(row) for row in rows))
    cells = [list(row) + [None] * (width - len(row)) for row in rows]
    return tx.seqdef(np.array(cells, dtype=object), alphabet=alphabet, right="DEL")


@pytest.mark.parametrize("case", sorted(TWED_CASES))
def test_twed_raw_distances_are_r_s(case: str) -> None:
    data, options, _ = TWED_CASES[case]
    rows, alphabet = DATA[data]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = np.asarray(
            tx.seqdist(_rows_with_empty(rows, alphabet), method="TWED", norm="none", **options),
            dtype=float,
        )
    np.testing.assert_allclose(raw, _golden(case, "none", len(rows)), rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("norm", BOUND_NORMS)
@pytest.mark.parametrize("case", sorted(TWED_CASES))
def test_twed_normalized_values_use_a_true_bound(case: str, norm: str) -> None:
    data, options, sm = TWED_CASES[case]
    rows, alphabet = DATA[data]
    indel = _twed_indel(options, sm)
    r_raw, r_norm = _golden(case, "none", len(rows)), _golden(case, norm, len(rows))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(
            tx.seqdist(_rows_with_empty(rows, alphabet), method="TWED", norm=norm, **options),
            dtype=float,
        )
    for i, j in itertools.product(range(len(rows)), repeat=2):
        m, n = len(rows[i]), len(rows[j])
        bound = _twed_bound(
            m, n, options["nu"], options["h"], _offdiag_max(sm), indel, norm == "YujianBo"
        )
        expected = _normalize(r_raw[i, j], bound, m * indel, n * indel, norm)
        assert actual[i, j] == pytest.approx(expected, rel=1e-12, abs=1e-14)
        assert actual[i, j] <= 1 + 1e-12
        if r_norm[i, j] > 1 + 1e-9:
            assert actual[i, j] < r_norm[i, j]


def test_twed_with_an_empty_sequence_normalizes_to_one() -> None:
    """R gives d = n * indel against an empty sequence (TWEDdistance.cpp:144-149), the only value it can take."""
    rows, alphabet = DATA["empty_aba"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for norm, expected in (("maxdist", 1.0), ("YujianBo", 1.0), ("gmean", 1.0)):
            value = np.asarray(
                tx.seqdist(
                    _rows_with_empty(rows, alphabet),
                    method="TWED",
                    subm="CONSTANT",
                    nu=0.5,
                    h=0.5,
                    norm=norm,
                )
            )[0, 1]
            assert value == expected, norm


def test_twed_bound_holds_for_every_pair() -> None:
    rows = [""] + _all_sequences(AB, 4)
    seqs = _rows_with_empty(rows, AB)
    lengths = [len(row) for row in rows]
    matrices = [CONSTANT2, SMALL, np.array([[0.0, 0.4], [2.6, 0.0]])]
    checked = 0
    for sm, nu, h, indel in itertools.product(
        matrices, (0.1, 0.5, 2.0), (0.0, 0.5, 2.0), (0.5, None)
    ):
        options = dict(method="TWED", subm=sm, nu=nu, h=h)
        if indel is not None:
            options["indel"] = indel
        used_indel = _twed_indel(options, sm)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            raw = np.asarray(tx.seqdist(seqs, norm="none", **options), dtype=float)
            normalized = {
                norm: np.asarray(tx.seqdist(seqs, norm=norm, **options), dtype=float)
                for norm in NORMS
            }
        bound = np.array(
            [
                [_twed_bound(m, n, nu, h, _offdiag_max(sm), used_indel, False) for n in lengths]
                for m in lengths
            ]
        )
        assert np.all(raw <= bound * (1 + 1e-12) + 1e-12), (nu, h, indel)
        for norm in BOUND_NORMS:
            assert np.all(normalized[norm] <= 1 + 1e-12), (norm, nu, h, indel)
        size = np.array(lengths)
        c_r = min(_upper_max(sm), 2 * used_indel)
        r_bound = np.abs(size[:, None] - size[None, :]) * (nu + h + c_r) + 2 * (
            c_r + nu
        ) * np.minimum.outer(size, size)
        _assert_r_gmean(normalized["gmean"], raw, r_bound, size * used_indel)
        checked += raw.size
    assert checked == 54 * len(rows) ** 2


# TX-REP-20261022-003 (batch-10 review I-2): with tpow < 0 a spell's indel cost
# indel + expcost * (t**tpow - 1) (OMPerdistance.h) is negative once expcost * (1 - t**tpow) exceeds the indel, and R
# then returns a negative distance, which no normalization can bound. traminex refuses such data, as it refuses a
# negative indel or substitution cost; a zero spell cost stays valid (omspell_x1_tpowm05 has one, for DDDD).
SM4 = np.array([[0.0, 4.0], [4.0, 0.0]])


# R's own cost diagnostic (a cost above twice the indel) precedes the refusal.
@pytest.mark.filterwarnings("ignore:.*substitution cost greater than twice the indel")
@pytest.mark.parametrize("norm", ("none",) + NORMS)
def test_omspell_refuses_a_negative_spell_indel_cost(norm: str) -> None:
    rows = ["AABBAABB", "BBAABBAA"]
    r_raw = _golden("omspell_negative_spell_indel", "none", 2)
    assert r_raw[0, 1] == pytest.approx(-15.6, rel=1e-12)  # R's distance is negative
    # Every spell lasts 2: 0.3 + 3 * (2**-2 - 1) = -1.95.
    with pytest.raises(ValueError, match="spell indel cost"):
        tx.seqdist(
            _seqs(rows, AB),
            method="OMspell",
            subm=SM4,
            indel=0.3,
            expcost=3.0,
            tpow=-2.0,
            norm=norm,
        )


@pytest.mark.filterwarnings("ignore:.*substitution cost greater than twice the indel")
def test_omspell_refusal_names_the_longest_spell_and_its_cost() -> None:
    with pytest.raises(ValueError, match=r"length 2 .*-1\.95"):
        tx.seqdist(
            _seqs(["AABBAABB", "BBAABBAA"], AB),
            method="OMspell",
            subm=SM4,
            indel=0.3,
            expcost=3.0,
            tpow=-2.0,
        )


# TX-REP-20261022-002 (batch-10 review I-1; owner decision of 2026-10-22, option A): gmean keeps R's bound for
# OMloc, OMspell and TWED, so it is R's value given R's raw distance, also where the corrected bound would have made it
# negative (OMloc A vs AAB, TWED A vs B with asymmetric costs) and for asymmetric costs generally.
REPO = Path(__file__).resolve().parents[1]
ALL_CASES = {**{c: "OMloc" for c in OMLOC_CASES}, **{c: "OMspell" for c in OMSPELL_CASES}}
ALL_CASES.update({c: "TWED" for c in TWED_CASES})


@pytest.mark.parametrize("case", sorted(ALL_CASES))
def test_gmean_is_r_s_value(case: str) -> None:
    method = ALL_CASES[case]
    if method == "OMloc":
        data, sm, expcost, context = OMLOC_CASES[case]
        options = _omloc_options(sm, expcost, context, "gmean", case == "omloc_constant_e025")
        rows, alphabet = DATA[data]
        seqs = _seqs(rows, alphabet)
    else:
        data, options, _ = (OMSPELL_CASES if method == "OMspell" else TWED_CASES)[case]
        options = dict(options, method=method, norm="gmean")
        rows, alphabet = DATA[data]
        seqs = _rows_with_empty(rows, alphabet)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(tx.seqdist(seqs, **options), dtype=float)
    np.testing.assert_allclose(actual, _golden(case, "gmean", len(rows)), rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize(
    ("case", "rows", "options", "r_value"),
    [
        (
            "omloc_gmean_a_aab",
            ["A", "AAB"],
            dict(method="OMloc", subm=SM4, expcost=0.0, context=0.5),
            1.0,
        ),
        (
            "twed_gmean_a_b",
            ["A", "B"],
            dict(method="TWED", subm=np.array([[0.0, 0.4], [2.5, 0.0]]), indel=2.0, nu=1.0, h=0.0),
            0.4,
        ),
    ],
    ids=["omloc-A-AAB", "twed-A-B-asymmetric"],
)
def test_gmean_keeps_r_s_value_where_the_corrected_bound_made_it_negative(
    case, rows, options, r_value
) -> None:
    assert _golden(case, "gmean", 2)[0, 1] == pytest.approx(r_value, rel=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(tx.seqdist(_seqs(rows, AB), norm="gmean", **options), dtype=float)
    assert actual[0, 1] == pytest.approx(r_value, rel=1e-12)


def test_the_documents_say_gmean_keeps_r_s_bound() -> None:
    def flat(path: str) -> str:
        return " ".join((REPO / path).read_text(encoding="utf-8").split())

    for path in (
        "docs/DISTANCES_NORMALIZATION.md",
        "CHANGELOG.md",
        "docs/API_PARITY_DECISIONS.csv",
    ):
        text = flat(path)
        assert "traminex −2.46" not in text and "traminex -2.46" not in text, path
        assert "gmean keeps R's bound" in text or "`gmean` keeps R's bound" in text, path


def test_the_maximum_cost_helpers_and_the_consumed_surface_note_the_gmean_bound() -> None:
    # Follow-up review m-3 and m-5.
    from traminex.distances import python_kernels

    assert "gmean of OMloc, OMspell and TWED" in " ".join(
        python_kernels._subm_max_offdiag.__doc__.split()
    )
    assert "except their gmean" in " ".join(python_kernels._subm_max_bound.__doc__.split())
    consumed = " ".join((REPO / "docs" / "CONSUMED_SURFACE.md").read_text(encoding="utf-8").split())
    assert "TX-REP-20261022-001 to -003" in consumed and "`gmean` keeps R's bound" in consumed
