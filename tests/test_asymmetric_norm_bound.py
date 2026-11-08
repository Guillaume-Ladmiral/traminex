"""Normalized OM-family distances stay within the maximum possible distance (TX-AUD-20261019-008).

R golden: ``scripts/r_goldens/make_goldens_asym_bound.R`` (``asym_bound.csv``) records R's values. R bounds a
distance with ``maxscost`` read from the strict upper triangle of ``sm`` (OMdistance.cpp:26-33), and seqdist.Rd
defines that bound as the maximum possible distance; with a larger cost below the diagonal R's maxdist reaches 3.
traminex takes the largest off-diagonal cost of the whole matrix (DIVERGENCE note of seqdist). Expected values:
OM's bound worked by hand and by its formula; for OM and OMslen the defining property (no normalized value above
1) over 200 random matrices with varied indels; for OMspell and TWED that property at the parameters of
``EXTRA`` only: with others their normalized values exceed 1 in R and traminex alike, symmetric costs included,
and so do OMloc's with expcost below 0.5 (batch-8 review I-1, I-2; TX-REP-20261021-003 to -005). Raw distances,
YujianBo, maxlength, OMloc and OMstran stay R's, and so does the gmean of OMspell and TWED: it keeps R's bound, read
from the upper triangle (owner decision of 2026-10-22, TX-REP-20261022-002), so it is checked against R's goldens and
R's formula here, not against 1.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex as tx

GOLDEN = pd.read_csv(R_GOLDENS_DIR / "asym_bound.csv", dtype={"value": str})
HAND = np.array([[0.0, 1.0], [3.0, 0.0]])
AUDIT_SM = np.array([[0, 1, 1, 1], [3, 0, 1, 1], [3, 1, 0, 1], [3, 1, 1, 0]], dtype=float)
AUDIT_ROWS = ["ABCD", "BBCD", "DCBA", "AAAA", "CCCC", "BADC"]
# Inputs of the largest cases of a seeded random search (the values compared are R's).
SEARCH = {
    "search_twed": (
        "TWED",
        np.array(
            [
                [0.0, 0.6733111123483121, 0.21390833763794703, 0.7168896906456412],
                [3.8856502611627395, 0.0, 0.7410285674363155, 0.43170017778345815],
                [2.6979008942483964, 3.1296502410905918, 0.0, 0.3394614040410233],
                [3.85890450310538, 3.58931482401455, 2.4626499838039315, 0.0],
            ]
        ),
        ["BCD", "CCDC", "DADBD", "CBAB", "DDAAB", "DDAB"],
    ),
    "search_om": (
        "OM",
        np.array(
            [
                [0.0, 0.6733111123483121, 0.21390833763794703, 0.7168896906456412],
                [3.8856502611627395, 0.0, 0.7410285674363155, 0.43170017778345815],
                [2.6979008942483964, 3.1296502410905918, 0.0, 0.3394614040410233],
                [3.85890450310538, 3.58931482401455, 2.4626499838039315, 0.0],
            ]
        ),
        ["BCD", "CCDC", "DADBD", "CBAB", "DDAAB", "DDAB"],
    ),
}
EXTRA = {
    "OM": {"indel": 2.0},
    "OMslen": {"indel": 2.0},
    "OMspell": {"indel": 2.0},
    "TWED": {"indel": 2.0, "h": 0.5, "nu": 0.5},
    "OMloc": {},
    "OMstran": {"indel": 2.0, "otto": 0.5},
}
CORRECTED = ("OM", "OMslen", "OMspell", "TWED")
BOUND_NORMS = ("maxdist", "gmean")
R_GMEAN = ("OMspell", "TWED")  # gmean keeps R's upper-triangle bound (TX-REP-20261022-002)


def _corrected_norms(method: str) -> tuple[str, ...]:
    return ("maxdist",) if method in R_GMEAN else BOUND_NORMS


def _quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args, **kwargs)


def _r(case: str, method: str, norm: str) -> np.ndarray:
    rows = GOLDEN[
        (GOLDEN["case"] == case) & (GOLDEN["method"] == method) & (GOLDEN["norm"] == norm)
    ]
    assert not rows.empty, (case, method, norm)
    values = rows["value"].astype(float).to_numpy()
    if (rows["j"] == 0).all():
        return values[np.argsort(rows["i"].to_numpy())]
    size = int(rows["i"].max())
    matrix = np.zeros((size, size))
    matrix[rows["i"].to_numpy() - 1, rows["j"].to_numpy() - 1] = values
    return matrix


def _sequences(rows: list[str]):
    return _quiet(tx.seqdef, [list(row) for row in rows], alphabet=list("ABCD"), right="DEL")


def _ours(x, method: str, sm: np.ndarray, norm: str, **kwargs) -> np.ndarray:
    return np.asarray(
        _quiet(tx.seqdist, x, method=method, subm=sm, norm=norm, **EXTRA[method], **kwargs), float
    )


@pytest.mark.parametrize(
    "case, rows, maxdist, gmean",
    [("hand_b_first", ["BB", "AA"], 1.0, 1.0), ("hand_a_first", ["AA", "BB"], 1 / 3, 0.5)],
)
def test_om_bound_worked_by_hand(case: str, rows: list[str], maxdist: float, gmean: float) -> None:
    # d(BB, AA) = 2 x 3 = 6 (a substitution B->A costs 3, less than two indels, 4); d(AA, BB) = 2 x 1 = 2.
    # The bound m = min(3, 2 x 2) x 2 = 6 (the whole off-diagonal maximum 3, not R's upper-triangle 1);
    # gmean = 1 - (m - d) / (2 sqrt(p q)) with p = q = indel x length = 4.
    x = _quiet(tx.seqdef, [list(row) for row in rows], alphabet=["A", "B"])
    for norm, expected in (("maxdist", maxdist), ("gmean", gmean)):
        value = float(_quiet(tx.seqdist, x, method="OM", subm=HAND, indel=2.0, norm=norm)[0, 1])
        assert value == pytest.approx(expected, abs=1e-12)
        assert value != pytest.approx(float(_r(case, "OM", norm)[0, 1]), abs=1e-9)


def test_om_bound_follows_its_formula_on_ragged_data() -> None:
    method, sm, rows = SEARCH["search_om"]
    x = _sequences(rows)
    raw = _ours(x, "OM", sm, "none")
    lengths = np.array([len(row) for row in rows], dtype=float)
    max_cost = min(sm[~np.eye(4, dtype=bool)].max(), 2 * 2.0)
    for i in range(len(rows)):
        for j in range(len(rows)):
            if i == j:
                continue
            bound = abs(lengths[i] - lengths[j]) * 2.0 + max_cost * min(lengths[i], lengths[j])
            assert _ours(x, "OM", sm, "maxdist")[i, j] == pytest.approx(
                raw[i, j] / bound, abs=1e-12
            )


def _cases():
    yield "audit", AUDIT_SM, AUDIT_ROWS, CORRECTED + ("OMloc", "OMstran")
    for case, (method, sm, rows) in SEARCH.items():
        yield case, sm, rows, (method,)


@pytest.mark.parametrize(
    "case, sm, rows, methods", list(_cases()), ids=lambda v: v if isinstance(v, str) else ""
)
def test_raw_distances_and_unaffected_norms_stay_r_s(case, sm, rows, methods) -> None:
    x = _sequences(rows)
    norms = ("none", "YujianBo", "maxlength") if case == "audit" else ("none",)
    for method in methods:
        for norm in norms:
            np.testing.assert_allclose(
                _ours(x, method, sm, norm), _r(case, method, norm), rtol=1e-12, atol=R_GOLDEN_ATOL
            )
        if method in ("OMloc", "OMstran"):
            for norm in BOUND_NORMS:
                np.testing.assert_allclose(
                    _ours(x, method, sm, norm),
                    _r(case, method, norm),
                    rtol=1e-12,
                    atol=R_GOLDEN_ATOL,
                )


@pytest.mark.parametrize(
    "case, sm, rows, methods", list(_cases()), ids=lambda v: v if isinstance(v, str) else ""
)
def test_normalized_values_stay_within_the_bound(case, sm, rows, methods) -> None:
    x = _sequences(rows)
    for method in (m for m in methods if m in CORRECTED):
        for norm in _corrected_norms(method):
            ours = _ours(x, method, sm, norm)
            assert np.nanmax(ours) <= 1 + 1e-12, (method, norm)
            r_values = _r(case, method, norm)
            if np.nanmax(r_values) > 1 + 1e-9:  # R's defect, recorded
                assert not np.allclose(ours, r_values)


@pytest.mark.parametrize(
    "case, sm, rows, methods", list(_cases()), ids=lambda v: v if isinstance(v, str) else ""
)
def test_gmean_of_omspell_and_twed_is_r_s(case, sm, rows, methods) -> None:
    x = _sequences(rows)
    checked = 0
    for method in (m for m in methods if m in R_GMEAN):
        np.testing.assert_allclose(
            _ours(x, method, sm, "gmean"), _r(case, method, "gmean"), rtol=1e-12, atol=R_GOLDEN_ATOL
        )
        checked += 1
    assert checked == sum(m in R_GMEAN for m in methods)


def _r_gmean(
    raw: np.ndarray, rows: list[str], method: str, sm: np.ndarray, options: dict
) -> np.ndarray:
    """R's gmean from R's bound: maxscost from the upper triangle, capped at two indels (OMdistance.cpp:23-34)."""
    indel = options["indel"]
    c = min(sm[np.triu_indices(4, k=1)].max(), 2 * indel)
    lengths = np.array([len(row) for row in rows], dtype=float)
    out = np.zeros_like(raw)
    for i in range(len(rows)):
        for j in range(len(rows)):
            m, n = lengths[i], lengths[j]
            if method == "TWED":
                bound = abs(n - m) * (options["nu"] + options["h"] + c) + 2 * (
                    c + options["nu"]
                ) * min(m, n)
            else:
                bound = abs(n - m) * indel + c * min(m, n)
            if raw[i, j] != 0:
                out[i, j] = 1 - (bound - raw[i, j]) / (2 * np.sqrt(m * indel) * np.sqrt(n * indel))
    return out


def test_reference_route_uses_the_same_bound() -> None:
    # refseq gives d(x_i, reference): d(BB, AA) = 6 and d(AA, AA) = 0; the bound is 6 (as above), so
    # maxdist is 1 and 0 (R: 3 and 0).
    x = _quiet(tx.seqdef, [list("BB"), list("AA")], alphabet=["A", "B"])
    to_second = np.asarray(
        _quiet(tx.seqdist, x, method="OM", subm=HAND, indel=2.0, norm="maxdist", refseq=2), float
    ).ravel()
    np.testing.assert_allclose(to_second, [1.0, 0.0], rtol=0, atol=1e-12)


@pytest.mark.parametrize("seed", range(200))
def test_random_asymmetric_costs_stay_within_the_bound(seed: int) -> None:
    """The plan's 200 seeded cases (25 had been run; batch-8 review I-2).

    OM and OMslen (both links) with an indel drawn in [0.25, 3]; OMspell and TWED at the
    parameters of ``EXTRA`` only (TX-REP-20261021-004, -005).
    """
    rng = np.random.default_rng(seed)
    sm = np.triu(rng.uniform(0.2, 1.0, (4, 4)), 1) + np.tril(rng.uniform(0.2, 4.0, (4, 4)), -1)
    rows = ["".join("ABCD"[k] for k in rng.integers(0, 4, rng.integers(3, 7))) for _ in range(6)]
    x = _sequences(rows)
    indel = float(rng.uniform(0.25, 3.0))
    runs = [
        ("OM", {"indel": indel}),
        ("OMslen", {"indel": indel, "link": "mean"}),
        ("OMslen", {"indel": indel, "link": "gmean"}),
        ("OMspell", EXTRA["OMspell"]),
        ("TWED", EXTRA["TWED"]),
    ]
    for method, options in runs:
        for norm in _corrected_norms(method):
            values = _quiet(tx.seqdist, x, method=method, subm=sm, norm=norm, **options)
            assert np.nanmax(values) <= 1 + 1e-12, (seed, method, options, norm)
        if method in R_GMEAN:
            raw = np.asarray(
                _quiet(tx.seqdist, x, method=method, subm=sm, norm="none", **options), float
            )
            gmean = np.asarray(
                _quiet(tx.seqdist, x, method=method, subm=sm, norm="gmean", **options), float
            )
            np.testing.assert_allclose(
                gmean, _r_gmean(raw, rows, method, sm, options), rtol=1e-12, atol=1e-12
            )
