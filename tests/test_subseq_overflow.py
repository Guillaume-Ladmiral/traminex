"""NMS, NMSMST and SVRspell where R's double sums overflow (TX-AUD-20261019-029, TX-AUD-20261011-047).

R computes distances between distinct sequences only, so identical sequences are at distance 0; traminex returned
NaN for them once the counts overflowed. Between different sequences R forms Ival + Jval - 2 Aval (NMSdistance.cpp)
in doubles and gives NaN (Inf - Inf) or Inf once a count or a duration power overflows, although the distance is
finite: an R defect traminex reproduced. traminex recomputes such a pair with power-of-two rescaling and gives
the definition's value; a distance beyond the double range is infinite.
``scripts/r_goldens/make_goldens_subseq_overflow.R`` records R's matrices (its NaN, Inf and 0).

Expected values come from the definitions, computed exactly with Python integers: the number of pairs of matching
subsequences of each length (an exact counter below; for one-state sequences of lengths m and n the sum over
lengths gives sqrt(C(2m,m) + C(2n,n) - 2 C(m+n,n))), and for single spells of lengths a and b raised to tpow,
NMSMST's shared time min(a, b) and SVRspell's product, so sqrt(|a - b|) and |a - b| for one state.
"""

from __future__ import annotations

import math
import warnings
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "subseq_overflow.csv", dtype=str
).set_index("case")


def _r_matrix(case: str, size: int) -> np.ndarray:
    values = [float(v) for v in GOLDEN.loc[case, "value"].split()]
    return np.array(values).reshape(size, size, order="F")


def _seqs(rows: list[str]) -> tx.STSList:
    width = max(len(row) for row in rows)
    data = [list(row) + [None] * (width - len(row)) for row in rows]
    return tx.seqdef(np.array(data, dtype=object), alphabet=["A", "B"], right="DEL")


def _matches(x: str, y: str) -> list[int]:
    """Exact number of pairs of equal subsequences of each length 1..min(len) (index tuples, as NMS counts)."""
    size = min(len(x), len(y))
    previous = [[1] * (len(y) + 1) for _ in range(len(x) + 1)]  # length 0
    counts = []
    for _ in range(size):
        current = [[0] * (len(y) + 1) for _ in range(len(x) + 1)]
        for i in range(1, len(x) + 1):
            for j in range(1, len(y) + 1):
                current[i][j] = current[i - 1][j] + current[i][j - 1] - current[i - 1][j - 1]
                if x[i - 1] == y[j - 1]:
                    current[i][j] += previous[i - 1][j - 1]
        counts.append(current[len(x)][len(y)])
        previous = current
    return counts


def _nms_square(x: str, y: str, weight: Fraction) -> Fraction:
    return weight * (sum(_matches(x, x)) + sum(_matches(y, y)) - 2 * sum(_matches(x, y)))


def _root(square: Fraction | int) -> float:
    """The square root of an exact non-negative rational, as a double (Inf beyond the double range)."""
    square = Fraction(square)
    scale = 2 ** max(
        0, (square.numerator.bit_length() - square.denominator.bit_length()) // 2 - 500
    )
    try:
        return math.sqrt(float(square / (scale * scale))) * scale
    except OverflowError:
        return math.inf


def _run(rows: list[str], **options) -> tuple[np.ndarray, list[warnings.WarningMessage]]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = np.asarray(tx.seqdist(_seqs(rows), **options), dtype=float)
    return result, [w for w in caught if issubclass(w.category, RuntimeWarning)]


def test_the_one_state_closed_form_is_the_exact_count() -> None:
    for m, n in ((6, 5), (9, 4)):
        exact = _nms_square("A" * m, "A" * n, Fraction(1))
        assert exact == math.comb(2 * m, m) + math.comb(2 * n, n) - 2 * math.comb(m + n, n)


def test_nms_counts_beyond_the_double_range() -> None:
    rows = ["A" * 515, "A" * 515, "A" * 514]
    actual, runtime_warnings = _run(rows, method="NMS")
    r = _r_matrix("nms_onestate", 3)
    assert r[0, 1] == 0 and not np.isfinite(r[0, 2]) and not np.isfinite(r[1, 2])
    assert actual[0, 1] == 0.0
    expected = _root(math.comb(1030, 515) + math.comb(1028, 514) - 2 * math.comb(1029, 514))
    np.testing.assert_allclose([actual[0, 2], actual[1, 2]], [expected, expected], rtol=1e-12)
    assert runtime_warnings == []


def test_nms_weighted_sums_beyond_the_double_range() -> None:
    rows = ["A" * 12, "A" * 12, "A" * 11, "B", "BA"]
    actual, runtime_warnings = _run(rows, method="NMS", kweights=[1e305] * 12)
    r = _r_matrix("nms_kweights", 5)
    weight = Fraction(1e305)
    assert r[0, 1] == 0 and actual[0, 1] == 0.0
    for i, j in ((0, 2), (0, 3), (0, 4), (2, 3), (2, 4)):
        assert not np.isfinite(r[i, j])
        assert actual[i, j] == pytest.approx(
            _root(_nms_square(rows[i], rows[j], weight)), rel=1e-12
        )
    # A pair R computes finitely is R's, bit for bit (Review Focus 1).
    assert np.isfinite(r[3, 4]) and actual[3, 4] == r[3, 4]
    assert runtime_warnings == []


@pytest.mark.parametrize(
    ("method", "tpow"), [("NMSMST", 250.0), ("SVRspell", 150.0), ("SVRspell", 250.0)]
)
def test_spell_durations_beyond_the_double_range(method: str, tpow: float) -> None:
    rows = ["A" * 50, "A" * 50, "A" * 40, "B" * 50, "AB", "BA"]
    actual, runtime_warnings = _run(rows, method=method, tpow=tpow)
    r = _r_matrix(f"{method.lower()}_tpow_{int(tpow)}", 6)
    a, b = 50 ** int(tpow), 40 ** int(tpow)
    if method == "NMSMST":  # shared time: min(a, b) for one state, 0 for two
        expected = {(0, 2): _root(a - b), (0, 3): _root(2 * a), (2, 3): _root(a + b)}
    else:  # product of the durations: one state (a - b)^2, two states a^2 + b^2
        expected = {
            (0, 2): _root((a - b) ** 2),
            (0, 3): _root(2 * a * a),
            (2, 3): _root(a * a + b * b),
        }
    assert r[0, 1] == 0 and actual[0, 1] == 0.0
    for (i, j), value in expected.items():
        assert not np.isfinite(r[i, j])
        assert actual[i, j] == pytest.approx(value, rel=1e-12)
    # The long spell dominates its distance to a short sequence.
    leading = _root(a) if method == "NMSMST" else _root(a * a)
    for i in (0, 2, 3):
        for j in (4, 5):
            assert not np.isfinite(r[i, j])
            dominant = leading if i != 2 else (_root(b) if method == "NMSMST" else _root(b * b))
            assert actual[i, j] == pytest.approx(dominant, rel=1e-12)
    assert np.isfinite(r[4, 5]) and actual[4, 5] == r[4, 5]
    assert runtime_warnings == []


PROX = np.array(
    [[1.0, 0.3, 0.2, 0.1], [0.3, 1.0, 0.3, 0.2], [0.2, 0.3, 1.0, 0.3], [0.1, 0.2, 0.3, 1.0]]
)


@pytest.mark.parametrize(
    "options",
    [
        dict(method="NMS"),
        dict(method="NMS", prox=PROX),
        dict(method="NMSMST", tpow=1.3),
        dict(method="SVRspell", tpow=1.3),
    ],
    ids=["NMS", "NMS-prox", "NMSMST", "SVRspell"],
)
def test_the_rescaled_path_rounds_as_the_plain_path(options: dict) -> None:
    """Weights times 2**1022 multiply every distance by exactly 2**511. Every pair then overflows and takes the
    rescaled path, so bit equality shows that path rounds as the plain one, i.e. as R's arithmetic would with an
    unlimited exponent range."""
    from test_r_goldens_spell_powers import X

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plain = np.asarray(tx.seqdist(X, kweights=[1.0] * 30, **options), dtype=float)
        scaled = np.asarray(tx.seqdist(X, kweights=[2.0**1022] * 30, **options), dtype=float)
    assert np.all(np.isfinite(plain)) and np.all(np.isfinite(scaled))
    np.testing.assert_array_equal(scaled, np.ldexp(plain, 511))
