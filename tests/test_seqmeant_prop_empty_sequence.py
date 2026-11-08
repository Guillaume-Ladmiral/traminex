"""``seqmeant(prop=True)`` is NaN when a sequence has no valid cell.

``seqmeant.R:10-26`` averages ``seqistatd(prop = TRUE)``, which divides each
row by its own total (``seqistatd.R:31-33``).  A sequence without any valid
cell (all missing, with ``with.missing=FALSE``) has a ``0/0 = NaN`` row, so
every weighted mean is ``NaN`` -- even when that sequence's weight is zero,
since ``NaN * 0`` is ``NaN``.  traminex divided with a zero fallback, so it
returned finite means (``0.111 0.333 0.222`` where R gives ``NaN NaN NaN``),
contradicting its own ``seqistatd(prop=True)`` (2026-09-22 audit S3-08).  The
marker-padded rows of the result stay zero, as documented.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-08), live TraMineR 2.2-14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.sequences import STSList

GOLDENS = Path(__file__).parent / "r_goldens" / "audit20260922_s3a.txt"
ROWS = [["A", None, "B", "B"], [None, None, None, None], ["B", None, "C", "C"]]
STATES = ["A", "B", "C"]


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDENS.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, *_ = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _r_values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", R[case]
    return np.array(payload.split(" "), dtype=float)


def _mm(weights=None) -> STSList:
    return tx.seqdef(np.array(ROWS, dtype=object), alphabet=STATES, right=None, weights=weights)


def _external(x: STSList) -> STSList:
    assert x.alphabet == ["A", "B", "C", "*"] and x.nr_code == 3 and x.void_code is None
    weights = None if x.weights is None else x.weights.copy()
    return STSList(x.data.copy(), alphabet=STATES, nr="*", void="%", nr_code=3, weights=weights)


def _close(actual, case: str) -> None:
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), _r_values(case), rtol=0, atol=1e-15, equal_nan=True
    )


@pytest.fixture(params=["in_alphabet", "external"])
def rep(request):
    return (lambda x: x) if request.param == "in_alphabet" else _external


def test_prop_means_are_nan_as_in_r(rep):
    result = tx.seqmeant(rep(_mm()), prop=True)
    _close(result.loc[STATES, "Mean"], "s3_08_meant_prop")
    # The padded nr row is not part of R's result and stays zero.
    assert result.loc["*", "Mean"] == 0.0


def test_prop_serr_columns_are_nan_as_in_r(rep):
    result = tx.seqmeant(rep(_mm()), prop=True, serr=True)
    _close(
        result.loc[STATES, ["Mean", "Var", "Stdev", "SE"]].to_numpy().ravel(order="F"),
        "s3_08_meant_prop_serr",
    )


def test_a_zero_weight_does_not_hide_the_nan(rep):
    result = tx.seqmeant(rep(_mm(weights=[2, 0, 1])), prop=True)
    _close(result.loc[STATES, "Mean"], "s3_08_meant_prop_zero_weight")


def test_with_missing_and_count_means_are_unchanged(rep):
    with_missing = tx.seqmeant(rep(_mm()), prop=True, include_missing=True)
    _close(with_missing.loc[STATES + ["*"], "Mean"], "s3_08_meant_prop_wm")
    counts = tx.seqmeant(rep(_mm()))
    _close(counts.loc[STATES, "Mean"], "s3_08_meant_counts")


def test_seqmeant_prop_agrees_with_seqistatd_prop():
    x = _mm()
    proportions = np.asarray(tx.seqistatd(x, prop=True), dtype=float)
    np.testing.assert_allclose(
        tx.seqmeant(x, prop=True)["Mean"].to_numpy(),
        proportions.mean(axis=0),
        rtol=0,
        atol=1e-15,
        equal_nan=True,
    )
