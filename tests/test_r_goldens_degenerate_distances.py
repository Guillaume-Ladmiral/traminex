"""seqdist on degenerate inputs (TX-AUD-20261011-039, -040; fixture
``scripts/r_goldens/make_goldens_degenerate_distances.R``).

Transition-rate costs need at least two states and at least one transition. R's
seqcost(TRATE), seqtrate and every seqdist route that builds such costs stop on a
one-state alphabet or on sequences of length 1, and so do traminex's seqcost and
seqtrate; but seqdist's own OM/TRATE and DHD routes returned numbers (zeros for
OM/TRATE on one state). They now refuse as seqcost and seqtrate do.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "degenerate_distances.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
RAGGED = np.array([["A"] * 4, ["A", "A", None, None], ["A", "A", "A", None]], dtype=object)
ONE = tx.seqdef(RAGGED, alphabet=["A"])
ONE_FULL = tx.seqdef(np.array([["A", "A"], ["A", "A"]], dtype=object), alphabet=["A"])
LEN1 = tx.seqdef(np.array([["A"], ["B"], ["A"]], dtype=object), alphabet=["A", "B"])
ONE_STATE = "TRATE costs require at least two states"
NO_TRANSITION = "absolute lag must be smaller than the sequence width"
TRANSITION_RATE_CASES = {
    "om_trate_one_state": (lambda: tx.seqdist(ONE, method="OM", subm="TRATE"), ONE_STATE),
    "omstran_trate_one_state": (
        lambda: tx.seqdist(ONE, method="OMstran", subm="TRATE", otto=0.5),
        ONE_STATE,
    ),
    "seqcost_trate_one_state": (lambda: tx.seqcost(ONE, "TRATE"), ONE_STATE),
    "dhd_one_state": (lambda: tx.seqdist(ONE_FULL, method="DHD"), ONE_STATE),
    "dhd_trate_one_state": (lambda: tx.seqdist(ONE_FULL, method="DHD", subm="TRATE"), ONE_STATE),
    "om_trate_len1": (lambda: tx.seqdist(LEN1, method="OM", subm="TRATE"), NO_TRANSITION),
    "dhd_len1": (lambda: tx.seqdist(LEN1, method="DHD"), NO_TRANSITION),
    "dhd_trate_len1": (lambda: tx.seqdist(LEN1, method="DHD", subm="TRATE"), NO_TRANSITION),
    "ham_trate_len1": (lambda: tx.seqdist(LEN1, method="HAM", subm="TRATE"), NO_TRANSITION),
    "seqcost_trate_len1": (lambda: tx.seqcost(LEN1, "TRATE"), NO_TRANSITION),
    "seqtrate_len1": (lambda: tx.seqtrate(LEN1), NO_TRANSITION),
}


@pytest.mark.parametrize("case", sorted(TRANSITION_RATE_CASES))
def test_transition_rate_costs_are_refused_where_r_stops(case: str) -> None:
    assert GOLDEN.loc[case, "ok"] == "FALSE"
    call, message = TRANSITION_RATE_CASES[case]
    with pytest.raises(ValueError, match=message):
        call()


# --- R defects on degenerate inputs that traminex already corrects
# (TX-AUD-20261011-038, -039, -040). No code changes: these tests pin the
# definitions' values against R's recorded ones.

LENGTHS = np.array([4, 2, 3])  # AAAA, AA, AAA: one state, so LCS(i, j) = min(l_i, l_j)
TWO = tx.seqdef(RAGGED, alphabet=["A", "B"])


def _r_matrix(case: str) -> np.ndarray:
    values = np.array([float(v) for v in GOLDEN.loc[case, "value"].split()])
    side = int(round(np.sqrt(values.size)))
    return values.reshape(side, side, order="F")


def _lcs_definitions() -> dict[str, np.ndarray]:
    lcs = np.minimum.outer(LENGTHS, LENGTHS).astype(float)
    total = np.add.outer(LENGTHS, LENGTHS).astype(float)
    raw = total - 2 * lcs
    return {
        "gmean": 1 - lcs / np.sqrt(np.multiply.outer(LENGTHS, LENGTHS)),
        "maxdist": raw / total,
        "YujianBo": np.divide(
            2 * raw, total + raw, out=np.zeros_like(raw), where=(total + raw) > 0
        ),
    }


@pytest.mark.parametrize("norm", ["gmean", "maxdist", "YujianBo"])
def test_lcs_normalisations_follow_the_definition_on_one_state(norm: str) -> None:
    expected = _lcs_definitions()[norm]
    with pytest.warns(UserWarning, match="max substitution cost less or equal to 0"):
        one = tx.seqdist(ONE, method="LCS", norm=norm)
    np.testing.assert_allclose(one, expected, atol=1e-12, rtol=0)
    np.testing.assert_allclose(
        tx.seqdist(TWO, method="LCS", norm=norm), expected, atol=1e-12, rtol=0
    )
    # R gives the definition with a two-state alphabet ...
    two_r = _r_matrix(f"lcs_{norm.lower()}_two_states")
    np.testing.assert_allclose(two_r, expected, atol=1e-12, rtol=0)
    # ... but 1 off the diagonal for gmean and maxdist with a one-state alphabet: its 1x1
    # cost matrix (an R defect); YujianBo agrees.
    r_one = _r_matrix(f"lcs_{norm.lower()}_one_state")
    if norm == "YujianBo":
        np.testing.assert_allclose(r_one, expected, atol=1e-12, rtol=0)
    else:
        np.testing.assert_array_equal(r_one, 1 - np.eye(3))


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
def test_chi2_and_euclid_give_zero_where_r_stops(method: str) -> None:
    from traminex.sequences import STSList

    s3 = tx.seqdef(
        np.array([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"]], dtype=object),
        alphabet=["A", "B"],
    )
    one_row = STSList(
        data=s3.data[:1].copy(),
        alphabet=list(s3.alphabet),
        nr_code=s3.nr_code,
        void_code=s3.void_code,
    )
    # A sequence against itself, and identical state distributions, are at distance 0.
    np.testing.assert_array_equal(tx.seqdist(one_row, method=method), [[0.0]])
    np.testing.assert_array_equal(tx.seqdist(ONE_FULL, method=method), np.zeros((2, 2)))
    for case in (f"{method.lower()}_one_row", f"{method.lower()}_one_state"):
        assert GOLDEN.loc[case, "ok"] == "FALSE"
        assert "at least two dimensions" in GOLDEN.loc[case, "error"]


def test_omstran_follows_r_s_algorithm_where_r_crashes() -> None:
    """R's OMstran (seqdist-OMstran.R) on AAAAA against A, constant costs, otto 0.5.

    The transition alphabet is {A->A}: AAAAA has 5 transitions, A has 1. R scales the
    indel by indelrate = max(indel)/max(sm), or 0 when max(sm) = 0 (:73-76), then gives
    each transition the indel otto * indel * indelrate / max(indel) (same states: no
    transition part). One state: sm = [0], so every cost is 0 and the distance is 0.
    Two states: sm = [[0, 2], [2, 0]] and seqcost's indel 1, so indelrate = 1/2, the
    transition indel is 0.5 * 0.5 = 0.25 and the distance is 4 * 0.25 = 1. R stops on
    both (its loop for (i in 1:(alphabet_size - 1)), :126, runs 1:0 for one element).
    """
    five = np.array([["A"] * 5, ["A", None, None, None, None]], dtype=object)
    one = tx.seqdef(five, alphabet=["A"])
    two = tx.seqdef(five, alphabet=["A", "B"])
    with pytest.warns(UserWarning):
        np.testing.assert_array_equal(
            tx.seqdist(one, method="OMstran", subm="CONSTANT", otto=0.5), np.zeros((2, 2))
        )
    np.testing.assert_allclose(
        tx.seqdist(two, method="OMstran", subm="CONSTANT", otto=0.5),
        [[0.0, 1.0], [1.0, 0.0]],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    for case in ("omstran_constant_one_state", "omstran_constant_two_states"):
        assert GOLDEN.loc[case, "ok"] == "FALSE"
