"""Validation gaps in the low-level distance API (2026-09-22 audit).

These surfaces are traminex's own (the low-level ``pairwise``, the array forms
of ``seqLLCS``/``seqLLCP``), so the standard is internal consistency with the
R-facing ``seqdist`` and with each other, not an R value.
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx
from traminex.distances import pairwise

X = np.array([[0, 1, 2, 1], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 2, 0]])


@pytest.mark.parametrize("bad", [[-1], [0.7], [4]])
def test_pairwise_ref_sets_validates_indices_like_seqdist(bad):
    """S1a-07: a negative index wrapped, a fractional one was truncated."""
    with pytest.raises((ValueError, IndexError), match="ref_sets"):
        pairwise(X, method="LCS", ref_sets=(bad, [1]))


def test_seqllcs_refuses_a_multi_row_array_as_it_refuses_a_multi_row_object():
    """S1a-09: a 2-D array silently used its first row."""
    with pytest.raises(ValueError, match="exactly one sequence"):
        tx.seqLLCS(X[:2], X[2])
    with pytest.raises(ValueError, match="exactly one sequence"):
        tx.seqLLCP(X[2], X[:2])
    assert tx.seqLLCS(X[:1], X[2]) == tx.seqLLCS(X[0], X[2])


@pytest.mark.parametrize("method", ["OMLOC", "OMSLEN", "OMSPELL", "TWED", "OMSTRAN"])
def test_a_cost_method_string_is_refused_by_name_at_the_low_level(method):
    """S1a-14: 'could not convert string to float: TRATE' named no remedy."""
    extra = {"nu": 0.5, "h": 0.5} if method == "TWED" else {}
    if method == "OMSTRAN":
        extra["otto"] = 0.5
    with pytest.raises(ValueError, match="seqdist"):
        pairwise(X, method=method, subm="TRATE", **extra)


def test_state_codes_beyond_32_bits_are_refused_instead_of_wrapping():
    """S1b-11: code 2**32 became 0 in an int32 cast, so seqdist returned 0
    between different sequences."""
    wide = np.array([[0, 2**32], [0, 0]], dtype=np.int64)
    with pytest.raises(ValueError, match="32-bit"):
        tx.seqdist(wide, method="HAM")


def test_dense_pairwise_dhd_refuses_a_misspelt_norm():
    """S1b-09: the dense default DHD accepted norm="bogus" silently. Valid
    spellings now normalise by TraMineR's bound, as seqdist does
    (TX-AUD-20261010-023; R values in test_pairwise_dhd_routes.py)."""
    normalised = pairwise(X, method="DHD", norm="auto")
    assert not np.array_equal(normalised, pairwise(X, method="DHD"))
    assert np.all((normalised >= 0) & (normalised <= 1))
    with pytest.raises(ValueError, match="norm"):
        pairwise(X, method="DHD", norm="bogus")


def test_legacy_dhd_smoothing_handles_zero_width_input():
    """S1b-15: np.pad raised on a zero-width frequency table."""
    empty = np.zeros((3, 0), dtype=int)
    result = pairwise(empty, method="DHD", variant="mean", smooth_window=3)
    np.testing.assert_array_equal(result, np.zeros((3, 3)))
