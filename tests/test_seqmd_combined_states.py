"""Combined states of ``seqMD`` (audit 2026-09-22, slice S2).

S2-08: ``seqMD`` sorted the combined alphabet by Python code points, so with
mixed-case state names ``B+Y`` came before ``a+x``, both in the combined
sequence object and in the rows and columns of ``what="cost"``. R builds the
combined object with ``seqdef``'s default alphabet (``seqMD.R:197``), i.e.
``seqstatl``'s ``sort()`` in the session collation (ICU under a UTF-8
locale): ``a+x, B+x, B+Y, c+x, c+Y``. Distances are unaffected (they do not
depend on the order of the alphabet).

S2-11: with channels of unequal width, ``seqMD`` padded the shorter channel
with its nonresponse marker whenever that channel had missing values
(``fill_with_miss``), giving combined states such as ``*+y`` and ``*+%`` where
R has ``%+y`` and a void cell: ``seqMD.R:169-172`` pads with void, and the
``fill.with.miss`` replacement (lines 175-178) touches only the channel's own
columns.

Source of truth: live TraMineR 2.2-14 under ``LC_COLLATE=C.UTF-8``; values
in ``tests/r_goldens/audit20260922_s2_cases.txt``
(``scripts/r_goldens/make_goldens_audit20260922_s2.R``).
"""

from __future__ import annotations

import ctypes.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import STSList, alphabet, seqdef, seqMD

FIXTURE = Path(__file__).resolve().parent / "r_goldens" / "audit20260922_s2_cases.txt"

needs_icu = pytest.mark.skipif(
    ctypes.util.find_library("icui18n") is None,
    reason="R-compatible collation needs the system ICU library",
)


def _r(case: str, field: str) -> str:
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        name, key, value = line.split("\t")
        if (name, key) == (case, field):
            return value
    raise AssertionError(f"no R record {case}/{field}")


def _decoded_rows(x: STSList) -> list[list[str]]:
    return [[x._state_for_code(int(code)) for code in row] for row in x.data]


def _mixed_case_channels() -> list[STSList]:
    first = seqdef(
        pd.DataFrame({"a": ["a", "B", "c"], "b": ["B", "c", "a"]}),
        alphabet=["a", "B", "c"],
    )
    second = seqdef(
        pd.DataFrame({"a": ["x", "x", "Y"], "b": ["Y", "x", "x"]}),
        alphabet=["x", "Y"],
    )
    return [first, second]


@needs_icu
def test_combined_alphabet_follows_r_collation():
    combined = seqMD(_mixed_case_channels())

    assert alphabet(combined) == _r("seqMD_case_order", "alphabet").split("|")
    expected_rows = [row.split("|") for row in _r("seqMD_case_order", "cells").split(";")]
    assert _decoded_rows(combined) == expected_rows


@needs_icu
def test_cost_rows_and_columns_follow_r_collation():
    costs = seqMD(_mixed_case_channels(), sm="CONSTANT", cweight=[1, 3], what="cost")

    assert costs["alphabet"] == _r("seqMD_case_cost", "alphabet").split("|")
    expected = np.array([float(value) for value in _r("seqMD_case_cost", "sm").split("|")])
    size = len(costs["alphabet"])
    np.testing.assert_array_equal(costs["sm"], expected.reshape(size, size, order="F"))


def _unequal_width_channels() -> list[STSList]:
    shorter = seqdef(
        pd.DataFrame({"a": ["A", "B"], "b": [None, "A"], "c": ["B", "B"]}),
        alphabet=["A", "B"],
        right=None,
    )
    longer = seqdef(
        pd.DataFrame(
            {
                "a": ["x", "y"],
                "b": ["y", "y"],
                "c": ["x", "y"],
                "d": ["y", "x"],
                "e": ["x", None],
            }
        ),
        alphabet=["x", "y"],
    )
    return [shorter, longer]


@needs_icu
@pytest.mark.parametrize(
    "case, reverse",
    [("seqMD_width_fill", False), ("seqMD_width_fill_reversed", True)],
)
def test_shorter_channel_is_padded_with_void_not_missing(case, reverse):
    channels = _unequal_width_channels()
    if reverse:
        channels.reverse()

    with pytest.warns(UserWarning, match="different length"):
        combined = seqMD(channels)

    assert alphabet(combined) == _r(case, "alphabet").split("|")
    assert _decoded_rows(combined) == [row.split("|") for row in _r(case, "cells").split(";")]
