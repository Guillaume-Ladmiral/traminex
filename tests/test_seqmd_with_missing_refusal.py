"""``seqMD(with_missing=False)`` refuses a channel with missing values (audit S2-12).

For ``what="cost"`` and ``what="diss"`` R stops when a channel declared
``with.missing = FALSE`` contains the nonresponse marker (``seqMD.R:216-223``:
"found missing values in channel i, set with.missing as TRUE for that
channel"). traminex only warned and went on to build costs that include the
missing state, so the user got numbers R refuses to produce. ``what="MDseq"``
returns before that check in R, so it keeps building the sequences.

Source of truth: ``SOURCES/TraMineR/R/seqMD.R`` and live TraMineR 2.2-14
(``reports/audit/full-codebase-audit-20260922/agents/S2/fix/s2_12_r.out``).
"""

from __future__ import annotations

import pandas as pd
import pytest

from traminex import STSList, seqdef, seqMD


def _channels() -> list[STSList]:
    first = seqdef(
        pd.DataFrame({"a": ["A", "A", "B"], "b": ["B", None, "B"], "c": ["B", "C", "A"]}),
        alphabet=["A", "B", "C"],
    )
    second = seqdef(
        pd.DataFrame({"a": ["x", "y", "y"], "b": ["x", "y", "x"], "c": ["y", "y", "x"]}),
        alphabet=["x", "y"],
    )
    return [first, second]


REFUSAL = "found missing values in channel 1, set with_missing as True for that channel"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"with_missing": False, "sm": "CONSTANT", "what": "cost"},
        {"with_missing": False, "method": "OM", "sm": "CONSTANT", "what": "diss"},
        {"with_missing": [False, False], "method": "HAM", "what": "diss"},
    ],
)
def test_declared_without_missing_but_missing_present_is_refused(kwargs):
    with pytest.warns(UserWarning, match="Bad with.missing value for domain 1"):
        with pytest.raises(ValueError, match=REFUSAL):
            seqMD(_channels(), **kwargs)


def test_mdseq_still_builds_the_sequences_with_a_warning():
    with pytest.warns(UserWarning, match="Bad with.missing value for domain 1"):
        combined = seqMD(_channels(), with_missing=False)
    assert isinstance(combined, STSList)


def test_declaring_the_missing_channel_true_is_accepted():
    costs = seqMD(_channels(), with_missing=[True, False], sm="CONSTANT", what="cost")
    assert costs["sm"].shape[0] == len(costs["alphabet"])
