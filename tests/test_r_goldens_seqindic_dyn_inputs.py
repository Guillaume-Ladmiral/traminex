"""seqindic_dyn checks indic as R's seqindic.dyn does (TX-AUD-20261010-051).

R stops for an indicator group name ("basic", "all", "complexity", ...) and
uses only the first value of a vector indic, with a message
(seqindic-dyn.R:10-22). traminex returned the first indicator's values for a
group name, raised TypeError for a list and IndexError for a name seqindic
does not know ("visitp", for which R fails inside its own code). Fixture:
``scripts/r_goldens/make_goldens_seqindic_dyn_inputs.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqindic_dyn_inputs.csv",
    dtype=str,
    keep_default_na=False,
)
STRINGS = ["AABBCCAABB", "ABCABCABCA", "AAAAABBBBB", "CCBBAACCBB", "ABBBCCCAAA", "BBBBBBCCCC"]
SEQS = tx.seqdef(
    np.array([list(s) for s in STRINGS], dtype=object),
    alphabet=["A", "B", "C"],
    id=[f"s{k}" for k in range(1, 7)],
)


def _r(case: str, kind: str) -> str:
    rows = GOLDEN[(GOLDEN["case"] == case) & (GOLDEN["kind"] == kind)]
    return rows["text"].iloc[0]


@pytest.mark.parametrize("group", ["basic", "all", "complexity"])
def test_an_indicator_group_name_is_refused_with_r_message(group: str) -> None:
    with pytest.raises(ValueError, match=_r(group, "error")):
        tx.seqindic_dyn(SEQS, indic=group)


def test_a_vector_indic_uses_its_first_value_with_r_message() -> None:
    message = _r("cplx,trans", "warning").removeprefix("[!!] ")
    with pytest.warns(UserWarning, match=message):
        vector = tx.seqindic_dyn(SEQS, indic=["cplx", "trans"])
    assert _r("cplx,trans", "equals_first") == "TRUE"
    np.testing.assert_array_equal(vector, tx.seqindic_dyn(SEQS, indic="cplx"))


def test_an_unknown_indicator_is_named() -> None:
    assert _r("visitp", "error") == "replacement has length zero"  # R fails internally
    with pytest.raises(ValueError, match="visitp"):
        tx.seqindic_dyn(SEQS, indic="visitp")


def test_the_result_has_one_column_per_r_window() -> None:
    r_columns = _r("cplx", "colnames").split("|")
    result = tx.seqindic_dyn(SEQS, indic="cplx")
    assert result.shape == (6, len(r_columns))
