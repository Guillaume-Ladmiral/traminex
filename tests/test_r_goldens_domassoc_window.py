"""seqdomassoc counts windowed co-occurrences within each person (TX-AUD-20261011-009).

With ``wrange``, R's seqxtabs.win windows the column-major unlist of each domain
in blocks of ncol/ncases as if each block were one case's sequence
(seqdomassoc.R:170,210,277-301): it counts co-occurrences across different
people, ignores the weights and uses the alphabet only with with.missing = TRUE.
traminex copied that. The corrected count takes, for each person, the state at
position p in the first domain against the states at p + wrange[0] ..
p + wrange[1] in the second, weighted by the person's weight. The definition's
values come from R's position-wise seqdomassoc on sequences listing every such
pair; R's own windowed values are kept to show the difference. Fixture:
``scripts/r_goldens/make_goldens_domassoc_window.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import load_golden

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "domassoc_window.csv",
    dtype=str,
    keep_default_na=False,
)


def _md_domains():
    alphabet = ["1", "2", "3"]
    domains = []
    for k in (1, 2, 3):
        codes = load_golden(f"md_input_s{k}").astype(int)
        domains.append(
            tx.seqdef(
                np.array([[alphabet[c] for c in row] for row in codes], dtype=object),
                alphabet=alphabet,
            )
        )
    return domains


TINY = [
    tx.seqdef([["A", "A", "A"], ["B", "B", "B"]], alphabet=["A", "B"]),
    tx.seqdef([["X", "X", "X"], ["Y", "Y", "Y"]], alphabet=["X", "Y"]),
]
MISSING = [
    tx.seqdef([["A", "A", "B"], ["B", None, "B"], [None, "B", "A"]], alphabet=["A", "B"]),
    tx.seqdef([["X", "X", "Y"], [None, "Y", "X"], ["Y", None, "Y"]], alphabet=["X", "Y"]),
]
CASES = {
    "tiny": (TINY, (0, 1), False),
    "md_wrange": (_md_domains(), (-1, 1), False),
    "missing": (MISSING, (0, 1), True),
}


def _r(case: str, pair: str, kind: str, section: str) -> dict[str, float]:
    rows = GOLDEN[
        (GOLDEN["case"] == case)
        & (GOLDEN["pair"] == pair)
        & (GOLDEN["kind"] == kind)
        & (GOLDEN["section"] == section)
    ]
    return {name: float(value) for name, value in zip(rows["name"], rows["value"])}


def _run(case: str):
    domains, wrange, with_missing = CASES[case]
    return tx.seqdomassoc(domains, wrange=wrange, with_missing=with_missing, cross_table=True)


def _table(result, pair: str) -> dict[str, float]:
    table = result.attrs["cross_tables"][pair]
    return {
        f"{row}|{column}": float(table.loc[row, column])
        for row in table.index
        for column in table.columns
    }


@pytest.mark.parametrize("case", sorted(CASES))
def test_windowed_tables_and_statistics_follow_the_definition(case: str) -> None:
    result = _run(case)
    for pair in result.index:
        assert _table(result, pair) == _r(case, pair, "definition", "table")
        expected = _r(case, pair, "definition", "stat")
        ours = result.loc[pair, list(expected)].to_numpy(dtype=float)
        np.testing.assert_allclose(ours, list(expected.values()), rtol=0, atol=1e-10)


@pytest.mark.parametrize("case", sorted(CASES))
def test_r_windowed_tables_mix_people(case: str) -> None:
    pairs = sorted(set(GOLDEN.loc[GOLDEN["case"] == case, "pair"]))
    assert any(
        _r(case, pair, "r_window", "table") != _r(case, pair, "definition", "table")
        for pair in pairs
    )


def test_weights_scale_each_persons_count() -> None:
    # Person 1 (weight 1) has 5 window pairs A-X, person 2 (weight 5) 5 pairs B-Y.
    weighted = [
        tx.seqdef([["A", "A", "A"], ["B", "B", "B"]], alphabet=["A", "B"], weights=[1, 5]),
        tx.seqdef([["X", "X", "X"], ["Y", "Y", "Y"]], alphabet=["X", "Y"], weights=[1, 5]),
    ]
    result = tx.seqdomassoc(weighted, wrange=(0, 1), cross_table=True)
    assert _table(result, "Dom1_with_Dom2") == {"A|X": 5.0, "A|Y": 0.0, "B|X": 0.0, "B|Y": 25.0}
    unweighted = tx.seqdomassoc(weighted, wrange=(0, 1), cross_table=True, weighted=False)
    assert _table(unweighted, "Dom1_with_Dom2") == {"A|X": 5.0, "A|Y": 0.0, "B|X": 0.0, "B|Y": 5.0}
