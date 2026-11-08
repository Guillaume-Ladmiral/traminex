"""HSPELL_to_STS gives its states the alphabet R's seqdef gives them (TX-AUD-20261010-031).

TraMineRextras::HSPELL_to_STS returns a character data frame and seqdef sorts
its alphabet as seqstatl does: numbers by value when every state is an R
number, R's collation otherwise. traminex sorted the states by code point
(1, 10, 11, 12, 2, 3, 9 and A, B, C, _x, a, b, c). Fixture:
``scripts/r_goldens/make_goldens_hspell_alphabet_order.R``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "hspell_alphabet_order.csv",
    dtype=str,
    keep_default_na=False,
)
STATUSES = {
    "numeric": ["1", "10", "11", "12", "2", "3", "9"],
    "text": ["_x", "a", "A", "b", "B", "c", "C"],
    "accents": ["e", "é", "f", "E", "ê", "d", "D"],
}


@pytest.mark.parametrize("case", sorted(STATUSES))
def test_hspell_alphabet_follows_r_seqdef(case: str) -> None:
    spells = pd.DataFrame(
        {"id": list("abcdefg"), "b1": [1] * 7, "e1": [2] * 7, "s1": STATUSES[case]}
    )
    seqs = tx.HSPELL_to_STS(spells, id_col="id", begin=["b1"], end=["e1"], status=["s1"])
    expected = GOLDEN[GOLDEN["case"] == case].sort_values("position")["state"].tolist()
    assert seqs.alphabet == expected
