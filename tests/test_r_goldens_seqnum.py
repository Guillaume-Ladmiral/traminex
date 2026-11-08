"""seqnum against live R (TX-AUD-20261010-021).

The ledger row for seqnum was EXACT/FULL, but no generator called R's seqnum.
R renumbers the states 0..n-1 in alphabet order; with ``with.missing = TRUE``
the missing code becomes the number n while the alphabet stays 0..n-1, and the
void symbol stays (seqnum.R:9-37). traminex keeps the marker slots at the end
of its alphabet and labels, so the real states are compared. Fixture:
``scripts/r_goldens/make_goldens_seqnum.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqnum_cases.csv",
    dtype=str,
    keep_default_na=False,
)


def _rows(strings: list[str]) -> np.ndarray:
    return np.array([list(s) for s in strings], dtype=object)


MISSING = ["A*BC", "*CBA", "BAC*", "AB**"]
CASES = {
    "unsorted_alphabet": (
        lambda: tx.seqdef(
            _rows(["ABCA", "CCBA", "BACC"]),
            alphabet=["C", "A", "B"],
            labels=["lab-C", "lab-A", "lab-B"],
        ),
        False,
    ),
    "missing_void": (
        lambda: tx.seqdef(_rows(MISSING), alphabet=["A", "B", "C"], missing="*"),
        False,
    ),
    "missing_void_with_missing": (
        lambda: tx.seqdef(_rows(MISSING), alphabet=["A", "B", "C"], missing="*"),
        True,
    ),
    "numeric_states": (lambda: tx.seqdef(_rows(["1231", "3321", "2133"])), False),
}


def _r(case: str, quantity: str) -> list[str]:
    rows = GOLDEN[(GOLDEN["case"] == case) & (GOLDEN["quantity"] == quantity)]
    return rows.sort_values("position", key=lambda c: c.astype(int))["value"].tolist()


@pytest.mark.parametrize("case", sorted(CASES))
def test_seqnum_matches_r(case: str) -> None:
    build, with_missing = CASES[case]
    seqs = build()
    assert [str(state) for state in seqs.alphabet][: len(_r(case, "input_alphabet"))] == _r(
        case, "input_alphabet"
    )
    result = tx.seqnum(seqs, with_missing=with_missing)
    markers = {result.nr, result.void}
    real = [str(state) for state in result.alphabet if state not in markers]
    assert real == _r(case, "alphabet")
    assert [str(result.nr)] == _r(case, "nr")
    assert [str(result.void)] == _r(case, "void")
    assert [str(label) for label in result.labels][: len(real)] == _r(case, "labels")
    canonical = result._canonicalized()
    cells = ["|".join(str(canonical.alphabet[code]) for code in row) for row in canonical.data]
    assert cells == _r(case, "cells")
