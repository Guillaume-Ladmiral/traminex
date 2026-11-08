"""seqedist with event sequences that have no events (TX-AUD-20261010-037).

TraMineR's seqedist.cpp leaves its column index at 1 when the row sequence is
empty, so the distance from an empty sequence to every later sequence is 0 and
the matrix depends on the row order; for two empty sequences its Yujian-Bo
normalisation divides 0 by 0. Both are R defects (ledger seqedist row,
DIVERGENCE). By definition the distance from an empty event sequence is the
cost of inserting every event of the other one, and two sequences at raw
distance 0 are at normalised distance 0, as TraMineR's own normaliser returns
(distancecalculator.h: ``if (rawdist==0) return 0``).

Fixture: ``scripts/r_goldens/make_goldens_seqedist_empty.R`` (R's values, kept
so the divergence stays visible). Cells between two non-empty sequences are not
affected by the defect and are compared with R.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "seqedist_empty.csv"


def _golden():
    cells: dict[tuple[str, str], dict[tuple[int, int], float]] = defaultdict(dict)
    sequences: dict[str, list[str]] = defaultdict(list)
    with GOLDEN.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["norm"] == "sequences":
                sequences[row["case"]].append(row["value"])
                continue
            cells[(row["case"], row["norm"])][(int(row["row"]), int(row["col"]))] = float(
                row["value"]
            )
    matrices = {}
    for key, values in cells.items():
        size = max(i for i, _ in values)
        matrix = np.empty((size, size))
        for (i, j), value in values.items():
            matrix[i - 1, j - 1] = value
        matrices[key] = matrix
    return matrices, dict(sequences)


R, R_SEQUENCES = _golden()

RECORD = tx.seqecreate(
    id=[1, 2, 2, 2, 3, 3, 4],
    timestamp=[5, 0, 2, 6, 1, 4, 3],
    event=["END", "X", "Y", "END", "Y", "X", "END"],
    end_event="END",
)
FRESH = tx.seqecreate(
    id=[1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 5],
    timestamp=[0, 1, 3, 2, 0, 2, 4, 1, 1, 1, 3],
    event=["A", "B", "END", "END", "B", "A", "END", "END", "A", "B", "END"],
    end_event="END",
)
REORDERED = tx.seqecreate(
    id=[1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5],
    timestamp=[2, 1, 0, 1, 3, 0, 2, 4, 1, 1, 3],
    event=["END", "END", "A", "B", "END", "B", "A", "END", "A", "B", "END"],
    end_event="END",
)
CASES = {
    "record": (RECORD, [1.0, 1.0, 1.0], 0.1),
    "fresh": (FRESH, [1.0, 2.0, 0.5], 0.2),
    "reordered": (REORDERED, [1.0, 2.0, 0.5], 0.2),
}
#: Cost of inserting every event of each sequence (idcost of X and Y are 1 in
#: the record case; A costs 1 and B 2 in the fresh cases). Empty sequences: 0.
INSERTION = {
    "record": [0.0, 2.0, 2.0, 0.0],
    "fresh": [3.0, 0.0, 3.0, 0.0, 3.0],
    "reordered": [0.0, 0.0, 3.0, 3.0, 3.0],
}


def _expected(case: str, norm: str) -> np.ndarray:
    """R's value between two non-empty sequences; the definition elsewhere."""
    cost = INSERTION[case]
    expected = R[(case, norm)].copy()
    for i, ci in enumerate(cost):
        for j, cj in enumerate(cost):
            if i == j or (ci > 0 and cj > 0):
                continue
            raw = ci + cj  # one side is empty: insert the other's events
            maxcost = ci + cj
            if norm == "none" or raw == 0:
                expected[i, j] = raw
            elif norm == "max":
                expected[i, j] = raw / maxcost
            else:
                expected[i, j] = 2 * raw / (maxcost + raw)
    return expected


def test_fixture_sequences_are_the_cases_tested() -> None:
    for case, (eseq, _, _) in CASES.items():
        assert eseq.as_character() == R_SEQUENCES[case]


@pytest.mark.parametrize("norm", ["none", "max", "YujianBo"])
@pytest.mark.parametrize("case", sorted(CASES))
def test_empty_sequence_distance_is_insertion_cost(case: str, norm: str) -> None:
    eseq, idcost, vparam = CASES[case]
    got = tx.seqedist(eseq, idcost=idcost, vparam=vparam, norm=norm)
    np.testing.assert_allclose(got, _expected(case, norm), rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_array_equal(got, got.T)


def test_result_does_not_depend_on_row_order() -> None:
    order = [1, 3, 0, 2, 4]  # fresh rows as the reordered case lists them
    fresh = tx.seqedist(FRESH, idcost=[1.0, 2.0, 0.5], vparam=0.2)
    reordered = tx.seqedist(REORDERED, idcost=[1.0, 2.0, 0.5], vparam=0.2)
    np.testing.assert_allclose(reordered, fresh[np.ix_(order, order)], rtol=0, atol=1e-12)


@pytest.mark.parametrize("norm", ["none", "max", "YujianBo"])
def test_two_empty_sequences_are_at_distance_zero(norm: str) -> None:
    got = tx.seqedist(RECORD, idcost=[1.0, 1.0, 1.0], vparam=0.1, norm=norm)
    assert got[0, 3] == got[3, 0] == 0.0


def test_r_values_recorded_and_differ() -> None:
    # R's first empty sequence is at distance 0 from every later sequence ...
    assert np.all(R[("record", "none")][0, :] == 0.0)
    assert R[("fresh", "none")][3, 4] == 0.0
    # ... and its Yujian-Bo distance between the two empty sequences is NaN.
    assert np.isnan(R[("record", "YujianBo")][0, 3])
