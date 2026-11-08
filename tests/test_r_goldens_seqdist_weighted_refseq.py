"""seqdist refseq on a weighted object for CHI2 and OMstran (TX-AUD-20261010-003).

R accepts every reference form for every method (seqdist.R:65-98); traminex
refused them for the two methods that use sequence weights. Fixture:
``scripts/r_goldens/make_goldens_seqdist_weighted_refseq.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "seqdist_weighted_refseq.csv"
ROWS = np.array([list("AABBC"), list("ABBCC"), list("CCBAA"), list("BBCCA")], dtype=object)
WEIGHTS = [1.0, 2.0, 0.5, 1.5]


def _golden() -> dict[str, np.ndarray]:
    cells: dict[str, dict[tuple[int, int], float]] = defaultdict(dict)
    with GOLDEN.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cells[row["case"]][(int(row["row"]), int(row["col"]))] = float(row["value"])
    out = {}
    for case, values in cells.items():
        shape = (max(i for i, _ in values), max(j for _, j in values))
        matrix = np.empty(shape)
        for (i, j), value in values.items():
            matrix[i - 1, j - 1] = value
        out[case] = matrix
    return out


R = _golden()
SEQS = tx.seqdef(ROWS, alphabet=["A", "B", "C"], weights=WEIGHTS)
METHODS = {
    "CHI2": {"method": "CHI2"},
    "OMstran": {"method": "OMstran", "subm": "CONSTANT", "indel": 1.0, "otto": 0.5},
}
REFERENCES = {
    "sets": [[1, 2], [3, 4]],
    "index2": 2,
    "mostfrequent": 0,
    # R's seqs[2, ]: the second sequence, carrying its weight.
    "sequence2": tx.seqdef(ROWS[1:2], alphabet=["A", "B", "C"], weights=[2.0]),
}
CASES = [
    (method, reference, weighted)
    for method in METHODS
    for reference in REFERENCES
    for weighted in (True, False)
]


def test_every_generated_case_is_exercised():
    expected = {f"{m}_{r}_{'weighted' if w else 'unweighted'}" for m, r, w in CASES}
    assert expected == set(R)


@pytest.mark.parametrize(("method", "reference", "weighted"), CASES)
def test_weighted_object_reference_matches_r(method: str, reference: str, weighted: bool):
    actual = np.asarray(
        tx.seqdist(SEQS, refseq=REFERENCES[reference], weighted=weighted, **METHODS[method]),
        dtype=float,
    )
    expected = R[f"{method}_{reference}_{'weighted' if weighted else 'unweighted'}"]
    if expected.shape[1] == 1:
        expected = expected[:, 0]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=R_GOLDEN_ATOL)


def test_explicit_weights_still_conflict_with_object_weights():
    with pytest.raises(ValueError, match="conflicts with the weights"):
        tx.seqdist(SEQS, method="CHI2", weights=WEIGHTS)
