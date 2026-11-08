"""``kweights`` for NMS, NMSMST and SVRspell follow ``seqdist.R:265-273``.

R recycles a scalar over the sequence length, pads a short vector with zeros,
refuses a negative (or NaN) weight, and an infinite weight gives NaN
distances. traminex crashed on a scalar (``IndexError``), accepted negative
weights, and clamped a NaN squared distance to 0 (2026-09-22 audit S1b-06).
Expected values: live TraMineR 2.2-14 via
``scripts/r_goldens/make_goldens_seqdist_kweights.R``.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDENS = Path(__file__).parent / "r_goldens"
SEQUENCES = [list("AABCC"), list("ABBCA"), list("CCAAB"), list("BACCC")]
KWEIGHTS = {
    "scalar": 0.5,
    "short": [1.0, 0.5],
    "full": [1.0, 0.8, 0.6, 0.4, 0.2],
    "infinite": [math.inf, 1.0, 1.0, 1.0, 1.0],
    "negative": [1.0, -0.5, 1.0, 1.0, 1.0],
    "not_a_number": [math.nan, 1.0, 1.0, 1.0, 1.0],
}
METHODS = {"NMS": "NMS", "NMSMST": "NMSMST", "SVRspell": "SVRSPELL"}


def _r_matrices() -> dict[tuple[str, str], np.ndarray]:
    cells: dict[tuple[str, str], dict] = defaultdict(dict)
    with (GOLDENS / "seqdist_kweights_cases.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["method"], row["case"])
            cells[key][(int(row["row"]), int(row["col"]))] = float(row["value"])
    out = {}
    for key, values in cells.items():
        matrix = np.empty((4, 4))
        for (i, j), value in values.items():
            matrix[i, j] = value
        out[key] = matrix
    return out


R = _r_matrices()
R_ERRORS = {
    (line.split("\t")[0], line.split("\t")[1])
    for line in (GOLDENS / "seqdist_kweights_outcomes.txt").read_text().splitlines()
    if not line.startswith("#") and line.split("\t")[2] == "error"
}


def _sequences():
    return tx.seqdef(np.array(SEQUENCES, dtype=object), alphabet=["A", "B", "C"])


@pytest.mark.parametrize("r_method", sorted(METHODS))
@pytest.mark.parametrize("case", ["scalar", "short", "full", "infinite"])
def test_accepted_kweights_match_r(r_method, case):
    result = tx.seqdist(_sequences(), method=METHODS[r_method], kweights=KWEIGHTS[case])

    np.testing.assert_allclose(
        np.asarray(result, dtype=float), R[(r_method, case)], rtol=1e-12, atol=1e-12, equal_nan=True
    )


@pytest.mark.parametrize("r_method", sorted(METHODS))
@pytest.mark.parametrize("case", ["negative", "not_a_number"])
def test_refused_kweights_are_refused(r_method, case):
    assert (r_method, case) in R_ERRORS
    with pytest.raises(ValueError, match="kweights"):
        tx.seqdist(_sequences(), method=METHODS[r_method], kweights=KWEIGHTS[case])
