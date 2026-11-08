"""Low-level CHI2/EUCLID ``weights`` (TX-AUD-20261011-008, TX-AUD-20261010-067).

The low-level ``pairwise`` checked only the length of ``weights``: NaN or infinite
weights gave an all-zero CHI2 matrix and a negative weight altered it, where
``seqdist`` refuses such weights. EUCLID read the weights too, although R never
weights EUCLID (seqdist-CHI2.R:7 sets ``weighted <- FALSE``): a zero weight
changed its distances. Fixture: ``scripts/r_goldens/make_goldens_chi2_euclid_weights.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx
from traminex.distances import pairwise

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "chi2_euclid_weights.csv"
WEIGHTS = {
    "unit": [1.0, 1.0, 1.0, 1.0],
    "uneven": [1.0, 5.0, 0.2, 3.0],
    "one_zero": [0.0, 1.0, 1.0, 1.0],
    "all_zero": [0.0, 0.0, 0.0, 0.0],
}
SEQS = tx.seqdef(
    np.array([list(s) for s in ("AABBC", "ABBCC", "CCBAA", "BBCCA")], dtype=object),
    alphabet=["A", "B", "C"],
)
CODES = np.asarray(SEQS.data, dtype=np.int32)


def _golden() -> dict[str, np.ndarray]:
    cells: dict[str, dict[tuple[int, int], float]] = defaultdict(dict)
    with GOLDEN.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cells[row["case"]][(int(row["row"]), int(row["col"]))] = float(row["value"])
    out = {}
    for case, values in cells.items():
        matrix = np.empty((4, 4))
        for (i, j), value in values.items():
            matrix[i - 1, j - 1] = value
        out[case] = matrix
    return out


R = _golden()


@pytest.mark.parametrize("name", sorted(WEIGHTS))
def test_low_level_chi2_matches_r(name: str) -> None:
    actual = pairwise(CODES, method="CHI2", weights=np.array(WEIGHTS[name]))
    np.testing.assert_allclose(actual, R[f"CHI2_{name}"], rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize("name", sorted(WEIGHTS))
def test_low_level_euclid_ignores_weights(name: str) -> None:
    weighted = pairwise(CODES, method="EUCLID", weights=np.array(WEIGHTS[name]))
    np.testing.assert_allclose(weighted, R[f"EUCLID_{name}"], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_array_equal(weighted, pairwise(CODES, method="EUCLID"))


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
@pytest.mark.parametrize(
    ("bad", "message"),
    [
        ([np.nan, 1.0, 1.0, 1.0], "finite"),
        ([np.inf, 1.0, 1.0, 1.0], "finite"),
        ([-1.0, 1.0, 1.0, 1.0], "non-negative"),
        ([1.0, 1.0, 1.0], "length"),
    ],
)
def test_low_level_weights_are_validated(method: str, bad: list[float], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        pairwise(CODES, method=method, weights=np.array(bad))
