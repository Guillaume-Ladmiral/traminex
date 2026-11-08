"""Low-level ``pairwise(method="DHD")`` validates costs and normalises as seqdist does.

TX-AUD-20261010-012: the explicit-cost route accepted negative, NaN and
non-zero-diagonal costs and raised a raw IndexError for an undersized array,
where R's checkcost refuses all four (TraMineR-checkcost.R:40-75).
TX-AUD-20261010-023: the default route ignored ``norm`` and the explicit route
divided by the sequence length, where R divides HAM and DHD distances by the sum
over positions of the largest cost (seqdist.R:645-649, 771-776).
Fixture: ``scripts/r_goldens/make_goldens_dhd_explicit_norms.R``.
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

GOLDENS = Path(__file__).resolve().parent / "r_goldens"
NORMS = ("none", "maxlength", "gmean", "maxdist", "YujianBo", "auto")
SM5 = np.abs(np.subtract.outer(np.arange(4), np.arange(4))).astype(float)
ARRAY = np.stack([SM5 * t for t in range(1, 7)])  # time-first (T, K, K)
SEQS = tx.seqdef(
    np.array([list(s) for s in ("ABBCDA", "BBCCDD", "DACBBA", "AADDCC", "CBADCB")], dtype=object),
    alphabet=["A", "B", "C", "D"],
)
CODES = np.asarray(SEQS.data, dtype=np.int32)


def _golden() -> dict[str, np.ndarray]:
    cells: dict[str, dict[tuple[int, int], float]] = defaultdict(dict)
    with (GOLDENS / "dhd_explicit_norms.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cells[row["case"]][(int(row["row"]), int(row["col"]))] = float(row["value"])
    out = {}
    for case, values in cells.items():
        matrix = np.empty((5, 5))
        for (i, j), value in values.items():
            matrix[i - 1, j - 1] = value
        out[case] = matrix
    return out


R = _golden()
R_ERRORS = dict(
    line.split("|", 1)
    for line in (GOLDENS / "dhd_explicit_norms_errors.txt").read_text(encoding="utf-8").splitlines()
)


def _low_level_norm(norm: str) -> str | None:
    return None if norm == "none" else norm


@pytest.mark.parametrize("norm", NORMS)
def test_low_level_explicit_dhd_matches_r(norm: str) -> None:
    actual = pairwise(CODES, method="DHD", subm=ARRAY, norm=_low_level_norm(norm))
    np.testing.assert_allclose(actual, R[f"explicit_{norm}"], rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize("norm", NORMS)
def test_low_level_default_dhd_matches_r(norm: str) -> None:
    actual = pairwise(CODES, method="DHD", norm=_low_level_norm(norm))
    np.testing.assert_allclose(actual, R[f"default_{norm}"], rtol=0, atol=R_GOLDEN_ATOL)


def test_low_level_static_dhd_matrix_uses_r_ham_bound() -> None:
    actual = pairwise(CODES, method="DHD", subm=SM5, norm="maxlength")
    np.testing.assert_allclose(actual, R["ham_sm5_maxlength"], rtol=0, atol=R_GOLDEN_ATOL)


def _invalid(kind: str) -> np.ndarray:
    costs = ARRAY.copy()
    if kind == "negative":
        costs[0, 0, 1] = costs[0, 1, 0] = -1.0
    elif kind == "nan":
        costs[0, 0, 1] = costs[0, 1, 0] = np.nan
    elif kind == "diagonal":
        costs[0, 0, 0] = 1.0
    elif kind == "undersized":
        costs = costs[:, :3, :3]
    elif kind == "negative_static":
        costs = -SM5
    return costs


@pytest.mark.parametrize("kind", ["negative", "nan", "diagonal", "undersized", "negative_static"])
def test_low_level_explicit_dhd_refuses_invalid_costs(kind: str) -> None:
    with pytest.raises(ValueError, match="DHD"):
        pairwise(CODES, method="DHD", subm=_invalid(kind))


def test_r_refuses_the_same_invalid_costs() -> None:
    assert R_ERRORS == {
        "negative": "[!] Negative substitution costs are not allowed",
        "nan": "missing value where TRUE/FALSE needed",
        "diagonal": "[!] Non-zero diagonal elements in substitution cost matrix (sm)!",
        "undersized": "[!] size of substitution cost matrix must be 4x4",
    }


@pytest.mark.parametrize("norm", ["maxlength", "auto", "bogus"])
def test_legacy_dhd_refuses_norm(norm: str) -> None:
    with pytest.raises(ValueError, match="legacy DHD"):
        pairwise(CODES, method="DHD", variant="mean", norm=norm)
    # seqdist checks the norm spelling first, as R does (seqdist.R:127; TX-REP-20261019-010),
    # so an unknown spelling is refused as such before the legacy rule.
    expected = "norm must be one of" if norm == "bogus" else "legacy DHD"
    with pytest.raises(ValueError, match=expected):
        tx.seqdist(SEQS, method="DHD", variant="mean", norm=norm)


def test_legacy_dhd_accepts_no_norm() -> None:
    plain = pairwise(CODES, method="DHD", variant="mean")
    np.testing.assert_array_equal(pairwise(CODES, method="DHD", variant="mean", norm=None), plain)
    np.testing.assert_array_equal(
        tx.seqdist(SEQS, method="DHD", variant="mean", norm="none"), plain
    )
