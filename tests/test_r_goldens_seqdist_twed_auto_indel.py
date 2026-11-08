"""seqdist(method="TWED") with the automatic indel (TX-AUD-20261010-011).

R's default ``indel = "auto"`` resolves for TWED to ``2 * max(sm) + nu + h`` over
the whole substitution matrix (seqdist.R:359-364, 433-439); traminex refused the
string. Fixture: ``scripts/r_goldens/make_goldens_seqdist_twed_auto_indel.R``.
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

GOLDEN = Path(__file__).resolve().parent / "r_goldens" / "seqdist_twed_auto_indel.csv"
SM5 = np.abs(np.subtract.outer(np.arange(4), np.arange(4))).astype(float)


def _golden() -> dict[str, np.ndarray]:
    cells: dict[str, dict[tuple[int, int], float]] = defaultdict(dict)
    with GOLDEN.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cells[row["case"]][(int(row["row"]), int(row["col"]))] = float(row["value"])
    out = {}
    for case, values in cells.items():
        size = max(i for i, _ in values)
        matrix = np.empty((size, size))
        for (i, j), value in values.items():
            matrix[i - 1, j - 1] = value
        out[case] = matrix
    return out


R = _golden()


def _seqs(strings: list[str], alphabet: list[str], **kwargs) -> tx.STSList:
    width = max(len(s) for s in strings)
    rows = [list(s) + [None] * (width - len(s)) for s in strings]
    return tx.seqdef(rows, alphabet=alphabet, **kwargs)


S4 = _seqs(["AABBC", "ABBCC", "CCBAA", "BBCCA"], ["A", "B", "C"])
S5 = _seqs(["ABBCDA", "BBCCDD", "DACBBA", "AADDCC", "CBADCB"], ["A", "B", "C", "D"])
S5R = _seqs(["ABBCDA", "BBCC", "DACBBA", "AAD"], ["A", "B", "C", "D"], right="DEL")

CASES = {
    "s4_constant_auto": (S4, dict(subm="CONSTANT", nu=0.5, h=0.5)),
    "s5_matrix_auto": (S5, dict(subm=SM5, nu=0.3, h=0.2)),
    "s5_trate_auto": (S5, dict(subm="TRATE", nu=0.3, h=0.2)),
    "s5_matrix_auto_yujianbo": (S5, dict(subm=SM5, nu=0.3, h=0.2, norm="YujianBo")),
    "s5r_matrix_auto": (S5R, dict(subm=SM5, nu=0.3, h=0.2)),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_twed_auto_indel_matches_r(case: str) -> None:
    seqs, options = CASES[case]
    actual = tx.seqdist(seqs, method="TWED", indel="auto", **options)
    np.testing.assert_allclose(actual, R[case], rtol=0, atol=R_GOLDEN_ATOL)


def test_twed_omitted_indel_equals_auto() -> None:
    omitted = tx.seqdist(S5, method="TWED", subm=SM5, nu=0.3, h=0.2)
    np.testing.assert_allclose(omitted, R["s5_matrix_auto"], rtol=0, atol=R_GOLDEN_ATOL)


def test_low_level_pairwise_still_refuses_auto() -> None:
    codes = np.asarray(S4.data, dtype=np.int32)
    with pytest.raises(ValueError, match="resolved by seqdist"):
        pairwise(codes, method="TWED", indel="auto", nu=0.5, h=0.5)
