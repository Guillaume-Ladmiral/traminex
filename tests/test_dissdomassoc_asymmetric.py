"""``dissdomassoc`` reads an asymmetric matrix from the triangle R reads.

Audit 2026-09-22, finding S4-09 (the ``dissdomassoc`` part). TraMineR 2.2-14
vectorises every domain matrix and ``jointdiss`` with ``as.dist``
(``SOURCES/TraMineR/R/dissdomassoc.R:55``, ``:63``), which takes the lower
triangle column by column. traminex took the upper triangle row by row: the
same numbers in the same order for a symmetric matrix, but different numbers
for an asymmetric one (TraMineR accepts asymmetric substitution costs with a
warning). Expected values come from live R through
``scripts/r_goldens/make_goldens_audit20260922_s4.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL

from traminex.diss import dissdomassoc

FIXTURE = Path(__file__).with_name("r_goldens") / "dissdomassoc_asymmetric.txt"


def _rows() -> dict[str, np.ndarray]:
    rows: dict[str, np.ndarray] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        key, value = line.split("\t")
        rows[key] = np.array([float(v) for v in value.split(" ")])
    return rows


ROWS = _rows()


def _square(key: str) -> np.ndarray:
    values = ROWS[key]
    return values.reshape(int(round(np.sqrt(values.size))), -1)


def test_asymmetric_domains_match_r():
    A, S, J = _square("A"), _square("S"), _square("J")
    assert not np.array_equal(A, A.T) and not np.array_equal(J, J.T)

    result = dissdomassoc(
        [A, S], jointdiss=J, what=("pearson", "spearman", "R2"), weights=ROWS["weights"]
    )

    for key in ("Pearson", "Spearman"):
        expected = _square(key).T  # R wrote the matrix column-major
        off_diagonal = ~np.eye(expected.shape[0], dtype=bool)
        np.testing.assert_allclose(
            result[key][off_diagonal], expected[off_diagonal], rtol=1e-12, atol=R_GOLDEN_ATOL
        )
    for key in ("Pearson.Rsquare", "Spearman.Rsquare"):
        np.testing.assert_allclose(result[key], ROWS[key], rtol=1e-12, atol=R_GOLDEN_ATOL)


def test_symmetric_input_reads_the_same_values_either_way():
    S = _square("S")
    lower = dissdomassoc([S, S * 2.0 + 1.0], what=("pearson",))
    upper = dissdomassoc([S.T, (S * 2.0 + 1.0).T], what=("pearson",))
    np.testing.assert_array_equal(lower["Pearson"], upper["Pearson"])
