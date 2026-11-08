"""``seqauto`` builds R's contingency table and applies R's ``chisq.test`` rules.

TraMineRextras ``seqauto.R`` line 12 cross-tabulates the current state against
the lagged one with ``table(srs[, nc], srs[, nc - i])``. ``table`` keeps every
level that either column takes, including a state seen only beside an ``NA`` in
the other column (for instance a state that occurs only at the first
position); such a level contributes an all-zero row or column, and
``chisq.test`` (R ``stats``) then divides 0 by an expected count of 0, so the
statistic, the p-value and Cramer's V are NaN. traminex cross-tabulated only
the complete pairs, dropped those levels and returned finite values (2026-09-22
audit S3-13).

The same port also missed ``chisq.test``'s treatment of a table with a single
row or column (found while verifying S3-13): R turns it into a goodness-of-fit
test against equal probabilities, so the p-value is that test's and Cramer's V
divides by ``min(dim) - 1 = 0`` (``Inf``), and it refuses a 1 x 1 table with
"'x' must at least have 2 elements". traminex returned NaN and a p-value of 1.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3b.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s3b.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, _warnings = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    return np.array([float(v) for v in payload.split(" ")])


CASES = {
    "level_seen_only_beside_na": (
        [["X", "A", "A", "B"], ["X", "B", "B", "A"], ["A", "A", "B", "B"], ["B", "A", "A", "A"]],
        ["A", "B", "X"],
        1,
    ),
    "zero_row_after_first_position": ([["B", "A", "A"], ["A", "A", "A"]], ["A", "B"], 1),
    "level_only_in_current_column": (
        [["A", "B", "A", "B"], ["B", "A", "B", "A"], ["A", "A", "B", "B"], ["B", "B", "A", "C"]],
        ["A", "B", "C"],
        1,
    ),
    "one_column": ([["A", "A", "B"], ["A", "A", "A"]], ["A", "B"], 1),
    "control_order2": (
        [["A", "B", "A", "B"], ["B", "A", "B", "A"], ["A", "A", "B", "B"], ["B", "B", "A", "A"]],
        ["A", "B"],
        2,
    ),
    "control_lag2_uneven": (
        [["B", "A", "A", "B"], ["A", "A", "A", "A"], ["A", "B", "A", "A"]],
        ["A", "B"],
        2,
    ),
    "unobserved_alphabet_state": (
        [["A", "B", "A", "B"], ["B", "A", "B", "A"], ["A", "A", "B", "B"], ["B", "B", "A", "A"]],
        ["A", "B", "Z"],
        1,
    ),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_seqauto_matches_r_table_and_chisq_test(name):
    rows, alphabet, order = CASES[name]
    x = tx.seqdef(np.array(rows, dtype=object), alphabet=alphabet)
    result = tx.seqauto(x, order=order)
    assert result.shape == (order, 2)
    np.testing.assert_allclose(
        result.ravel(), _values(f"s313_{name}"), rtol=0, atol=1e-10, equal_nan=True
    )


def test_seqauto_skips_missing_cells_like_r():
    rows = [["A", None, "A", "B"], ["B", "A", "B", "A"], ["A", "A", "B", "B"], ["B", "B", "A", "A"]]
    x = tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B"], right=None)
    np.testing.assert_allclose(
        tx.seqauto(x, order=1).ravel(), _values("s313_missing_cells"), rtol=0, atol=1e-10
    )


def test_seqauto_refuses_a_one_by_one_table_as_r_does():
    kind, message = R["s313_one_by_one"]
    assert kind == "error"
    x = tx.seqdef(np.array([["A", "A", "A"], ["A", "A", "A"]], dtype=object), alphabet=["A", "B"])
    with pytest.raises(ValueError, match=message):
        tx.seqauto(x, order=1)
