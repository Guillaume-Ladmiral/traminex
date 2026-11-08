"""``seqecmpgroup`` computes its chi-square statistics exactly as R does.

``seqecmpgroup.R`` tabulates ``xtabs(ww ~ group + presence)`` -- a
``tapply(w, by, sum)`` -- and runs ``stats::chisq.test`` on it. In R 4.6.1
``chisq.test`` takes ``sr <- rowSums(x)``, ``sc <- colSums(x)``,
``n <- sum(x)``, ``E <- outer(sr, sc)/n``, a single scalar
``YATES <- min(0.5, abs(x - E))`` for a 2 x 2 table (0 otherwise) and
``STATISTIC <- sum((abs(x - E) - YATES)^2/E)``; ``sum``, ``rowSums`` and
``colSums`` accumulate in long double. traminex corrected each cell
separately and summed with NumPy in double precision and another order, so
statistics that tie in exact arithmetic differed from R's in the last bits
(and lost R's ~1e-31 residue on zero statistics), and seqecmpgroup, which
orders subsequences by statistic, put tied subsequences in a different order
from R (2026-09-22 audit S5-16). The table, statistic, expected counts,
residuals and row order now match R bit for bit.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5b.R``
(TraMineR 2.2-14 on R 4.6.1).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.events._mining import _chisq_test, _group_presence_table, _r_sum

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5b_chisq_arithmetic.txt"


def _load() -> dict[str, list[list[str]]]:
    rows: dict[str, list[list[str]]] = defaultdict(list)
    for line in GOLDEN.read_text().splitlines():
        if not line.startswith("#"):
            kind, *fields = line.split("\t")
            rows[kind].append(fields)
    return rows


R = _load()


def _numbers(text: str) -> np.ndarray:
    return np.array(text.split(" "), dtype=float)


@pytest.mark.parametrize("fields", R["r_sum"], ids=[fields[0] for fields in R["r_sum"]])
def test_r_sum_adds_in_r_order_and_precision(fields):
    """``sum``/``rowSums``/``colSums`` on matrices whose total depends on both."""
    _case, nrow, values, total, row_sums, column_sums = fields
    matrix = _numbers(values).reshape((int(nrow), -1), order="F")
    assert _r_sum(matrix) == float(total)
    rows = [_r_sum(matrix[row, :]) for row in range(matrix.shape[0])]
    columns = [_r_sum(matrix[:, column]) for column in range(matrix.shape[1])]
    np.testing.assert_array_equal(rows, _numbers(row_sums))
    np.testing.assert_array_equal(columns, _numbers(column_sums))


@pytest.mark.parametrize("fields", R["table"], ids=[fields[0] for fields in R["table"]])
def test_table_and_chisq_test_match_r_bit_for_bit(fields):
    _case, groups, presence, weights, cells, statistic, expected, residuals = fields
    codes = np.array(groups.split(" "), dtype=int) - 1
    nlevels = int(codes.max()) + 1
    table = _group_presence_table(
        _numbers(weights), codes, np.array(presence.split(" "), dtype=float), nlevels
    )
    np.testing.assert_array_equal(table, _numbers(cells).reshape((nlevels, 2), order="F"))
    stat, _pval, resid, expect, _sr = _chisq_test(table)
    assert stat == float(statistic)
    np.testing.assert_array_equal(expect, _numbers(expected).reshape((nlevels, 2), order="F"))
    np.testing.assert_array_equal(resid, _numbers(residuals).reshape((nlevels, 2), order="F"))


def _eseq(weighted: bool):
    rows = R["corpus"]
    eseq = tx.seqecreate(
        id=[int(fields[0]) for fields in rows],
        timestamp=[float(fields[1]) for fields in rows],
        event=[fields[2] for fields in rows],
    )
    if weighted:
        eseq.weights = _numbers(R["weights"][0][0])
    return eseq


@pytest.mark.parametrize("weighted", [False, True], ids=["unweighted", "weighted"])
@pytest.mark.parametrize("grouping", ["two", "three"])
def test_seqecmpgroup_rows_and_order_match_r(grouping, weighted):
    eseq = _eseq(weighted)
    index = np.arange(1, len(eseq) + 1)
    group = (
        np.where(index % 2 == 0, "x", "y")
        if grouping == "two"
        else np.array(["x", "y", "z"])[index % 3]
    )
    result = tx.seqecmpgroup(tx.seqefsub(eseq, min_support=2, max_k=2), group=group)
    case = f"{grouping}_{'weighted' if weighted else 'unweighted'}"
    expected = [fields[1:] for fields in R["row"] if fields[0] == case]
    assert result.subseq == [label for _, label, _ in expected]
    values = np.vstack([_numbers(numbers) for _, _, numbers in expected])
    columns = list(result.data.columns)
    for position, column in enumerate(columns):
        observed = result.data[column].to_numpy(dtype=float)
        if column == "p.value":
            # R's pchisq and SciPy's chi2.sf agree to rounding, not bit for bit.
            np.testing.assert_allclose(observed, values[:, position], rtol=1e-12, atol=0)
        else:
            np.testing.assert_array_equal(observed, values[:, position], err_msg=column)
