"""``seqecmpgroup`` groups sequences the way R's ``factor(group)`` does.

``seqecmpgroup.R:62`` converts the grouping with ``factor(group)``: numbers
sort as numbers and are named as R prints them (``2 10``, ``1.5 1e+05``,
``FALSE TRUE``), a factor keeps its declared level order and drops unused
levels, and a missing group value is left out of every group-by-presence table
while still counting in ``length(group)``. traminex turned every value into
text and sorted the text, so numeric groups came out ``10 2`` with their
``Freq.*``/``Resid.*`` columns swapped, a pandas categorical lost its order,
``b B a`` sorted as ``B a b`` instead of R's ``a b B``, and a missing value
became a group called ``"None"`` -- a different test altogether (2026-09-22
audit S5-02). Python's missing values (``None``, ``NaN``, ``pd.NA``) play the
part of R's ``NA``.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5.R``
(TraMineR 2.2-14). Rows and their order are R's exactly; only the p-values,
from SciPy's ``chi2.sf`` rather than R's ``pchisq``, are compared to rounding
(the chi-square arithmetic itself is R's since audit S5-16).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5_cmpgroup_factor.txt"


def _load() -> dict[str, list[list[str]]]:
    rows: dict[str, list[list[str]]] = defaultdict(list)
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        kind, *fields = line.split("\t")
        rows[kind].append(fields)
    return rows


R = _load()


def _subsequences() -> tx.SubSeqElist:
    rows = R["corpus"]
    eseq = tx.seqecreate(
        id=[int(fields[0]) for fields in rows],
        timestamp=[float(fields[1]) for fields in rows],
        event=[fields[2] for fields in rows],
    )
    return tx.seqefsub(eseq, min_support=3, max_k=2)


def _groups(n: int) -> dict[str, object]:
    index = np.arange(1, n + 1)
    missing: list[object] = ["x" if i % 3 == 0 else "y" for i in index]
    for position in (2, 5, 11):
        missing[position - 1] = None
    return {
        "numeric": np.where(index % 3 == 0, 10, 2),
        "character": [("b", "B", "a")[(1 + i % 3) - 1] for i in index],
        "factor": pd.Categorical(
            ["x" if i % 2 == 0 else "y" for i in index], categories=["z", "y", "x"]
        ),
        "missing": missing,
        "logical": index % 2 == 0,
        "fractional": np.where(index % 2 == 0, 1.5, 1e5),
    }


FSUB = _subsequences()
GROUPS = _groups(len(FSUB.eseq))


@pytest.mark.parametrize("tag", list(GROUPS))
def test_group_levels_and_columns_follow_r_factor(tag):
    result = tx.seqecmpgroup(FSUB, group=GROUPS[tag])
    (labels,) = [fields[1] for fields in R["labels"] if fields[0] == tag]
    (columns,) = [fields[1] for fields in R["columns"] if fields[0] == tag]
    assert result.labels == labels.split(" ")
    assert list(result.data.columns) == columns.split(" ")


@pytest.mark.parametrize("tag", list(GROUPS))
def test_group_statistics_match_r(tag):
    result = tx.seqecmpgroup(FSUB, group=GROUPS[tag])
    expected = [fields[1:] for fields in R["row"] if fields[0] == tag]
    assert result.subseq == [label for _, label, _ in expected]
    r_rows = np.vstack([np.array(values.split(" "), dtype=float) for _, _, values in expected])
    for position, column in enumerate(result.data.columns):
        observed = result.data[column].to_numpy(dtype=float)
        if column == "p.value":
            np.testing.assert_allclose(observed, r_rows[:, position], rtol=1e-12, atol=0)
        else:
            np.testing.assert_array_equal(observed, r_rows[:, position], err_msg=column)


def test_no_subsequence_passing_pvalue_limit_gives_an_empty_result():
    """A recorded divergence (2026-09-22 audit S5-07), not R parity.

    ``seqecmpgroup.R:88`` keeps ``1:sum(p.value <= pvalue.limit)`` rows; when
    none passes, ``1:0`` is ``c(1, 0)`` and R still returns the top
    subsequence, whose p-value fails the limit. traminex keeps only the
    subsequences that pass, so it returns none.
    """
    ((limit, rows, labels, pvalues),) = R["pvalue_limit"]
    assert int(rows) == 1 and float(pvalues) > float(limit)  # R's 1:0 artefact
    result = tx.seqecmpgroup(FSUB, group=GROUPS["numeric"], pvalue_limit=float(limit))
    top = tx.seqecmpgroup(FSUB, group=GROUPS["numeric"])
    assert top.subseq[0] == labels
    assert len(result.subseq) == 0 and len(result.data) == 0
