"""with_missing keeps the missing level for every domain, as R does (TX-AUD-20261011-013).

R's position-wise table has every level of the stslist factors and removes void
and, without with.missing, the missing state (seqdomassoc.R:206-208): with
with.missing = TRUE the missing row and column stay for every domain, also one
without missing cells, and their zero cells enter the df adjustment. traminex
added the missing level only to a domain with missing cells (df 1 where R has
4). Fixture: ``scripts/r_goldens/make_goldens_domassoc_missing_level.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "domassoc_missing_level.csv",
    dtype=str,
    keep_default_na=False,
)
WITH_NA = tx.seqdef(
    np.array(
        [
            ["A", "A", "B"],
            ["B", None, "B"],
            [None, "B", "A"],
            ["A", "B", "A"],
            ["B", "B", "B"],
            ["A", "A", "B"],
        ],
        dtype=object,
    ),
    alphabet=["A", "B"],
)
WITHOUT_NA = tx.seqdef(
    np.array(
        [
            ["X", "Y", "X"],
            ["Y", "Y", "Y"],
            ["Y", "X", "Y"],
            ["X", "X", "Y"],
            ["Y", "Y", "X"],
            ["X", "X", "X"],
        ],
        dtype=object,
    ),
    alphabet=["X", "Y"],
)
ALSO_NA = tx.seqdef(
    np.array(
        [
            ["X", "Y", "X"],
            [None, "Y", "Y"],
            ["Y", "X", "Y"],
            ["X", None, "Y"],
            ["Y", "Y", "X"],
            ["X", "X", "X"],
        ],
        dtype=object,
    ),
    alphabet=["X", "Y"],
)
CASES = {
    "one_domain_without_missing": [WITH_NA, WITHOUT_NA],
    "both_domains_with_missing": [WITH_NA, ALSO_NA],
}


def _r(case: str, struct_zero: bool, section: str) -> dict[str, float]:
    rows = GOLDEN[
        (GOLDEN["case"] == case)
        & (GOLDEN["struct_zero"] == str(struct_zero).upper())
        & (GOLDEN["section"] == section)
    ]
    return {name: float(value) for name, value in zip(rows["name"], rows["value"])}


@pytest.mark.parametrize("struct_zero", [True, False])
@pytest.mark.parametrize("case", sorted(CASES))
def test_missing_level_is_kept_for_every_domain(case: str, struct_zero: bool) -> None:
    result = tx.seqdomassoc(
        CASES[case], with_missing=True, cross_table=True, struct_zero=struct_zero
    )
    table = next(iter(result.attrs["cross_tables"].values()))
    ours = {
        f"{row}|{column}": float(table.loc[row, column])
        for row in table.index
        for column in table.columns
    }
    assert ours == _r(case, struct_zero, "table")
    expected = _r(case, struct_zero, "stat")
    np.testing.assert_allclose(
        result.iloc[0][list(expected)].to_numpy(dtype=float),
        list(expected.values()),
        rtol=0,
        atol=1e-10,
    )
