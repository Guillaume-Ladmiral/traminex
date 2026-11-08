"""seqstatl decides numeric states with R's as.numeric (TX-AUD-20261010-050).

R's seqstatl returns the states as numbers, in numeric order, when every
distinct value converts with as.numeric. traminex used Python's float(), which
reads "1_1" as 11 (digit-group underscores) and refuses "0x10" (hex), so it
returned states that are not in the data or the wrong kind. Fixture:
``scripts/r_goldens/make_goldens_seqstatl_r_numbers.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqstatl_r_numbers.csv",
    dtype=str,
    keep_default_na=False,
)
# The matrices of the generator (R fills them by column).
CASES = {
    "underscores": np.array([["1_1", "2_1"], ["1_2", "1_1"]], dtype=object),
    "hex": np.array([["0x10", "2"], ["2", "0x10"]], dtype=object),
    "spellings": np.array([[" 2", "0x1A", "3"], ["1e2", "Inf", " 2"]], dtype=object),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_seqstatl_reads_numbers_as_r_does(case: str) -> None:
    rows = GOLDEN[GOLDEN["case"] == case].sort_values("position", key=lambda c: c.astype(int))
    result = tx.seqstatl(CASES[case])
    if rows["kind"].iloc[0] == "numeric":
        assert result == [float(value) for value in rows["value"]]
        assert all(isinstance(value, float) for value in result)
    else:
        assert result == rows["value"].tolist()
