"""seqpm matches POSIX character classes as R's grep does (TX-AUD-20261010-054).

R's seqpm greps the pattern in the seqconc strings with R's extended regular
expressions, where bracket classes such as ``[[:upper:]]`` work. Python's
``re`` has no such classes and read ``[[:upper:]]`` as a set of the characters
"[", ":", "u", ..., so it silently matched nothing. Inside a bracket a
backslash is an ordinary character for R but an escape for ``re``, so
``[\\]A]`` must not read as a bracket holding "]" and "A". Fixture:
``scripts/r_goldens/make_goldens_seqpm_posix.R`` (LC_CTYPE C.UTF-8).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqpm_posix.csv",
    dtype=str,
    keep_default_na=False,
)
DATA = {
    "plain": tx.seqdef(
        np.array([list("ABCCAB"), list("ABAACC"), list("CCCCCC")], dtype=object),
        alphabet=["A", "B", "C"],
    ),
    "accents": tx.seqdef(
        np.array([list("AéB"), list("BBé"), list("BAB")], dtype=object), alphabet=["A", "B", "é"]
    ),
    "brackets": tx.seqdef(
        np.array([list("A]1b"), list("\\Ac2"), list("BBBx")], dtype=object),
        alphabet=["1", "2", "A", "B", "\\", "]", "b", "c", "x"],
    ),
}


@pytest.mark.parametrize(
    "row", GOLDEN.to_dict("records"), ids=lambda r: f"{r['data']}:{r['pattern']}"
)
def test_seqpm_matches_as_r(row: dict) -> None:
    result = tx.seqpm(DATA[row["data"]], row["pattern"])
    assert int(result["MTab"]["nbocc"].iloc[0]) == int(row["nbocc"])
    expected = [int(index) for index in row["index"].split("|")] if row["index"] else []
    assert result["MIndex"].tolist() == expected
