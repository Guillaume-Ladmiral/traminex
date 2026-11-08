"""OMslen's link is checked before computing (TX-AUD-20261011-051).

R refuses a link outside none/mean/gmean for every method (seqdist.R:185-186)
and "none" for OMslen ("'link' must be specified", :229-233), whatever the
number of sequences. traminex checked the link only inside the pair loop, so
with one sequence or identical rows (no pair computed) an invalid link passed,
at the low level and, for "none", through seqdist. Fixture
``scripts/r_goldens/make_goldens_omslen_link.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.distances import python_kernels as pk
from traminex.sequences import STSList

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "omslen_link.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
S3 = tx.seqdef(
    np.array([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"]], dtype=object), alphabet=["A", "B"]
)


def _subset(rows: list[int]) -> STSList:
    return STSList(
        data=S3.data[rows].copy(),
        alphabet=list(S3.alphabet),
        nr_code=S3.nr_code,
        void_code=S3.void_code,
    )


DATA = {"one": _subset([0]), "identical": _subset([0, 0]), "distinct": S3}
CODES = {name: np.asarray(seqs.data, dtype=np.int64) for name, seqs in DATA.items()}


@pytest.mark.parametrize("data", sorted(DATA))
@pytest.mark.parametrize("link", ["none", "bogus", "mean", "gmean"])
def test_seqdist_matches_r(link: str, data: str) -> None:
    row = GOLDEN.loc[f"{link}_{data}"]
    if row["ok"] == "FALSE":
        with pytest.raises(ValueError, match="link"):
            tx.seqdist(DATA[data], method="OMslen", subm="CONSTANT", link=link)
        return
    expected = np.array([float(v) for v in row["value"].split()])
    actual = np.asarray(tx.seqdist(DATA[data], method="OMslen", subm="CONSTANT", link=link))
    np.testing.assert_allclose(actual.ravel(order="F"), expected, atol=1e-12, rtol=0)


@pytest.mark.parametrize("data", sorted(CODES))
@pytest.mark.parametrize("link", ["none", "bogus"])
@pytest.mark.parametrize(
    "pairwise", [pk.pairwise, tx.distances.pairwise], ids=["python_kernels", "dispatcher"]
)
def test_the_low_level_layer_checks_the_link_first(pairwise, link: str, data: str) -> None:
    with pytest.raises(ValueError, match="link"):
        pairwise(CODES[data], method="OMSLEN", link=link, subm=2.0, indel=1.0)
