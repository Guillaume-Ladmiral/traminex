"""NMS, NMSMST and SVRspell add their weighted counts in R's order (TX-REP-20261022-001).

R adds ``kweights[k] * count[k]`` one subsequence length after another and skips zero weights
(``NMSdistance.cpp``, ``SUBSEQdistance::distance``). traminex summed with ``np.sum``, whose pairwise summation
rounds differently, so distances differed from R's in the last bit. R's values are the NMSMST, SVRspell and NMS
rows of ``scripts/r_goldens/make_goldens_spell_powers.R`` (hexadecimal).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_r_goldens_spell_powers import X, _r_values  # noqa: E402

CASES = {
    "NMSMST_tpow_1.3": dict(method="NMSMST", tpow=1.3),
    "SVRspell_tpow_1.3": dict(method="SVRspell", tpow=1.3),
    "NMS": dict(method="NMS"),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_subsequence_distances_are_r_s_bit_for_bit(case: str) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        actual = np.asarray(tx.seqdist(X, **CASES[case]), dtype=float).ravel(order="F")
    np.testing.assert_array_equal(actual, _r_values("dist", case))
