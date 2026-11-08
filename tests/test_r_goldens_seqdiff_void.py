"""R-golden tests for seqdiff on ragged / missing data (RP3/RP4 audit F2).

The original seqdiff golden used a no-missing fixture, so the with_missing path
was unexercised. With with_missing=True, R keeps the missing state nr as a state
and DELETES the structural void of ragged sequences before the within-window LCS.
traminex now matches (the void code flows through to the distance layer). See
scripts/r_goldens/make_goldens_seqdiff_void.R.
"""

from __future__ import annotations

import numpy as np
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import seqdef, seqdiff

_ALPH = ["1", "2", "3", "4"]

# Matched fixture: interior "*" -> nr, trailing "*" (right="DEL" in R) -> void "%".
_M = [
    ["1", "2", "3", "4", "1", "2"],
    ["2", "1", "*", "3", "4", "1"],
    ["3", "3", "1", "2", "2", "3"],
    ["4", "1", "2", "*", "2", "1"],
    ["1", "1", "2", "3", "4", "2"],
    ["2", "3", "4", "1", "1", "3"],
    ["3", "2", "1", "4", "3", "%"],
    ["4", "4", "1", "1", "1", "%"],
    ["1", "3", "2", "4", "3", "%"],
    ["2", "2", "3", "3", "4", "%"],
    ["3", "1", "2", "2", "1", "%"],
    ["4", "2", "1", "3", "2", "%"],
]


def _fixture():
    w = 0.5 + 0.1 * np.arange(1, len(_M) + 1)
    return seqdef(np.array(_M, dtype=object), alphabet=_ALPH, nr="*", void="%", weights=w)


@pytest.mark.parametrize("tag,weighted", [("uw", False), ("w", True)])
def test_seqdiff_with_missing_ragged_matches_live_r(tag, weighted):
    x = _fixture()
    g = np.array([1, 2] * (len(_M) // 2))
    out = seqdiff(x, g, method="LCS", norm=None, with_missing=True, weighted=weighted)
    np.testing.assert_allclose(
        out["stat"].to_numpy(dtype=float),
        load_golden(f"seqdiff_void_{tag}_stat"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        out["discrepancy"].to_numpy(dtype=float),
        load_golden(f"seqdiff_void_{tag}_disc"),
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
