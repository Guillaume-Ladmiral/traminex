"""R-golden parity for seqmeant (T6) vs TraMineR 2.2-14.

Mean time spent in each state (+ Var/Stdev/SE), and the prop=True path
(per-sequence proportions averaged). R's rows are statl (real alphabet, + nr
when with.missing), void never a row. Goldens:
scripts/r_goldens/make_goldens_seqmeant.R.
"""

from __future__ import annotations

import numpy as np
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.stats import seqmeant


def _state_rows(x, include_missing: bool):
    return [
        i
        for i in range(len(x.alphabet))
        if i != x.void_code and (include_missing or i != x.nr_code)
    ]


def _check(x, tag: str, include_missing: bool):
    rows = _state_rows(x, include_missing)
    ms = seqmeant(x, weighted=True, include_missing=include_missing, prop=False, serr=True)
    got = ms.to_numpy()[rows][:, :4]  # Mean, Var, Stdev, SE
    np.testing.assert_allclose(got, load_golden(f"seqmeant_{tag}_ms"), rtol=0, atol=R_GOLDEN_ATOL)

    pr = seqmeant(x, weighted=True, include_missing=include_missing, prop=True, serr=False)
    got_p = pr.to_numpy()[rows][:, :1]
    np.testing.assert_allclose(
        got_p, load_golden(f"seqmeant_{tag}_prop"), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqmeant_a_unweighted():
    _check(fixture_a(), "A", False)


def test_seqmeant_a_weighted():
    _check(fixture_a(weighted=True), "Aw", False)


def test_seqmeant_b_with_and_without_missing():
    _check(fixture_b(), "B_wm0", False)
    _check(fixture_b(), "B_wm1", True)


def test_seqmeant_c_ragged_with_and_without_missing():
    _check(fixture_c(), "C_wm0", False)
    _check(fixture_c(), "C_wm1", True)
