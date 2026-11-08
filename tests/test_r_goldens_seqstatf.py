"""R-golden parity for seqstatf (T10) vs TraMineR 2.2-14.

seqstatf is the OVERALL state frequency (weighted total time in each state),
re-targeted from the previous first-position implementation. Index is statl
(real alphabet, + nr when include_missing). Goldens:
scripts/r_goldens/make_goldens_seqstatf.R.
"""

from __future__ import annotations

import numpy as np
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import seqstatf


def _check(x, tag: str, include_missing: bool):
    f = seqstatf(x, weighted=True, include_missing=include_missing)
    got = f.to_numpy()  # [Freq, Percent], rows = statl
    np.testing.assert_allclose(got, load_golden(f"seqstatf_{tag}"), rtol=0, atol=R_GOLDEN_ATOL)


def test_seqstatf_a_unweighted():
    _check(fixture_a(), "A", False)


def test_seqstatf_a_weighted():
    _check(fixture_a(weighted=True), "Aw", False)


def test_seqstatf_b_with_and_without_missing():
    _check(fixture_b(), "B_wm0", False)
    _check(fixture_b(), "B_wm1", True)


def test_seqstatf_c_ragged_with_and_without_missing():
    _check(fixture_c(), "C_wm0", False)
    _check(fixture_c(), "C_wm1", True)
