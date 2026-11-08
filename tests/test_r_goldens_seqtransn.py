"""R-golden parity for seqtransn (T4) vs TraMineR 2.2-14.

seqtransn is the per-sequence transition count (length(seqdss) - 1), re-pointed
from the previous dataset-level K x K matrix. Pins counts, norm=True (divide by
the non-void length - 1) and the pweight path. Goldens:
scripts/r_goldens/make_goldens_seqtransn.R.
"""

from __future__ import annotations

import numpy as np
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.stats import seqtransn


def _check(x, tag: str, with_missing: bool):
    golden = load_golden(f"seqtransn_{tag}")
    counts = seqtransn(x, with_missing=with_missing, norm=False)
    normed = seqtransn(x, with_missing=with_missing, norm=True)
    np.testing.assert_allclose(counts, golden[:, 0], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(normed, golden[:, 1], rtol=0, atol=R_GOLDEN_ATOL)


def test_seqtransn_a():
    _check(fixture_a(), "A", False)


def test_seqtransn_b_with_and_without_missing():
    _check(fixture_b(), "B_wm0", False)
    _check(fixture_b(), "B_wm1", True)


def test_seqtransn_c_ragged_with_and_without_missing():
    _check(fixture_c(), "C_wm0", False)
    _check(fixture_c(), "C_wm1", True)


def test_seqtransn_pweight():
    np.testing.assert_allclose(
        seqtransn(fixture_a(), pweight=True),
        load_golden("seqtransn_pweight_A")[:, 0],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    np.testing.assert_allclose(
        seqtransn(fixture_c(), pweight=True),
        load_golden("seqtransn_pweight_C")[:, 0],
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
