"""R-golden parity for seqmodst (T7) vs TraMineR 2.2-14.

Modal-state sequence: the per-position modal state codes, the modal-state
frequency at each position, and the number of sequences exactly matching the
modal sequence. Goldens: scripts/r_goldens/make_goldens_seqmodst.R.

traminex's modal codes are indices over the STSList alphabet (real states, nr,
void); R's seqasnum(m, with.missing=TRUE) codes the real alphabet then nr, which
coincides because nr is appended right after the real states.
"""

from __future__ import annotations

import numpy as np
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import seqmodst


def _check(x, tag: str, include_missing: bool):
    m = seqmodst(x, weighted=True, include_missing=include_missing)
    cf = load_golden(f"seqmodst_{tag}_cf")
    np.testing.assert_allclose(
        np.asarray(m.data[0], dtype=float), cf[0], rtol=0, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(
        np.asarray(m.frequencies_, dtype=float), cf[1], rtol=0, atol=R_GOLDEN_ATOL
    )
    occ = load_golden(f"seqmodst_{tag}_occ")
    assert float(m.occurrences_) == occ[0, 0]


def test_seqmodst_a_unweighted():
    _check(fixture_a(), "A", False)


def test_seqmodst_a_weighted():
    _check(fixture_a(weighted=True), "Aw", False)


def test_seqmodst_b_with_missing():
    _check(fixture_b(), "B_wm1", True)


def test_seqmodst_c_ragged_with_missing():
    _check(fixture_c(), "C_wm1", True)
