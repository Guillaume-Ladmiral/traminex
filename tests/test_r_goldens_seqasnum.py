"""R-golden parity for seqasnum / seqnum (T11) vs TraMineR 2.2-14.

seqasnum is the numeric recoding of the sequences (0-based state codes; NaN for
void always and for nr unless with_missing). seqnum returns the same sequences
over a numeric alphabet (codes unchanged, so seqasnum round-trips). Both are
re-targets from the previous distinct-state-count implementation. Goldens:
scripts/r_goldens/make_goldens_seqasnum.R (NaN -> -1 sentinel).
"""

from __future__ import annotations

import numpy as np
from corestats_fixtures import fixture_a, fixture_b, fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.stats import seqasnum, seqnum


def _with_sentinel(out: np.ndarray) -> np.ndarray:
    out = np.asarray(out, dtype=float).copy()
    out[np.isnan(out)] = -1
    return out


def _check(x, tag: str, with_missing: bool):
    got = _with_sentinel(seqasnum(x, with_missing=with_missing))
    np.testing.assert_allclose(got, load_golden(f"seqasnum_{tag}"), rtol=0, atol=R_GOLDEN_ATOL)


def test_seqasnum_a():
    _check(fixture_a(), "A", False)


def test_seqasnum_b_with_and_without_missing():
    _check(fixture_b(), "B_wm0", False)
    _check(fixture_b(), "B_wm1", True)


def test_seqasnum_c_ragged_with_and_without_missing():
    _check(fixture_c(), "C_wm0", False)
    _check(fixture_c(), "C_wm1", True)


def test_seqnum_numeric_alphabet_preserves_codes():
    """seqnum recodes the alphabet to 0..K-1 without changing the codes (seqasnum)."""
    x = fixture_c()
    num = seqnum(x)
    # real-state labels become their numeric codes; void symbol preserved.
    for i, lab in enumerate(num.alphabet):
        if i == x.void_code:
            assert lab == x.void
        elif i == x.nr_code:
            assert lab == x.nr  # with_missing=False keeps the nr symbol
        else:
            assert lab == str(i)
    # codes (and therefore seqasnum) are unchanged by the relabel.
    np.testing.assert_array_equal(_with_sentinel(seqasnum(num)), _with_sentinel(seqasnum(x)))
