"""R-golden tests for seqemlt (EMLT) -- B-2 / E-2.

Anchors traminex's EMLT to TraMineRextras 0.6.9 ``seqemlt`` (see
scripts/r_goldens/make_goldens_seqemlt.R). The sign-stable quantities are
goldened exactly (situation correlation matrix, profile squared-distance matrix,
and princomp standard deviations); the sequence coordinates -- which come from
``princomp`` on the square, rank-deficient Benzecri matrix -- are compared up to
a per-component sign on the non-degenerate components (plan decision 0.2.3).
"""

from __future__ import annotations

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import seqdef, seqemlt

SEQMAT = [
    list("AAAAA"),
    list("BBBBB"),
    list("CCCCC"),
    list("ABCAB"),
    list("CABCA"),
    list("BCABC"),
    list("ACBAC"),
    list("BACBA"),
]


def _result():
    sts = seqdef(SEQMAT, alphabet=["A", "B", "C"])
    return seqemlt(sts, a=1, b=1, weighted=False)


def test_seqemlt_correlation_matches_r():
    r = _result()
    golden = load_golden("seqemlt_corr")
    np.testing.assert_allclose(r.correlation, golden, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqemlt_profile_distance_matches_r():
    r = _result()
    golden = load_golden("seqemlt_distc")
    np.testing.assert_allclose(r.distance_matrix, golden, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqemlt_sdev_matches_r():
    r = _result()
    golden = load_golden("seqemlt_sdev").ravel()
    # princomp standard deviations are eigenvalues of the situation correlation
    # matrix (itself R-faithful to atol=1e-10); the eigendecomposition limits the
    # agreement to ~1e-8, so this is a 6-decimal check on the spectrum.
    np.testing.assert_allclose(r.sdev, golden, rtol=0, atol=1e-6)


def test_seqemlt_coordinates_match_r():
    # princomp coordinates are individually sign-/rotation-arbitrary (sign per
    # component, plus rotation within near-degenerate components and the final
    # rank-deficient zero component). The sequence Gram matrix coord @ coord.T is
    # invariant to all of that and captures the clustering-relevant content; it
    # matches R essentially exactly.
    r = _result()
    golden = load_golden("seqemlt_coords")
    actual = np.asarray(r.coordinates)
    np.testing.assert_allclose(actual @ actual.T, golden @ golden.T, rtol=0, atol=R_GOLDEN_ATOL)
