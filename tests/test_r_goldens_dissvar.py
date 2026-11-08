"""R-golden tests for the weighted discrepancy ``dissvar`` (Phase 2.1 / A5).

Anchors traminex's ``dissvar`` to TraMineR 2.2-14 ``dissvar(diss, weights=,
squared=)`` on the shared Phase 2 fixture (see
``scripts/r_goldens/make_goldens_phase2_inputs.R`` /
``make_goldens_dissvar.R``).
"""

from __future__ import annotations

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.diss import dissvar, weighted_inertia, weighted_pairwise_sum


def _inputs():
    D = load_golden("phase2_input_D")
    wreal = load_golden("phase2_input_weights_real").ravel()
    wint = load_golden("phase2_input_weights_int").ravel()
    return D, wreal, wint


def test_dissvar_matches_r_golden():
    D, wreal, wint = _inputs()
    golden = load_golden("dissvar_phase2").ravel()
    actual = np.array(
        [
            dissvar(D),
            dissvar(D, weights=wreal),
            dissvar(D, weights=wreal, squared=True),
            dissvar(D, squared=True),
            dissvar(D, weights=wint),
        ]
    )
    np.testing.assert_allclose(actual, golden, rtol=0, atol=R_GOLDEN_ATOL)


def test_dissvar_unweighted_equals_unit_weights():
    # R's unweighted branch (sum(diss)/(2n^2)) must agree with the weighted
    # path under unit weights — same identity as in the C routine.
    D, _, _ = _inputs()
    n = D.shape[0]
    np.testing.assert_allclose(dissvar(D), dissvar(D, weights=np.ones(n)), rtol=0, atol=1e-12)


def test_weighted_inertia_consistency():
    # weighted_inertia is the SS building block: dissvar == inertia / totw,
    # and the subset form must equal the masked-matrix form.
    D, wreal, _ = _inputs()
    totw = wreal.sum()
    np.testing.assert_allclose(
        dissvar(D, weights=wreal),
        weighted_inertia(D, None, wreal) / totw,
        rtol=0,
        atol=1e-12,
    )
    idx = np.array([0, 3, 7, 8, 12])
    sub = D[np.ix_(idx, idx)]
    wsub = wreal[idx]
    np.testing.assert_allclose(
        weighted_pairwise_sum(D, idx, wreal),
        weighted_pairwise_sum(sub, None, wsub),
        rtol=0,
        atol=1e-12,
    )
    iu, ju = np.triu_indices(idx.size, k=1)
    manual = float(np.sum(wsub[iu] * wsub[ju] * sub[iu, ju]))
    np.testing.assert_allclose(weighted_pairwise_sum(D, idx, wreal), manual, rtol=0, atol=1e-12)
