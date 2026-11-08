"""R-golden parity for seqrep_grp (RP6a) vs TraMineRextras 0.6.9 seqrep.grp.

seqrep.grp runs seqrep() per group with a shared dmax = max(diss) and returns the
per-group quality Statistics. Reuses the make_goldens_dissrep.R seqrep fixture
(10 seqs, LCS norm=none) split into two groups of five. Goldens:
scripts/r_goldens/make_goldens_seqrep_grp.R.
"""

from __future__ import annotations

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex import seqdef, seqdist
from traminex.cluster import seqrep_grp

_ALPH = ["1", "2", "3", "4"]
_STAT_COLS = ["na", "na(%)", "nb", "nb(%)", "SD", "MD", "DC", "V", "Q"]


def _fixture():
    codes = load_golden("seqrep_input_codes").astype(int)
    w = load_golden("seqrep_input_weights").ravel()
    chars = np.array([[_ALPH[c] for c in row] for row in codes], dtype=object)
    return seqdef(chars, alphabet=_ALPH, weights=w)


def _group():
    return load_golden("seqrep_grp_input_group").ravel().astype(int)


def _check(weighted, suff):
    x = _fixture()
    D = load_golden("seqrep_input_diss")
    stats = seqrep_grp(x, _group(), diss=D, ret="stat", criterion="density", weighted=weighted)
    keys = sorted(stats)  # group levels in sorted order -> g1, g2
    assert len(keys) == 2
    for k, key in enumerate(keys, start=1):
        got = stats[key][_STAT_COLS].to_numpy(dtype=float)
        exp = load_golden(f"seqrep_grp_{suff}_g{k}_stats")
        np.testing.assert_allclose(got, exp, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqrep_grp_weighted():
    _check(True, "w")


def test_seqrep_grp_unweighted():
    _check(False, "uw")


def test_seqrep_grp_requires_two_per_group():
    import pytest

    x = _fixture()
    D = load_golden("seqrep_input_diss")
    bad = np.array([1] + [2] * 9)  # group 1 has a single case
    with pytest.raises(ValueError, match="2 or more"):
        seqrep_grp(x, bad, diss=D)


def test_seqrep_grp_ret_modes():
    x = _fixture()
    D = load_golden("seqrep_input_diss")
    reps = seqrep_grp(x, _group(), diss=D, ret="rep")
    both = seqrep_grp(x, _group(), diss=D, ret="both")
    assert set(reps) == {"1", "2"}
    assert "rep_indices" in reps["1"]
    assert isinstance(both, tuple) and len(both) == 2


def test_seqrep_grp_internal_ragged_distance_matches_sts_precomputed_route():
    rows = [
        ["C", "B", "B", "B"],
        ["A", "C"],
        ["C", "B"],
        ["B"],
        ["C", "C", "C", "B"],
        ["A", "A", "B", "B", "C", "A"],
        ["A", "A", "C"],
        ["B", "A", "A"],
        ["B"],
        ["B", "B", "A", "B", "A", "C"],
    ]
    x = seqdef(
        rows,
        alphabet=["A", "B", "C"],
        weights=np.asarray([1.0, 2.5, 1.5, 3.0, 1.0, 2.0, 1.0, 4.0, 2.0, 1.0]),
    )
    group = np.asarray([1, 1, 1, 1, 1, 2, 2, 2, 2, 2])
    aware_distance = seqdist(x, method="LCS")

    internal = seqrep_grp(
        x,
        group,
        ret="rep",
        criterion="dist",
        weighted=True,
        method="LCS",
    )
    precomputed = seqrep_grp(
        x,
        group,
        diss=aware_distance,
        ret="rep",
        criterion="dist",
        weighted=True,
    )

    assert internal.keys() == precomputed.keys()
    for key in internal:
        np.testing.assert_array_equal(internal[key]["rep_indices"], precomputed[key]["rep_indices"])
        np.testing.assert_array_equal(
            internal[key]["representatives"], precomputed[key]["representatives"]
        )
        np.testing.assert_allclose(
            internal[key]["statistics"][_STAT_COLS],
            precomputed[key]["statistics"][_STAT_COLS],
            rtol=0,
            atol=0,
        )
