"""Regression coverage for audit findings TX-AUD-20261005-022 and -023.

The four-case default and nrep=2 expectations were independently recorded by
the pinned R oracle in evidence/final-ledger-boundary-R-v2.log. The LCS matrix
below follows 2 * (3 - longest common subsequence length), by hand.
"""

import numpy as np
import pytest

from traminex import seqdef, seqrep_grp

DISS = np.array([[0, 2, 4, 4], [2, 0, 2, 2], [4, 2, 0, 2], [4, 2, 2, 0]], float)
STAT_COLUMNS = ["na", "na(%)", "nb", "nb(%)", "SD", "MD", "DC", "V", "Q"]


@pytest.fixture
def sequences():
    return seqdef([list(row) for row in ("AAB", "ABB", "BBA", "BBB")], alphabet=["A", "B"])


@pytest.mark.parametrize("group_kwargs", [{}, {"group": None}], ids=["omitted", "None"])
def test_default_group_contains_every_case(sequences, group_kwargs):
    stats = seqrep_grp(sequences, diss=DISS, **group_kwargs)
    assert list(stats) == ["1"]
    np.testing.assert_allclose(
        stats["1"][STAT_COLUMNS],
        [[4, 100, 1, 25, 10, 2.5, 4, 1, -150], [4, 100, 1, 25, 10, 2.5, 4, 1, -150]],
        rtol=0,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    "group",
    ["g", ["g", "g"], ["g"] * 5, [["g"] * 4], [["g"], ["g"], ["g"], ["g"]]],
    ids=["scalar", "short", "long", "row-matrix", "column-matrix"],
)
def test_group_must_be_one_dimensional_and_match_case_count(sequences, group):
    with pytest.raises(ValueError, match="group|grp"):
        seqrep_grp(sequences, group, diss=DISS)


@pytest.mark.parametrize("computed", [False, True], ids=["supplied", "computed"])
def test_fixed_representative_count_is_forwarded(sequences, computed):
    distances = {"method": "LCS"} if computed else {"diss": DISS}
    reps = seqrep_grp(sequences, ["g"] * 4, nrep=2, ret="rep", **distances)["g"]
    np.testing.assert_array_equal(reps["rep_indices"], [0, 1])
    np.testing.assert_array_equal(reps["representatives"], [[0, 0, 1], [0, 1, 1]])
    assert reps["quality"] == pytest.approx(0, abs=1e-12)


@pytest.mark.parametrize("computed", [False, True], ids=["supplied", "computed"])
@pytest.mark.parametrize("decreasing,want", [(True, [1]), (False, [3])])
def test_explicit_score_and_direction_select_requested_trajectory(
    sequences, computed, decreasing, want
):
    distances = {"method": "LCS"} if computed else {"diss": DISS}
    reps = seqrep_grp(
        sequences,
        ["g"] * 4,
        criterion="prob",
        score=[10, 30, 20, 0],
        decreasing=decreasing,
        nrep=1,
        ret="rep",
        **distances,
    )["g"]
    np.testing.assert_array_equal(reps["rep_indices"], want)
    np.testing.assert_array_equal(reps["scores"], [10, 30, 20, 0])


@pytest.mark.parametrize("computed", [False, True], ids=["supplied", "computed"])
def test_neighborhood_radius_changes_density_representative(sequences, computed):
    distances = {"method": "LCS"} if computed else {"diss": DISS}
    reps = seqrep_grp(sequences, ["g"] * 4, pradius=0.6, ret="rep", **distances)["g"]
    # Radius 2.4 includes all four cases only around ABB.
    np.testing.assert_array_equal(reps["rep_indices"], [1])
    assert reps["statistics"].loc["Total", "nb"] == 4
    assert reps["metadata"]["pradius"] == pytest.approx(2.4)


@pytest.mark.parametrize("computed", [False, True], ids=["supplied", "computed"])
def test_unknown_keyword_is_rejected_instead_of_discarded(sequences, computed):
    distances = {"method": "LCS"} if computed else {"diss": DISS}
    with pytest.raises((TypeError, ValueError), match="unknown|unexpected|unsupported|nrepp"):
        seqrep_grp(sequences, ["g"] * 4, nrepp=2, **distances)


@pytest.mark.parametrize("computed", [False, True], ids=["supplied", "computed"])
def test_user_dmax_cannot_replace_shared_group_normalization(sequences, computed):
    distances = {"method": "LCS"} if computed else {"diss": DISS}
    with pytest.raises((TypeError, ValueError), match="dmax"):
        seqrep_grp(sequences, ["g"] * 4, dmax=100, **distances)


def test_missing_distance_control_and_representative_count_are_both_honored():
    sequences = seqdef(
        [["A", None, "B"], ["A", "B", "B"], ["B", None, "A"], ["B", "B", "B"]],
        alphabet=["A", "B"],
    )
    with pytest.raises(ValueError, match="with_missing"):
        seqrep_grp(sequences, ["g"] * 4, method="LCS", with_missing=False)
    reps = seqrep_grp(sequences, ["g"] * 4, method="LCS", with_missing=True, nrep=2, ret="rep")["g"]
    assert reps["rep_indices"].size == 2
