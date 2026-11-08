"""Parity regressions for TraMineR's accepted asymmetric substitution costs."""

import re

import numpy as np
import pytest

import traminex as tx
from traminex.distances import pairwise

#: TraMineR 2.2-14's checkcost warning (TraMineR-checkcost.R:116), recorded by
#: scripts/r_goldens/make_goldens_audit20260922_s10.R; the 2.2-13 text was
#: "The substitution cost matrix is not symmetric." (audit S10-01).
NON_SYMMETRIC = re.escape("Non-symmetric substitution cost matrix!")


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ pos\ \d+,\ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("method", ["OM", "HAM", "DHD"])
@pytest.mark.parametrize(
    "rows,directional_cost",
    [
        (["A", "B", "A"], 1.0),
        (["B", "A", "B"], 2.0),
    ],
)
def test_asymmetric_costs_match_traminer_stable_distinct_sequence_orientation(
    method, rows, directional_cost
):
    """Duplicate rows share the direction chosen by their first distinct representative."""
    sequence = tx.seqdef(np.asarray(rows, dtype=object)[:, None], alphabet=["A", "B"])
    costs = np.array([[0.0, 1.0], [2.0, 0.0]])
    kwargs = {"subm": costs[None, :, :] if method == "DHD" else costs, "norm": "none"}
    if method == "OM":
        kwargs["indel"] = 2.0

    with pytest.warns(UserWarning, match=NON_SYMMETRIC):
        distance = np.asarray(tx.seqdist(sequence, method=method, **kwargs))

    expected = np.array(
        [
            [0.0, directional_cost, 0.0],
            [directional_cost, 0.0, directional_cost],
            [0.0, directional_cost, 0.0],
        ]
    )
    np.testing.assert_array_equal(distance, expected)


def test_asymmetric_position_dependent_om_costs_group_duplicate_rows():
    rows = np.array([[0, 0], [1, 0], [0, 0]], dtype=np.int32)
    costs = np.repeat(np.array([[[0.0, 1.0], [2.0, 0.0]]]), 2, axis=0)

    with pytest.warns(UserWarning, match=NON_SYMMETRIC):
        distance = pairwise(rows, method="OM", indel=2.0, subm_by_pos=costs)

    expected = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    np.testing.assert_array_equal(distance, expected)


@pytest.mark.parametrize("method", ["OM", "HAM", "DHD"])
def test_direct_pairwise_uses_stable_distinct_sequence_orientation(method):
    rows = np.array([[1], [0], [1]], dtype=np.int32)
    costs = np.array([[0.0, 1.0], [2.0, 0.0]])
    kwargs = {"subm": costs[None, :, :] if method == "DHD" else costs}
    if method != "OM":
        kwargs["norm"] = "none"
    if method == "OM":
        kwargs["indel"] = 2.0

    with pytest.warns(UserWarning, match=NON_SYMMETRIC):
        distance = pairwise(rows, method=method, **kwargs)

    expected = np.array([[0.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 0.0]])
    np.testing.assert_array_equal(distance, expected)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("method", ["OM", "HAM"])
def test_explicit_deduplication_preserves_first_distinct_sequence_orientation(method):
    sequence = tx.seqdef(np.asarray(["B", "A", "B"], dtype=object)[:, None], alphabet=["A", "B"])
    costs = np.array([[0.0, 1.0], [2.0, 0.0]])
    kwargs = {"subm": costs, "norm": "none", "deduplicate": True}
    if method == "OM":
        kwargs["indel"] = 2.0

    with pytest.warns(UserWarning, match=NON_SYMMETRIC):
        distance = tx.seqdist(sequence, method=method, **kwargs)

    expected = np.array([[0.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 0.0]])
    np.testing.assert_array_equal(distance, expected)
