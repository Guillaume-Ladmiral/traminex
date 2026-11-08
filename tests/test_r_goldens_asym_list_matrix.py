"""A substitution matrix given as nested lists is the same matrix as an array (TX-AUD-20261019-003).

R golden: ``scripts/r_goldens/make_goldens_asym_list_matrix.R`` (``asym_list_matrix.csv``). On
sequences of different lengths R's OM applies ``sm[i, j]`` with ``i`` the first sequence's state
and warns "Non-symmetric substitution cost matrix!". traminex applied a nested-list matrix the
other way round on that route, without the warning.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex as tx
from traminex.distances import pairwise, pairwise_between_sets

GOLDEN = pd.read_csv(R_GOLDENS_DIR / "asym_list_matrix.csv")
ROWS = [
    ["A", "B", "C", None],
    ["B", "A", None, None],
    ["A", "B", "C", None],
    ["C", "C", "B", "A"],
    ["B", "A", None, None],
]
SM = np.array([[0, 1, 5], [3, 0, 1], [1, 4, 0]], dtype=float)
INDELS = {"indel_scalar": 2.0, "indel_vector": [1.0, 2.0, 1.5]}
SPELLINGS = {
    "array": lambda m: m,
    "nested_list": lambda m: m.tolist(),
    "tuple": lambda m: tuple(tuple(row) for row in m.tolist()),
}


def _r(case: str) -> tuple[np.ndarray, str]:
    rows = GOLDEN[GOLDEN["case"] == case]
    matrix = np.zeros((5, 5))
    for _, row in rows.iterrows():
        matrix[int(row["i"]) - 1, int(row["j"]) - 1] = float(row["distance"])
    return matrix, str(rows["warning"].iloc[0])


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ at\ least\ one\ indel\ cost\ does\ not\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("case", sorted(INDELS))
@pytest.mark.parametrize("spelling", sorted(SPELLINGS))
def test_seqdist_applies_r_direction_and_warns_for_every_spelling(case: str, spelling: str) -> None:
    expected, message = _r(case)
    assert message == "Non-symmetric substitution cost matrix!"
    x = tx.seqdef(ROWS, alphabet=["A", "B", "C"], right="DEL")
    indel = INDELS[case]
    indel_arg = np.asarray(indel) if isinstance(indel, list) else indel
    with pytest.warns(UserWarning, match=message):
        actual = tx.seqdist(x, method="OM", subm=SPELLINGS[spelling](SM), indel=indel_arg)
    np.testing.assert_allclose(np.asarray(actual, float), expected, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize("spelling", sorted(SPELLINGS))
def test_low_level_pairwise_void_route_reads_every_spelling_alike(spelling: str) -> None:
    expected, message = _r("indel_scalar")
    x = tx.seqdef(ROWS, alphabet=["A", "B", "C"], right="DEL")
    with pytest.warns(UserWarning, match=message):
        actual = pairwise(
            np.asarray(x.data),
            method="OM",
            subm=SPELLINGS[spelling](SM),
            indel=2.0,
            void=x.void_code,
        )
    np.testing.assert_allclose(np.asarray(actual, float), expected, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize("spelling", ["nested_list", "tuple"])
def test_between_sets_ragged_route_reads_a_nested_list_as_the_array(spelling: str) -> None:
    # Batch-7 review I-3: rows of different lengths reach the row engine directly. The
    # contract here is traminex's (one matrix, however spelled); how the two sets are
    # oriented against R's refseq=list(a, b) is finding TX-AUD-20261019-009 (batch 8).
    x = tx.seqdef(ROWS, alphabet=["A", "B", "C"], right="DEL")
    rows = [row[row != x.void_code] for row in np.asarray(x.data)]
    with pytest.warns(UserWarning, match="Non-symmetric substitution cost matrix!"):
        from_array = pairwise_between_sets(rows[:2], rows[2:], method="OM", subm=SM, indel=2.0)
    with pytest.warns(UserWarning, match="Non-symmetric substitution cost matrix!"):
        from_spelling = pairwise_between_sets(
            rows[:2], rows[2:], method="OM", subm=SPELLINGS[spelling](SM), indel=2.0
        )
    np.testing.assert_array_equal(np.asarray(from_spelling), np.asarray(from_array))


@pytest.mark.parametrize("spelling", ["nested_list", "tuple"])
def test_dense_route_reads_a_nested_list_as_the_array(spelling: str) -> None:
    dense = [["A", "B", "C"], ["B", "A", "A"], ["C", "C", "B"]]
    x = tx.seqdef(dense, alphabet=["A", "B", "C"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from_array = tx.seqdist(x, method="OM", subm=SM, indel=2.0)
        from_spelling = tx.seqdist(x, method="OM", subm=SPELLINGS[spelling](SM), indel=2.0)
    np.testing.assert_array_equal(np.asarray(from_spelling), np.asarray(from_array))
