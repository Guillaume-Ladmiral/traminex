"""Accepted asymmetric costs: identical sequences get identical distances.

TraMineR removes duplicated sequences before its kernels (``seqdist.R:524``),
so with asymmetric substitution costs -- where the direction of a comparison
matters -- every pair is evaluated in first-occurrence order and identical
sequences always get identical distances; ``checkcost`` also warns
"Non-symmetric substitution cost matrix!". traminex's OMloc, OMslen, OMspell,
OMstran and TWED engines compared pairs in row order and did not warn, so two
copies of one sequence could get different distances to a third (2026-09-22
audit S1b-03). Expected values: live TraMineR 2.2-14 via
``scripts/r_goldens/make_goldens_seqdist_asymmetric_duplicates.R``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDENS = Path(__file__).parent / "r_goldens"
SEQUENCES = [list("AABC"), list("CCBA"), list("AABC"), list("BACC"), list("CCBA")]
SM = np.array([[0.0, 1.0, 3.0], [2.0, 0.0, 1.0], [1.0, 4.0, 0.0]])
METHODS = {
    "OMloc": ("OMLOC", {}),
    "OMslen": ("OMSLEN", {"indel": 1.5}),
    "OMspell": ("OMSPELL", {"indel": 1.5}),
    "OMstran": ("OMSTRAN", {"indel": 1.5, "otto": 0.5}),
    "TWED": ("TWED", {"indel": 1.5, "nu": 0.5, "h": 0.5}),
}


def _r_matrices() -> dict[str, np.ndarray]:
    cells: dict[str, dict] = defaultdict(dict)
    with (GOLDENS / "seqdist_asymmetric_duplicates.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            cells[row["method"]][(int(row["row"]), int(row["col"]))] = float(row["value"])
    out = {}
    for method, values in cells.items():
        matrix = np.empty((5, 5))
        for (i, j), value in values.items():
            matrix[i, j] = value
        out[method] = matrix
    return out


R = _r_matrices()


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ \ at\ least,\ one\ substitution\ cost\ doesn't\ respect\ the\ triangle\ inequality\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
@pytest.mark.parametrize("r_method", sorted(METHODS))
def test_duplicated_sequences_get_rs_first_occurrence_distances(r_method):
    method, extra = METHODS[r_method]
    data = tx.seqdef(np.array(SEQUENCES, dtype=object), alphabet=["A", "B", "C"])

    with pytest.warns(UserWarning, match="symmetric"):
        result = np.asarray(tx.seqdist(data, method=method, subm=SM, **extra), dtype=float)

    np.testing.assert_allclose(result, R[r_method], rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(result[0], result[2])
    np.testing.assert_array_equal(result[1], result[4])
