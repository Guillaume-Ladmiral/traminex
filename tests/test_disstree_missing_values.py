"""``disstree`` with missing or infinite dissimilarities and weights, as in R.

Audit 2026-09-22, finding S4-16. TraMineR 2.2-14 compares candidate splits
with ``if (SCres > currentSCres)`` (``SOURCES/TraMineR/R/disstree.R:476``). A
missing dissimilarity between two cases, a missing weight, or a predictor level
whose cases all weigh 0 makes a node's residual discrepancy NaN, so that
comparison is NA and R stops ("missing value where TRUE/FALSE needed").
traminex instead skipped the predictor and silently returned a tree with the
root alone.

When no comparison is made (``max_depth=1``, or a predictor with a single
level) R returns the root alone with a NaN discrepancy. R's C kernels never
read the diagonal (their loops run over ``j > i``), so a missing value there is
ignored, where traminex's sums made the whole tree NaN. Infinite
dissimilarities are ordinary arithmetic in both. Expected outcomes come from
live R through ``scripts/r_goldens/make_goldens_audit20260922_s4b.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import disstree

FIXTURE = Path(__file__).with_name("r_goldens") / "disstree_missing_values.txt"


def _cases() -> dict[str, dict[str, list[str]]]:
    cases: dict[str, dict[str, list[str]]] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        case, key, *values = line.split("\t")
        cases.setdefault(case, {})[key] = values
    return cases


CASES = _cases()
CASES.pop("pos")


def _numbers(text: str) -> np.ndarray:
    return np.array([np.nan if token == "NA" else float(token) for token in text.split(" ")])


def _fit(case: str):
    row = CASES[case]
    D = _numbers(row["D"][0])
    D = D.reshape(int(round(np.sqrt(D.size))), -1)
    weights = None if row["weights"][0] == "none" else _numbers(row["weights"][0])
    return disstree(
        D,
        pd.DataFrame({"v": row["v"][0].split(" ")}),
        weights=weights,
        R=1,
        min_size=0,
        max_depth=int(row["max_depth"][0]),
        weight_permutation="diss",
    )


REFUSED = sorted(case for case, row in CASES.items() if row["outcome"][0] == "error")
GROWN = sorted(case for case, row in CASES.items() if row["outcome"][0] == "tree")


def test_the_fixture_covers_both_outcomes():
    assert len(REFUSED) == 5 and len(GROWN) == 4


@pytest.mark.parametrize("case", REFUSED)
def test_undefined_split_comparison_is_refused(case):
    assert CASES[case]["outcome"][1] == "missing value where TRUE/FALSE needed"
    with pytest.raises(ValueError, match="residual discrepancy is NaN"):
        _fit(case)


@pytest.mark.parametrize("case", GROWN)
def test_tree_matches_r(case):
    row = CASES[case]
    tree = _fit(case)
    assert tree.fitted.tolist() == [int(v) for v in row["fitted"][0].split(" ")]
    np.testing.assert_allclose(
        tree.root.info["vardis"], _numbers(row["root_vardis"][0])[0], rtol=1e-12
    )
    np.testing.assert_allclose(
        tree.adjustment.stat["t0"].to_numpy(),
        _numbers(row["adjustment_t0"][0]),
        rtol=1e-12,
        atol=R_GOLDEN_ATOL,
    )


@pytest.mark.parametrize("case", ["nan_off_diagonal_depth_one", "nan_off_diagonal_one_level"])
def test_a_root_only_tree_with_a_nan_discrepancy_is_quiet(case):
    """S4-18: the split schedule averaged tied ranks over ``np.unique(vals)``,
    which contains NaN; ``vals == nan`` selects nothing, so two NumPy
    ``RuntimeWarning``s ("Mean of empty slice") leaked where R prints nothing.
    R's ``rank`` (``na.last = TRUE``) gives each NaN its own trailing rank, so
    the schedule was already right; only the warnings were wrong."""
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        tree = _fit(case)
    assert tree.root.kids is None
    assert tree.root.info["splitschedule"] == 1
