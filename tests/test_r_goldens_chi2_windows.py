"""CHI2/EUCLID step and breaks (TX-AUD-20261011-048).

R refuses a step that is not a positive integer (seqdist.R:189-190) and documents
breaks as "List of ordered pairs of integers" of start and end positions
(seqdist.Rd:209-214), but does not check them: a non-integer, unordered or
out-of-range pair gives whatever ``seq(start, end)`` happens to select (``(2.7, 5)``
uses positions 2-4, ``(5, 3)`` the window 3-5). traminex truncated each position
with ``int()``, so ``(2.7, 5)`` used 2-5. traminex now refuses every pair outside
the documented contract (a recorded correction of R, as for refusals of other
undocumented inputs) and keeps R's values for the documented ones. R stops on
step = 0 unless breaks replaces it; traminex used the whole sequence, and its
low-level ``pairwise`` truncated a fractional step (batch-6 review M-1, I-3).
Fixture ``scripts/r_goldens/make_goldens_chi2_windows.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.distances import python_kernels as pk

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "chi2_windows.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
S5 = tx.seqdef(
    np.array([list("ABABA"), list("BBAAA"), list("AABAA"), list("BABBB")], dtype=object),
    alphabet=["A", "B"],
)
VALID = {
    "step_2": {"step": 2},
    "step_2_numeric": {"step": 2.0},
    "breaks_2_to_4": {"breaks": [(2, 4)]},
    "breaks_2_to_5": {"breaks": [(2, 5)]},
    "breaks_3_to_5": {"breaks": [(3, 5)]},
    "breaks_1_to_2_and_3_to_5": {"breaks": [(1, 2), (3, 5)]},
}
OUTSIDE_THE_CONTRACT = {
    "breaks_2_7_to_5": {"breaks": [(2.7, 5)]},
    "breaks_5_to_3": {"breaks": [(5, 3)]},
    "breaks_0_to_5": {"breaks": [(0, 5)]},
    "breaks_1_to_2_5_and_2_5_to_5": {"breaks": [(1, 2.5), (2.5, 5)]},
}
R_STOPS = {
    "step_0": {"step": 0},
    "step_minus_1": {"step": -1},
    "step_0_overlap": {"step": 0, "overlap": True},
    "step_2_5_with_breaks_2_to_4": {"step": 2.5, "breaks": [(2, 4)]},
}
CODES = np.asarray(S5.data, dtype=np.int64)
LOW_LEVEL = pytest.mark.parametrize(
    "pairwise", [pk.pairwise, tx.distances.pairwise], ids=["python_kernels", "dispatcher"]
)


def _r(case: str) -> np.ndarray:
    assert GOLDEN.loc[case, "ok"] == "TRUE"
    return np.array([float(v) for v in GOLDEN.loc[case, "value"].split()]).reshape(4, 4, order="F")


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
@pytest.mark.parametrize("name", sorted(VALID))
def test_documented_windows_equal_r(method: str, name: str) -> None:
    actual = tx.seqdist(S5, method=method, **VALID[name])
    np.testing.assert_allclose(actual, _r(f"{method.lower()}_{name}"), atol=1e-12, rtol=0)


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
def test_a_fractional_step_is_refused_as_r_refuses_it(method: str) -> None:
    assert GOLDEN.loc[f"{method.lower()}_step_2_5", "ok"] == "FALSE"
    with pytest.raises(ValueError, match="'step' must be a positive integer"):
        tx.seqdist(S5, method=method, step=2.5)


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
@pytest.mark.parametrize("name", sorted(OUTSIDE_THE_CONTRACT))
def test_windows_outside_the_contract_are_refused(method: str, name: str) -> None:
    _r(f"{method.lower()}_{name}")  # R computes a by-product value; recorded, not matched
    with pytest.raises(ValueError, match=r"breaks\[\d\]"):
        tx.seqdist(S5, method=method, **OUTSIDE_THE_CONTRACT[name])


@pytest.mark.parametrize("pair", [(2.0, 4), (np.int64(2), 4), (2, np.float64(4.0))])
def test_integer_valued_positions_are_accepted(pair) -> None:
    np.testing.assert_allclose(
        tx.seqdist(S5, method="CHI2", breaks=[pair]), _r("chi2_breaks_2_to_4"), atol=1e-12, rtol=0
    )


def test_a_boolean_position_is_refused() -> None:
    with pytest.raises(ValueError, match=r"breaks\[0\]"):
        tx.seqdist(S5, method="CHI2", breaks=[(True, 4)])


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
@pytest.mark.parametrize("name", sorted(R_STOPS))
def test_a_step_r_stops_on_is_refused(method: str, name: str) -> None:
    assert GOLDEN.loc[f"{method.lower()}_{name}", "ok"] == "FALSE"
    with pytest.raises(ValueError, match="'step' must be a positive integer"):
        tx.seqdist(S5, method=method, **R_STOPS[name])


@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
def test_a_zero_step_is_ignored_with_breaks_as_r_ignores_it(method: str) -> None:
    with pytest.warns(UserWarning, match="'step' has been ignored because 'breaks' is specified"):
        actual = tx.seqdist(S5, method=method, step=0, breaks=[(2, 4)])
    np.testing.assert_allclose(
        actual, _r(f"{method.lower()}_step_0_with_breaks_2_to_4"), atol=1e-12, rtol=0
    )


@LOW_LEVEL
@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
@pytest.mark.parametrize("step", [2.5, 0, -1, True])
def test_the_low_level_layer_refuses_these_steps(pairwise, method: str, step) -> None:
    with pytest.raises(ValueError, match="'step' must be a positive integer"):
        pairwise(CODES, method=method, step=step)


@LOW_LEVEL
@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
@pytest.mark.parametrize("step", [2, 2.0, np.int64(2)])
def test_the_low_level_layer_takes_an_integer_valued_step(pairwise, method: str, step) -> None:
    np.testing.assert_allclose(
        pairwise(CODES, method=method, step=step),
        _r(f"{method.lower()}_step_2"),
        atol=1e-12,
        rtol=0,
    )


@LOW_LEVEL
@pytest.mark.parametrize("method", ["CHI2", "EUCLID"])
def test_the_low_level_layer_ignores_a_zero_step_with_breaks(pairwise, method: str) -> None:
    np.testing.assert_allclose(
        pairwise(CODES, method=method, step=0, breaks=[(2, 4)]),
        _r(f"{method.lower()}_step_0_with_breaks_2_to_4"),
        atol=1e-12,
        rtol=0,
    )
