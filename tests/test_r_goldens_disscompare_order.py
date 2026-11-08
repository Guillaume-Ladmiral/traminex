"""dissCompare draws the larger group first, as seqCompare does (TX-AUD-20261011-022).

seqCompare makes the larger group the first because it draws that group's
sample from 1:max(n) (seqCompare.R:143-157); dissCompare keeps the first factor
level first (dissCompare.R:117-131), so with the smaller group first its
"first" sample takes cases of the other group and the sampled statistic depends
on the level order (R: mean LRT 5.8 smaller-first vs 9.6 larger-first over 20
seeds). traminex copied that. It now draws the larger group first, so the same
seed gives the same result for either order. traminex's draws are its own
NumPy stream (ledger dissCompare row), so the values are compared with R's
in distribution only. Fixture: ``scripts/r_goldens/make_goldens_disscompare_order.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

HERE = Path(__file__).resolve().parent / "r_goldens"
GOLDEN = pd.read_csv(HERE / "disscompare_order.csv")
DISS = np.loadtxt(HERE / "disscompare_order_diss.csv", delimiter=",")
SMALLER_FIRST = np.where(np.arange(1, 91) % 3 == 0, "a", "b")  # 30 "a" (first level), 60 "b"
LARGER_FIRST = np.where(SMALLER_FIRST == "a", "z", "b")  # the same split, "b" first


def _lrt(group, seed: int) -> float:
    return float(tx.dissCompare(DISS, group, s=20, seed=seed, stat="LRT")["LRT"].iloc[0])


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_the_group_order_does_not_change_the_sampled_statistic(seed: int) -> None:
    assert _lrt(SMALLER_FIRST, seed) == pytest.approx(_lrt(LARGER_FIRST, seed), abs=1e-12)


def test_r_sampled_statistic_depends_on_the_group_order() -> None:
    smaller = GOLDEN.loc[GOLDEN["order"] == "smaller_first", "LRT"].mean()
    larger = GOLDEN.loc[GOLDEN["order"] == "larger_first", "LRT"].mean()
    assert larger / smaller > 1.4


def _welch_z(a: np.ndarray, b: np.ndarray) -> float:
    """Welch's statistic for the difference of two sample means."""
    return float(
        abs(a.mean() - b.mean()) / np.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size)
    )


def test_traminex_follows_r_larger_first_draws_in_distribution() -> None:
    # traminex's 20 seeded draws must be within Welch z 3 of R's larger-first draws (about a 0.3 % two-sided
    # level) and beyond 3 of R's smaller-first draws (measured 1.93 and 4.59; R's two orders are 7.63 apart). The
    # seeds are fixed, so the test is deterministic (TX-REP-20261019-003).
    ours = np.array([_lrt(SMALLER_FIRST, seed) for seed in range(1, 21)])
    r_smaller = GOLDEN.loc[GOLDEN["order"] == "smaller_first", "LRT"].to_numpy(dtype=float)
    r_larger = GOLDEN.loc[GOLDEN["order"] == "larger_first", "LRT"].to_numpy(dtype=float)
    assert _welch_z(ours, r_larger) < 3.0
    assert _welch_z(ours, r_smaller) > 3.0


def test_unsampled_statistic_is_unchanged() -> None:
    first = tx.dissCompare(DISS, SMALLER_FIRST, s=0, stat="LRT")["LRT"].iloc[0]
    second = tx.dissCompare(DISS, LARGER_FIRST, s=0, stat="LRT")["LRT"].iloc[0]
    assert first == pytest.approx(second, abs=1e-9)
