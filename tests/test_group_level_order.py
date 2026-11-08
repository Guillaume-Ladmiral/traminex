"""Group-level order of ``dissassoc``, ``disscenter`` and ``dissvar_grp``.

Audit 2026-09-22, finding S4-11. These three functions took their group levels
from ``np.unique``, so a pandas categorical ``low < mid < high`` came out as
``high, low, mid``, while ``seqdiff`` already honoured the declared order.
TraMineR 2.2-14 and TraMineRextras 0.6.9 build ``factor()`` from the group
(``dissassocweighted.R:294``, ``disscenter.R:47``, ``dissvar-grp.R:22``), which
keeps a factor's level order and drops unused levels. Expected values come from
live R through ``scripts/r_goldens/make_goldens_audit20260922_s4.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import dissvar_grp
from traminex.diss import dissassoc, disscenter

FIXTURE = Path(__file__).with_name("r_goldens") / "group_level_order.txt"


def _rows() -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        key, value = line.split("\t")
        rows[key] = value.split(" ")
    return rows


ROWS = _rows()
D = np.array([float(v) for v in ROWS["D"]]).reshape(8, 8)
WEIGHTS = np.array([float(v) for v in ROWS["weights"]])
CATEGORICAL = pd.Categorical(ROWS["group"], categories=ROWS["categories"])
GROUPS = {"Categorical": CATEGORICAL, "Series": pd.Series(CATEGORICAL)}


def _numbers(key: str) -> np.ndarray:
    return np.array([float(v) for v in ROWS[key]])


@pytest.mark.parametrize("kind", sorted(GROUPS))
def test_dissassoc_groups_follow_the_declared_order(kind):
    result = dissassoc(D, GROUPS[kind], weights=WEIGHTS, R=1, weight_permutation="diss")
    assert list(result.groups.index) == ROWS["dissassoc_rows"]
    expected = _numbers("dissassoc_groups").reshape(2, -1).T  # R: column-major
    np.testing.assert_allclose(result.groups.to_numpy(), expected, rtol=1e-12, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize("kind", sorted(GROUPS))
def test_disscenter_follows_the_declared_order(kind):
    allcenter = disscenter(D, group=GROUPS[kind], allcenter=True, weights=WEIGHTS)
    assert list(allcenter.columns) == ROWS["allcenter_cols"]
    expected = _numbers("allcenter").reshape(len(ROWS["allcenter_cols"]), -1).T
    np.testing.assert_allclose(allcenter.to_numpy(), expected, rtol=1e-12, atol=R_GOLDEN_ATOL)

    medoids = disscenter(D, group=GROUPS[kind], medoids_index="first", weights=WEIGHTS)
    assert [str(level) for level in medoids] == ROWS["medoid_names"]
    assert [int(index) + 1 for index in medoids.values()] == [int(v) for v in ROWS["medoids"]]


@pytest.mark.parametrize("kind", sorted(GROUPS))
def test_dissvar_grp_follows_the_declared_order(kind):
    result = dissvar_grp(D, GROUPS[kind])
    assert list(result) == ROWS["dissvar_grp_names"]
    np.testing.assert_allclose(
        list(result.values()), _numbers("dissvar_grp"), rtol=1e-12, atol=R_GOLDEN_ATOL
    )
