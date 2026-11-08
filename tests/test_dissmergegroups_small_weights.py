"""``dissmergegroups`` with a group whose weighted size is below 1.

Audit 2026-09-22, finding S4-12. The average silhouette width (ASW) cannot be
computed when a group's weighted size is below 1: WeightedCluster 2.0's
``wcClusterQuality`` warns and returns NA. traminex returned the input grouping
with quality NaN and no warning. TraMineR 2.2-14 ``dissmergegroups``
(``SOURCES/TraMineR/R/dissmergegroups.R:124-165``) compares that NA in its
merge loop and stops ("missing value where TRUE/FALSE needed"); when no merge is
attempted (``min.group`` already reached) it returns the grouping with quality
NA and WeightedCluster's warning. Expected outcomes come from live R through
``scripts/r_goldens/make_goldens_audit20260922_s4.R``.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from traminex import dissmergegroups

FIXTURE = Path(__file__).with_name("r_goldens") / "dissmergegroups_small_weights.txt"


def _rows() -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        key, *values = line.split("\t")
        rows[key] = values
    return rows


ROWS = _rows()
POSITIONS = np.array([float(v) for v in ROWS["positions"][0].split(" ")])
D = np.abs(POSITIONS[:, None] - POSITIONS[None, :])
GROUP = np.array([int(v) for v in ROWS["group"][0].split(" ")])
WEIGHTS = np.array([float(v) for v in ROWS["weights"][0].split(" ")])
OPTIONS = {"weights": WEIGHTS, "crit": 0.2, "ref": "max", "small": 0.05}


def test_merging_on_an_undefined_quality_is_refused():
    assert ROWS["loop_runs"][0] == "error"  # R stops
    with pytest.raises(ValueError, match="weighted size below one"):
        dissmergegroups(D, GROUP, min_group=2, **OPTIONS)


def test_no_merge_returns_nan_quality_with_a_warning():
    kind, groups, quality, warned = ROWS["no_loop"]
    assert kind == "value" and quality == "NA" and "less than one observation" in warned
    with pytest.warns(RuntimeWarning, match="weighted size below one"):
        result = dissmergegroups(D, GROUP, min_group=4, **OPTIONS)
    assert result["group"].tolist() == [int(v) for v in groups.split(" ")]
    assert math.isnan(result["quality"])
