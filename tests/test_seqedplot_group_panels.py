"""seqedplot builds its group panels with R's ``factor(group)``.

What was wrong (2026-09-22 audit S6-06): ``seqedplot(group=...)`` turned every
value into text and sorted the text, so numeric groups 1, 2, 10 came out as
panels ``1, 10, 2``; missing values (``None``, ``NaN``) became panels named
``None`` and ``nan``; and supplying fewer ``axes`` than groups silently dropped
the extra panels. The docstring promised "R factor level order".

Source of truth: TraMineRextras 0.6.9 ``seqedplot.R:53-66`` -- ``group <-
factor(group)`` and one panel per ``levels(group)``, titled by the level.
``factor`` sorts numbers numerically and leaves missing values out of every
level. Expected levels are live R output in ``tests/r_goldens/plot_group_levels.tsv``
(``scripts/r_goldens/make_goldens_audit20260922_s6.R``).
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from traminex.events import seqecreate
from traminex.events._duration import (
    _seqedplot_ages_range,
    _seqedplot_hazard,
    _subset_seqelist,
)
from traminex.visualization.plots import seqedplot

FIXTURE = Path(__file__).with_name("r_goldens") / "plot_group_levels.tsv"

# The same group vectors as make_goldens_audit20260922_s6.R (section 4).
GROUPS = {
    "numeric": [2, 10, 2, 10, 1, 1],
    # A float array, as a numeric group with missing values usually is: R
    # still titles its levels "1" and "2", not "1.0" and "2.0".
    "numeric_na": np.array([2, np.nan, 1, 2, np.nan, 1]),
    "character_na": ["a", None, "b", "a", "b", None],
    "factor_declared": pd.Categorical(["b", "a", "b", None, "a", "b"], categories=["c", "b", "a"]),
}


def _r_factor_levels():
    with FIXTURE.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["case"]: row["levels"].split("|") for row in rows if row["builder"] == "factor"}


R_LEVELS = _r_factor_levels()


def _events():
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4", "s4", "s5", "s5"]
            + ["s6", "s6"],
            "timestamp": [1, 3, 5, 2, 4, 1, 6, 1, 2, 4, 3, 5, 5, 2],
            "event": ["A", "A", "B", "A", "B", "B", "A", "A", "A", "A", "B", "B", "A", "B"],
        }
    )
    return seqecreate(tse)


@pytest.mark.parametrize("case", sorted(GROUPS))
def test_panels_follow_r_factor_levels(case):
    fig, axes = seqedplot(_events(), type="hazard", group=GROUPS[case], breaks=5)

    assert [ax.get_title() for ax in axes] == R_LEVELS[case]
    plt.close(fig)


def test_missing_group_values_are_left_out_of_every_panel():
    events = _events()
    group = GROUPS["numeric_na"]

    fig, axes = seqedplot(events, type="hazard", group=group, breaks=5)

    ages = _seqedplot_ages_range(events, None)
    for ax, level in zip(axes, (1, 2)):
        rows = [index for index, value in enumerate(group) if value == level]
        _, series = _seqedplot_hazard(_subset_seqelist(events, rows), None, ages, 5)
        assert [line.get_ydata().tolist() for line in ax.lines] == [
            np.asarray(values, dtype=float).tolist() for _, values in series
        ]
    plt.close(fig)


def test_supplied_axes_must_match_the_panel_count():
    fig, ax = plt.subplots()

    with pytest.raises(ValueError, match="3 group panels"):
        seqedplot(_events(), type="hazard", group=GROUPS["numeric"], breaks=5, axes=[ax])
    plt.close(fig)


def test_group_must_have_one_value_per_sequence():
    with pytest.raises(ValueError, match="one value for each event sequence"):
        seqedplot(_events(), type="hazard", group=[1, 2], breaks=5)
    plt.close("all")
