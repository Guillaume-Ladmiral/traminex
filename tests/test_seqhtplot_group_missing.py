"""seqHtplot(group=...) leaves sequences with a missing group out.

What was wrong (2026-09-22 audit S6-07): a ``NaN`` group value added a line
labelled ``nan`` made only of NaN entropies, a ``None`` value raised
``TypeError`` from ``np.unique``, and a group of the wrong length was not
refused. The sibling ``seqplot_tentrop`` already drops missing groups.

Source of truth: TraMineR 2.2-14 ``seqplot.R`` builds the groups with
``TraMineR:::group``, which maps a missing value to ``NA`` and so to no level;
``seqplot.R`` then draws one panel per level. traminex overlays the group
lines on one axis (documented) and keeps R ``factor`` level order and labels
(``tests/r_goldens/plot_group_levels.tsv``, ``factor`` rows); ``group()``'s
text order of numeric levels (audit S6-08) is an R quirk traminex corrects
(TX-REP-20261017-002, ledger row seqHtplot).
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from traminex import seqdef, seqstatd
from traminex.plotting import seqHtplot

FIXTURE = Path(__file__).with_name("r_goldens") / "plot_group_levels.tsv"

ROWS = [
    ["A", "A", "B", "B"],
    ["A", "B", "B", "C"],
    ["C", "C", "B", "A"],
    ["A", "A", "A", "B"],
    ["B", "C", "C", "C"],
    ["A", "C", "B", "B"],
]
# The same vectors as make_goldens_audit20260922_s6.R (section 4).
GROUPS = {
    "numeric_na": np.array([2, np.nan, 1, 2, np.nan, 1]),
    "character_na": ["a", None, "b", "a", "b", None],
}


def _r_factor_levels():
    with FIXTURE.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["case"]: row["levels"].split("|") for row in rows if row["builder"] == "factor"}


R_LEVELS = _r_factor_levels()


@pytest.mark.parametrize("case", sorted(GROUPS))
def test_missing_group_values_draw_no_line(case):
    x = seqdef(ROWS, alphabet=["A", "B", "C"])
    group = GROUPS[case]

    fig, ax = seqHtplot(x, group=group)

    assert [text.get_text() for text in ax.get_legend().get_texts()] == R_LEVELS[case]
    labels = [str(value) for value in group]
    for line, level in zip(ax.lines, R_LEVELS[case]):
        rows = [index for index, label in enumerate(labels) if label in (level, f"{level}.0")]
        subset = seqdef([ROWS[index] for index in rows], alphabet=["A", "B", "C"])
        expected = np.asarray(seqstatd(subset)["Entropy"], dtype=float)
        np.testing.assert_allclose(line.get_ydata(), expected, rtol=0, atol=1e-12)
    assert len(ax.lines) == len(R_LEVELS[case])
    plt.close(fig)


def test_group_must_have_one_value_per_sequence():
    x = seqdef(ROWS, alphabet=["A", "B", "C"])

    with pytest.raises(ValueError, match="one value for each sequence"):
        seqHtplot(x, group=[1, 2, 1])
    plt.close("all")
