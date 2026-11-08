"""Survival plots treat a palette shorter than the alphabet like every other plot.

What was wrong (2026-09-22 audit S6-14): ``STSList`` accepts a ``cpal`` shorter
than its alphabet. The index, distribution and frequency plots then warn and
fall back to the default palette (the 2026-09-02 decision pinned by
``test_ensure_cmap_warns_when_a_declared_palette_is_too_short``), and since
the S6-11 fix the legend helpers do the same, but ``seqsurvplot`` and
``seqsplot`` raised ``IndexError`` from ``seqsurv``'s palette lookup.

Source of truth: that decision -- one shared fallback (``ensure_cmap``) -- so a
state is drawn in the same colour by every renderer. Whether ``STSList`` should
refuse such a palette outright is a maintainer question (audit S6-14, see the
fix report).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import to_rgba

from traminex import seqsurv
from traminex.plotting import seqsplot, seqsurvplot
from traminex.sequences import STSList
from traminex.visualization.palette import _state_color, ensure_cmap

DATA = np.array([[0, 0, 1, 1, 2, 3], [3, 3, 2, 2, 1, 0], [0, 1, 1, 2, 2, 3]])


def _short_palette():
    return STSList(DATA.copy(), alphabet=list("ABCD"), cpal=["red", "blue"])


def _fallback_colours(x):
    with pytest.warns(UserWarning, match="cpal"):
        cmap = ensure_cmap(x, None)
    return {state: to_rgba(_state_color(cmap, code)) for code, state in enumerate("ABCD")}


def test_seqsurv_leaves_a_short_palette_unresolved():
    result = seqsurv(_short_palette())

    assert result.cpal is None


def test_seqsurvplot_falls_back_like_the_other_plots():
    x = _short_palette()
    expected = _fallback_colours(x)

    with pytest.warns(UserWarning, match="cpal"):
        fig, ax = seqsurvplot(x)

    assert {line.get_label(): to_rgba(line.get_color()) for line in ax.lines} == {
        state: expected[state] for state in "ABCD"
    }
    plt.close(fig)


def test_seqsplot_falls_back_like_the_other_plots():
    x = _short_palette()
    expected = _fallback_colours(x)

    with pytest.warns(UserWarning, match="cpal"):
        fig, ax = seqsplot(x)

    assert {line.get_label(): to_rgba(line.get_color()) for line in ax.lines} == {
        state: expected[state] for state in "ABCD"
    }
    plt.close(fig)
