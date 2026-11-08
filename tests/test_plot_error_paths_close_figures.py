"""A plotting call that fails closes the figures it opened.

What was wrong (2026-09-22 audit S6-09): several renderers created their figure
before validating arguments or applying caller styles, so a failing call left
an open figure in pyplot's registry. Twenty failing ``plot_seqimplic`` calls in a
loop left twenty figures open -- a leak in notebooks and batch jobs.

Source of truth: the renderers' contract -- a call that raises must not change
pyplot's figure registry, except for figures the caller supplied, which stay
open. The cases are the audit's reproductions (``probe04_axes_leaks.py``).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import traminex
from traminex import seqdef, seqdist
from traminex.plotting import (
    plot_SAMM,
    plot_seqimplic,
    seqplot_rf,
    seqplot_tentrop,
    seqrfplot,
    seqsplot,
    sequence_index_plot,
)
from traminex.stats_p2 import seqimplic

ROWS = [
    ["A", "A", "B", "B", "C", "C"],
    ["A", "B", "B", "C", "C", "C"],
    ["B", "B", "B", "C", "D", "D"],
    ["A", "A", "A", "A", "B", "B"],
    ["C", "C", "C", "D", "D", "D"],
    ["A", "A", "B", "B", "C", "C"],
    ["D", "D", "D", "D", "D", "D"],
    ["A", "A", "B", "B", "C", "C"],
]
GROUP = np.array(["g1", "g2", "g1", "g2", "g1", "g2", "g1", "g2"])


@pytest.fixture(scope="module")
def data():
    x = seqdef(ROWS, alphabet=["A", "B", "C", "D"])
    return {
        "x": x,
        "D": seqdist(x, method="LCS"),
        "implic": seqimplic(x, GROUP),
        "samm": traminex.seqsamm(x, 2),
    }


FAILING_CALLS = [
    ("plot_seqimplic xtstep=0", lambda d: plot_seqimplic(d["implic"], xtstep=0), ValueError),
    ("plot_seqimplic axes='bad'", lambda d: plot_seqimplic(d["implic"], axes="bad"), ValueError),
    ("plot_seqimplic short ltext", lambda d: plot_seqimplic(d["implic"], ltext=["a"]), ValueError),
    (
        "seqplot_tentrop decreasing ylim",
        lambda d: seqplot_tentrop(d["x"], GROUP, ylim=(1, 0)),
        ValueError,
    ),
    ("seqplot_tentrop xtstep=0", lambda d: seqplot_tentrop(d["x"], GROUP, xtstep=0), ValueError),
    (
        "seqplot_tentrop short xtlab",
        lambda d: seqplot_tentrop(d["x"], GROUP, xtlab=["a"]),
        ValueError,
    ),
    ("seqsplot short xtlab", lambda d: seqsplot(d["x"], xtlab=["a"]), ValueError),
    (
        "sequence_index_plot text highlight",
        lambda d: sequence_index_plot(d["x"], highlight=["x"]),
        TypeError,
    ),
    (
        "seqrfplot bad artist style",
        lambda d: seqrfplot(d["x"], d["D"], k=2, artist_style={"hatch": 1.5j}),
        TypeError,
    ),
    (
        "seqplot_rf unknown style",
        lambda d: seqplot_rf(d["x"], d["D"], k=2, bogus_property=1),
        AttributeError,
    ),
    ("plot_SAMM r without diss", lambda d: plot_SAMM(d["samm"], type="r"), ValueError),
    (
        "plot_SAMM rf wrong diss",
        lambda d: plot_SAMM(d["samm"], type="rf", diss=np.zeros((1, 1))),
        ValueError,
    ),
]


@pytest.mark.parametrize(
    ("name", "call", "error"), FAILING_CALLS, ids=[case[0] for case in FAILING_CALLS]
)
def test_failed_call_leaves_no_new_figure_open(data, name, call, error):
    plt.close("all")

    with pytest.raises(error):
        call(data)

    assert plt.get_fignums() == []


def test_repeated_failures_do_not_accumulate_figures(data):
    plt.close("all")
    for _ in range(20):
        with pytest.raises(ValueError):
            plot_seqimplic(data["implic"], xtstep=0)
    assert plt.get_fignums() == []


def test_a_caller_supplied_figure_stays_open_after_a_failure(data):
    plt.close("all")
    fig, ax = plt.subplots()

    with pytest.raises(ValueError):
        seqsplot(data["x"], xtlab=["a"], ax=ax)

    assert plt.get_fignums() == [fig.number]
    plt.close(fig)
