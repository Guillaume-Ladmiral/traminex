"""RP7a: seqpcplot (parallel-coordinate event plot) + plot_eseq/plot_seqelist.

Each event sequence is a polyline through its events: x = event rank, y = the
event code (1-based). #polylines == #sequences; polyline i visits Eseq[i].codes.
plot_eseq / plot_seqelist delegate to seqpcplot.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba

from traminex import seqecreate, seqpcfilter
from traminex.visualization.plots.events_viz import plot_eseq, plot_seqelist, seqpcplot


def _syn_eseq():
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"],
            "timestamp": [0, 1, 1, 3, 0, 2, 4, 0, 1, 2],
            "event": ["A", "B", "C", "A", "A", "B", "C", "A", "B", "C"],
        }
    )
    return seqecreate(tse)


def test_seqpcplot_one_polyline_per_event_sequence():
    el = _syn_eseq()
    fig, ax = seqpcplot(el)
    assert len(ax.lines) == len(el)
    for line, s in zip(ax.lines, el.sequences):
        # x is the rank of the event's distinct time, events of one time by
        # code, as R's seqpcplot (tests/test_r_goldens_seqpcplot_positions.py).
        events = sorted(zip(np.asarray(s.times, dtype=float), np.asarray(s.codes, dtype=float)))
        distinct = sorted({time for time, _ in events})
        np.testing.assert_array_equal(line.get_ydata(), [code for _, code in events])
        np.testing.assert_array_equal(
            line.get_xdata(), [distinct.index(time) + 1 for time, _ in events]
        )
    # y ticks == the event dictionary
    assert [t.get_text() for t in ax.get_yticklabels()] == list(el.dictionary)
    plt.close(fig)


def test_plot_seqelist_delegates_to_seqpcplot():
    el = _syn_eseq()
    fig, ax = plot_seqelist(el)
    assert len(ax.lines) == len(el)
    plt.close(fig)


def test_plot_eseq_draws_single_sequence():
    el = _syn_eseq()
    fig, ax = plot_eseq(el[0])
    assert len(ax.lines) == 1
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), np.asarray(el[0].codes, dtype=float))
    plt.close(fig)


def test_seqpcplot_recycles_palette_in_plotted_sequence_order():
    fig, ax = seqpcplot(_syn_eseq(), cpal=("#ff0000", "#0000ff"))

    assert [line.get_color() for line in ax.lines] == [
        to_rgba("#ff0000"),
        to_rgba("#0000ff"),
        to_rgba("#ff0000"),
    ]
    plt.close(fig)


def test_seqpcplot_restarts_palette_after_filtering_in_plotted_order():
    fig, ax = seqpcplot(
        _syn_eseq(),
        filt=seqpcfilter("minfreq", level=0.5),
        cpal=np.array(["orange", "purple"]),
    )

    assert len(ax.lines) == 2
    assert [line.get_color() for line in ax.lines] == [to_rgba("orange"), to_rgba("purple")]
    plt.close(fig)


def test_seqpcplot_uses_caller_owned_axes_and_scalar_color():
    fig, caller_ax = plt.subplots()

    returned_fig, returned_ax = seqpcplot(_syn_eseq(), cpal="black", ax=caller_ax)

    assert returned_fig is fig
    assert returned_ax is caller_ax
    assert [line.get_color() for line in returned_ax.lines] == [to_rgba("black")] * 3
    plt.close(fig)


def test_seqpcplot_rejects_invalid_matplotlib_color():
    with pytest.raises(ValueError, match="not-a-matplotlib-color"):
        seqpcplot(_syn_eseq(), cpal=["not-a-matplotlib-color"])
    plt.close("all")


def test_plot_delegates_preserve_palette():
    el = _syn_eseq()
    fig1, ax1 = plot_eseq(el[0], cpal="green")
    fig2, ax2 = plot_seqelist(el, cpal=["red", "blue"])

    assert [line.get_color() for line in ax1.lines] == [to_rgba("green")]
    assert [line.get_color() for line in ax2.lines] == [
        to_rgba("red"),
        to_rgba("blue"),
        to_rgba("red"),
    ]
    plt.close(fig1)
    plt.close(fig2)


def test_seqpcplot_labels_an_object_without_labels_by_its_alphabet():
    """An ``STSList`` built without labels is drawn with its alphabet on the y axis, as R's seqdef
    labels default to the states (778febfe; TX-AUD-20261019-056)."""
    from traminex import STSList

    states = STSList(data=np.array([[0, 1, 1, 2], [0, 0, 2, 2]]), alphabet=["a", "b", "c"])
    assert states.labels is None

    fig, ax = seqpcplot(states)

    assert [label.get_text() for label in ax.get_yticklabels()] == ["a", "b", "c"]
    plt.close(fig)
