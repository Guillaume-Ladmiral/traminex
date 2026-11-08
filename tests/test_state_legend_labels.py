"""State legends name states by their labels and never list void, as in R.

What was wrong (2026-09-22 audit S6-11): ``state_legend``,
``state_legend_figure`` and the ``state_distribution_plot`` legend listed the
internal codebook -- the state *codes* (``A``, ``B``, ...) followed by the missing
and void markers (``*``, ``%``) -- instead of the state labels. Structural void
is never a substantive state (``docs/VISUALIZATION_SERIES_CONTRACT.md``).

Source of truth: TraMineR 2.2-14 ``seqlegend.R`` (``ltext <- attr(seqdata,
"labels")``, plus a ``"missing"`` entry in ``missing.color`` when
``with.missing = "auto"`` and the data hold a missing value) and ``seqplot.R``'s
legend, which adds ``"missing"`` only when the plot uses missing values (index
plots always do; ``type = "d"`` only with ``with.missing = TRUE``). A live R run
on this fixture (evidence ``r/s6-11_seqlegend_text.out``) drew
``Alpha|Beta|Gamma|Delta|missing`` for ``seqlegend`` and index plots,
``Alpha|Beta|Gamma|Delta`` for ``seqplot(type = "d")`` and for data without
missing values, and never a void entry.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

from traminex import seqdef
from traminex.plotting import (
    sequence_index_plot,
    state_distribution_plot,
    state_legend,
    state_legend_figure,
)
from traminex.visualization.palette import _state_color, ensure_cmap

ROWS = [
    ["A", "A", "B", "B", "C", "C"],
    ["A", "B", "B", "C", "C", "C"],
    ["B", "B", None, "C", "D", "D"],
    ["A", "A", "A", "A", "B", "B"],
    ["C", "C", "C", "D", "D", None],
    ["A", None, "B", "B", "C", "C"],
    ["D", "D", "D", "D", "D", "D"],
    ["A", "A", "B", "B", "C", "C"],
]
LABELS = ["Alpha", "Beta", "Gamma", "Delta"]
WITH_MISSING = LABELS + ["missing"]


def _sequences(rows=ROWS):
    return seqdef(rows, alphabet=["A", "B", "C", "D"], labels=LABELS)


def _texts(legend):
    return [text.get_text() for text in legend.get_texts()]


def _expected_colours(x, states):
    canonical = x._canonicalized()
    cmap = ensure_cmap(canonical, None)
    codes = [canonical.alphabet.index(state) for state in states]
    return [to_rgba(_state_color(cmap, code)) for code in codes]


def test_state_legend_uses_labels_and_adds_missing_but_not_void():
    x = _sequences()
    fig, ax = plt.subplots()

    patches, labels = state_legend(x, ax=ax)

    assert labels == WITH_MISSING
    assert _texts(ax.get_legend()) == WITH_MISSING
    assert [to_rgba(patch.get_facecolor()) for patch in patches] == _expected_colours(
        x, ["A", "B", "C", "D", "*"]
    )
    plt.close(fig)


def _without_missing_values():
    # The rows without a missing value, keeping the declared missing code (as
    # a subset of a sequence object does), like R's seqdef on those rows.
    import numpy as np

    from traminex.sequences import STSList

    full = _sequences()
    rows = np.array([0, 1, 3, 6, 7])
    return STSList(
        data=full.data[rows],
        alphabet=list(full.alphabet),
        labels=list(full.labels),
        cpal=list(full.cpal),
        nr=full.nr,
        void=full.void,
        nr_code=full.nr_code,
        void_code=full.void_code,
    )


def test_state_legend_without_missing_values_has_no_missing_entry():
    x = _without_missing_values()
    assert x.nr_code is not None
    fig, ax = plt.subplots()

    _, labels = state_legend(x, ax=ax)

    assert labels == LABELS
    plt.close(fig)

    fig, ax = state_distribution_plot(x, include_missing=True)
    assert _texts(ax.get_legend()) == LABELS
    plt.close(fig)


def test_state_legend_figure_lists_the_same_entries():
    x = _sequences()

    fig, ax = state_legend_figure(x)

    assert _texts(ax.get_legend()) == WITH_MISSING
    plt.close(fig)


def test_index_plot_legend_lists_labels_and_missing():
    fig, ax = sequence_index_plot(_sequences(), legend=True)

    assert _texts(ax.get_legend()) == WITH_MISSING
    plt.close(fig)


def test_distribution_plot_legend_adds_missing_only_when_it_is_plotted():
    x = _sequences()

    fig, ax = state_distribution_plot(x)
    assert _texts(ax.get_legend()) == LABELS
    plt.close(fig)

    fig, ax = state_distribution_plot(x, include_missing=True)
    assert _texts(ax.get_legend()) == WITH_MISSING
    plt.close(fig)


def test_distribution_plot_still_draws_one_bar_series_per_codebook_state():
    # The legend changes; the drawn series (one per codebook column) do not.
    x = _sequences()

    fig, ax = state_distribution_plot(x)

    assert len(ax.containers) == len(x.alphabet)
    plt.close(fig)


def test_short_palette_falls_back_in_the_legend_as_in_the_plots():
    # 2026-09-22 audit S6-14: with a cpal shorter than the alphabet the plots
    # warned and fell back to the default palette, while state_legend raised
    # IndexError and state_legend_figure raised ValueError.
    import numpy as np
    import pytest

    from traminex.sequences import STSList

    x = STSList(np.array([[0, 1, 2, 3], [3, 2, 1, 0]]), alphabet=list("ABCD"), cpal=["red", "blue"])
    with pytest.warns(UserWarning, match="cpal"):
        fallback = ensure_cmap(x, None)
    expected = [to_rgba(_state_color(fallback, code)) for code in range(4)]

    fig, ax = plt.subplots()
    with pytest.warns(UserWarning, match="cpal"):
        patches, _ = state_legend(x, ax=ax)
    assert [to_rgba(patch.get_facecolor()) for patch in patches] == expected
    plt.close(fig)

    with pytest.warns(UserWarning, match="cpal"):
        fig, ax = state_legend_figure(x)
    assert [to_rgba(patch.get_facecolor()) for patch in ax.get_legend().get_patches()] == expected
    plt.close(fig)
