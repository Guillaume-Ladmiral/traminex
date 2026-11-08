"""A user palette colours the real states; missing cells keep the missing colour (TX-AUD-20261019-020).

R appends ``missing.color`` after a user ``cpal`` when the data hold the nonresponse code and
colours the states of ``c(alphabet, nr)`` with ``c(cpal, missing.color)`` in order
(plot.stslist.R:83-90; seqplot.R:433-443): with one colour per real state, missing cells are
drawn in the missing colour; with one colour more, in that colour. Void cells match no state and
stay blank. traminex recycled the user colours over every code, so missing cells took the first
state's colour.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402
from traminex.visualization.palette import ensure_cmap  # noqa: E402

USER = ["red", "blue"]
MISSING = mcolors.to_hex("darkgrey")


def _data():
    return tx.seqdef(
        np.array([["A", None, "B"], ["A", None, "B"], ["B", "A", "A"]], dtype=object),
        alphabet=["A", "B"],
    )


def _handles(legend) -> list:
    """The legend's handles; matplotlib 3.6.0, the declared floor, names them ``legendHandles``."""
    return legend.legend_handles if hasattr(legend, "legend_handles") else legend.legendHandles


def test_seqfplot_draws_missing_cells_in_the_missing_colour() -> None:
    fig, ax = tx.seqfplot(_data(), cpal=USER)
    cells = [mcolors.to_hex(patch.get_facecolor()) for patch in ax.patches[:3]]
    plt.close(fig)
    assert cells == [mcolors.to_hex("red"), MISSING, mcolors.to_hex("blue")]


def test_index_plot_draws_missing_cells_in_the_missing_colour() -> None:
    fig, ax = tx.sequence_index_plot(_data(), cmap=USER)
    image = ax.images[0]
    rgba = image.cmap(image.norm(image.get_array()))
    plt.close(fig)
    assert mcolors.to_hex(rgba[0, 1]) == MISSING


def test_distribution_legend_and_modal_bars_use_the_missing_colour() -> None:
    fig, ax = tx.state_distribution_plot(_data(), include_missing=True, cmap=USER)
    handles = [mcolors.to_hex(h.get_facecolor()) for h in _handles(ax.get_legend())]
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    plt.close(fig)
    assert dict(zip(labels, handles, strict=True))["missing"] == MISSING
    fig, ax = tx.seqmsplot(_data(), include_missing=True, cpal=USER)
    bars = [mcolors.to_hex(patch.get_facecolor()) for patch in ax.patches]
    plt.close(fig)
    assert bars[1] == MISSING  # position 2: missing is the modal state


def test_one_colour_more_than_the_states_colours_the_missing_cells() -> None:
    canonical = _data()._canonicalized()
    colours = ensure_cmap(canonical, ["red", "blue", "green"]).colors
    assert mcolors.to_hex(colours[canonical.nr_code]) == mcolors.to_hex("green")


def test_void_cells_stay_white() -> None:
    ragged = tx.seqdef(
        np.array([["A", "B", None], ["B", "A", "A"]], dtype=object), alphabet=["A", "B"]
    )
    canonical = ragged._canonicalized()
    colours = ensure_cmap(canonical, USER).colors
    assert canonical.void_code is not None
    assert mcolors.to_hex(colours[canonical.void_code]) == mcolors.to_hex("white")


@pytest.mark.parametrize("cmap", [None, "viridis"])
def test_the_object_palette_and_named_colormaps_are_unchanged(cmap) -> None:
    canonical = _data()._canonicalized()
    assert len(ensure_cmap(canonical, cmap).colors) == len(canonical.alphabet)


def _legend_colours(ax) -> dict[str, str]:
    legend = ax.get_legend()
    return {
        text.get_text(): mcolors.to_hex(handle.get_facecolor())
        for text, handle in zip(legend.get_texts(), _handles(legend), strict=True)
    }


@pytest.mark.parametrize(
    "colours, missing", [(USER, MISSING), (["red", "blue", "green"], mcolors.to_hex("green"))]
)
def test_index_plot_legends_use_the_colours_of_the_cells(colours, missing) -> None:
    """R's seqIplot legend takes the user cpal (TX-REP-20261021-002).

    Live TraMineR 2.2-14 passes legend() the fill ``red blue darkgrey`` for two colours and
    ``red blue green`` for three; the index plots drew the states in the object's palette.
    """
    expected = {"A": mcolors.to_hex("red"), "B": mcolors.to_hex("blue"), "missing": missing}
    fig, ax = tx.sequence_index_plot(_data(), cmap=colours, legend=True)
    plain = _legend_colours(ax)
    plt.close(fig)
    fig, ax, _ = tx.sequence_index_plot_grouped(_data(), np.array([1, 1, 2]), cmap=colours)
    grouped = _legend_colours(ax)
    plt.close(fig)
    assert plain == expected
    assert grouped == expected


@pytest.mark.parametrize(
    "colours, missing", [(USER, MISSING), (["red", "blue", "green"], mcolors.to_hex("green"))]
)
def test_data_with_void_and_missing_cells_keep_both_marker_colours(colours, missing) -> None:
    """One colour per real state, or one more, on data holding nonresponse and void together.

    R colours ``c(alphabet, nr)`` with ``c(cpal, missing.color)`` and leaves void blank, so
    the real states take the user colours, nr the extra colour or the missing colour, and
    void white (TX-AUD-20261019-020; Review Focus 5 of the batch-8 plan).
    """
    x = tx.seqdef(
        np.array([["A", None, "B", None], ["B", "A", "A", "A"]], dtype=object), alphabet=["A", "B"]
    )
    canonical = x._canonicalized()
    assert canonical.nr_code is not None and canonical.void_code is not None
    hexes = [mcolors.to_hex(colour) for colour in ensure_cmap(canonical, colours).colors]
    assert [hexes[canonical.alphabet.index(state)] for state in ("A", "B")] == [
        mcolors.to_hex("red"),
        mcolors.to_hex("blue"),
    ]
    assert hexes[canonical.nr_code] == missing
    assert hexes[canonical.void_code] == mcolors.to_hex("white")


def test_a_short_list_keeps_one_colour_per_state() -> None:
    """A list shorter than the real states is recycled over them, so each state keeps one colour.

    R's Rd asks for one colour per state (plot.stslist.Rd); given fewer, R recycles
    ``c(cpal, missing.color)`` cell by cell, so ``AAA`` is drawn red, blue, darkgrey (live R
    trace, batch-8 review M-3). traminex's rule for this input is its own.
    """
    x = tx.seqdef(
        np.array([["A", "A", "A"], ["B", "B", "B"], ["C", None, "C"]], dtype=object),
        alphabet=["A", "B", "C"],
    )
    canonical = x._canonicalized()
    hexes = [mcolors.to_hex(colour) for colour in ensure_cmap(canonical, USER).colors]
    red, blue = mcolors.to_hex("red"), mcolors.to_hex("blue")
    assert [hexes[canonical.alphabet.index(state)] for state in "ABC"] == [red, blue, red]
    assert hexes[canonical.nr_code] == MISSING
    fig, ax = tx.sequence_index_plot(x, cmap=USER)
    image = ax.images[0]
    rgba = image.cmap(image.norm(image.get_array()))
    plt.close(fig)
    assert {mcolors.to_hex(cell) for cell in rgba[0]} == {red}
