"""A colormap argument gives every renderer the same, distinct state colours.

What was wrong (2026-09-22 audit S6-13): ``state_distribution_plot(cmap=...)``,
``seqfplot(cpal=...)`` and the other bar and band renderers looked state code
``k`` up as entry ``k`` of a named colormap's 256-entry table, so
``cmap="viridis"`` drew six states in six near-identical dark purples
(``#440154`` .. ``#46085c``), while the index plot (fixed by S6-01) spread the
same states over the whole map. The same state therefore changed colour from one
plot to the next. A colour list longer or shorter than the alphabet was split the
same inconsistent way.

Source of truth: ``_state_code_norm`` (S6-01) defines the index plot's mapping --
state ``c`` of ``K`` takes the colour at the centre of the ``c``-th of ``K``
equal bands of a colormap, and entry ``c`` of a ``K``-colour list. Every renderer
must resolve the same colour for the same state.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import to_rgba

from traminex import seqdef, seqtab
from traminex.plotting import seqfplot, sequence_index_plot, state_distribution_plot

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
K = 4


def _sequences():
    return seqdef(ROWS, alphabet=["A", "B", "C", "D"])


def _index_colours(cmap):
    fig, ax = sequence_index_plot(_sequences(), cmap=cmap)
    image = ax.images[0]
    codes = np.asarray(image.get_array())
    rgba = image.to_rgba(codes)
    colours = {int(code): tuple(rgba[codes == code][0]) for code in np.unique(codes)}
    plt.close(fig)
    return colours


def _bar_colours(cmap):
    fig, ax = state_distribution_plot(_sequences(), cmap=cmap)
    colours = {
        code: tuple(container.patches[0].get_facecolor())
        for code, container in enumerate(ax.containers)
    }
    plt.close(fig)
    return {code: colours[code] for code in range(K)}


def _band_colours(cpal):
    x = _sequences()
    fig, ax = seqfplot(x, cpal=cpal)
    table = seqtab(x, idxs=0, normalize=True)
    rows = table.attrs["representatives"].data[:10]
    length = x.seqlen
    colours = {}
    for band, row in enumerate(rows):
        for position, code in enumerate(row):
            colours.setdefault(
                int(code), tuple(ax.patches[band * length + position].get_facecolor())
            )
    plt.close(fig)
    return colours


def _by_code(colours):
    return np.array([colours[code] for code in range(K)], dtype=float)


@pytest.mark.parametrize("cmap", ["viridis", plt.get_cmap("plasma")], ids=["name", "object"])
def test_colormap_spreads_states_and_every_renderer_agrees(cmap):
    base = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    expected = np.array([to_rgba(base((code + 0.5) / K)) for code in range(K)])

    for colours in (_index_colours(cmap), _bar_colours(cmap), _band_colours(cmap)):
        np.testing.assert_allclose(_by_code(colours), expected, rtol=0, atol=1e-12)
    assert len({tuple(row) for row in expected}) == K


def test_a_longer_colour_list_gives_state_k_entry_k_everywhere():
    palette = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"]
    expected = np.array([to_rgba(palette[code]) for code in range(K)])

    for colours in (_index_colours(palette), _bar_colours(palette), _band_colours(palette)):
        np.testing.assert_allclose(_by_code(colours), expected, rtol=0, atol=1e-12)


def test_a_shorter_colour_list_cycles_identically_everywhere():
    palette = ["#e41a1c", "#377eb8"]
    expected = np.array([to_rgba(palette[code % 2]) for code in range(K)])

    for colours in (_index_colours(palette), _bar_colours(palette), _band_colours(palette)):
        np.testing.assert_allclose(_by_code(colours), expected, rtol=0, atol=1e-12)
