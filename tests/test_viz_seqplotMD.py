"""RP7a: seqplotMD (multidomain) structure test.

One panel per domain; each ``type="i"`` panel draws that domain's raw state
matrix (an index plot). Asserted by checking each panel's imshow array == the
domain's data.
"""

from __future__ import annotations

from copy import deepcopy

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from corestats_fixtures import fixture_a
from matplotlib.colors import to_rgba

from traminex.format import seqdef
from traminex.visualization.plots.multidomain import seqplotMD


def _dom2():
    seqs = np.array(
        [
            ["X", "X", "Y", "Y", "Z", "Z"],
            ["Y", "Y", "Y", "Z", "Z", "X"],
            ["Z", "Z", "Z", "Z", "Z", "Z"],
            ["X", "Y", "Z", "X", "Y", "Z"],
        ],
        dtype=object,
    )
    return seqdef(seqs, alphabet=["X", "Y", "Z"])


def test_seqplotMD_one_index_panel_per_domain():
    dom1 = fixture_a()
    dom2 = _dom2()
    fig, axes = seqplotMD([dom1, dom2])
    assert len(axes) == 2
    assert np.array_equal(np.asarray(axes[0].images[0].get_array()), dom1.data)
    assert np.array_equal(np.asarray(axes[1].images[0].get_array()), dom2.data)
    assert axes[0].get_title() == "Dom1" and axes[1].get_title() == "Dom2"
    plt.close(fig)


def test_seqplotMD_requires_two_domains():
    with pytest.raises(ValueError):
        seqplotMD([fixture_a()])


def test_seqplotMD_applies_each_index_palette_and_all_domain_names():
    palettes = [
        ["#110000", "#220000", "#330000", "#440000"],
        ["#001100", "#002200", "#003300"],
    ]
    fig, axes = seqplotMD(
        [fixture_a(), _dom2()],
        cpal=palettes,
        dnames=["Work", "Family"],
        type="i",
        legend=True,
    )

    assert len(axes) == 2
    assert [ax.get_title() for ax in axes] == ["Work", "Family"]
    for ax, palette in zip(axes, palettes):
        assert list(ax.images[0].get_cmap().colors) == [to_rgba(color) for color in palette]
        assert ax.get_legend() is not None
    plt.close(fig)


def test_seqplotMD_applies_each_distribution_palette():
    palettes = [
        ["#100000", "#200000", "#300000", "#400000"],
        ["#001000", "#002000", "#003000"],
    ]
    channels = [fixture_a(), _dom2()]
    fig, axes = seqplotMD(channels, cpal=palettes, type="d", include_missing=True)

    assert len(axes) == 2
    for ax, channel, palette in zip(axes, channels, palettes):
        sequence_length = channel.data.shape[1]
        actual = [ax.patches[k * sequence_length].get_facecolor() for k in range(len(palette))]
        assert actual == [to_rgba(color) for color in palette]
    plt.close(fig)


@pytest.mark.parametrize("cpal", ["red", [["red"]], [["red"], ["blue"], ["green"]]])
def test_seqplotMD_rejects_invalid_outer_palette_shape(cpal):
    with pytest.raises(
        ValueError,
        match="cpal must be a collection with exactly one palette per domain",
    ):
        seqplotMD([fixture_a(), _dom2()], cpal=cpal)


@pytest.mark.parametrize("dnames", ["AB", ["only one"], ["one", "two", "three"]])
def test_seqplotMD_rejects_invalid_domain_name_cardinality(dnames):
    with pytest.raises(ValueError, match="dnames must contain exactly one name per domain"):
        seqplotMD([fixture_a(), _dom2()], dnames=dnames)


def test_seqplotMD_repeated_calls_do_not_mutate_stored_or_requested_palettes():
    channels = [fixture_a(), _dom2()]
    stored_before = [deepcopy(channel.cpal) for channel in channels]
    first = [["red", "green", "blue", "black"], ["cyan", "magenta", "yellow"]]
    second = [["white", "gray", "pink", "brown"], ["navy", "olive", "teal"]]
    requested_before = deepcopy([first, second])

    fig1, axes1 = seqplotMD(channels, cpal=first)
    fig2, axes2 = seqplotMD(channels, cpal=second)

    assert [channel.cpal for channel in channels] == stored_before
    assert [first, second] == requested_before
    assert axes1[0].images[0].get_cmap().colors[0] == to_rgba("red")
    assert axes2[0].images[0].get_cmap().colors[0] == to_rgba("white")
    plt.close(fig1)
    plt.close(fig2)


# With sortv "from.start"/"from.end", R orders the multidomain sequences once
# and draws every domain panel in that order (TX-AUD-20261010-008). R's drawn
# order: scripts/r_goldens/make_goldens_seqplotmd_sort_order.R.
def _sort_order_domains():
    domain1 = seqdef(
        np.array([list("AAB"), list("BAA"), list("ABB"), list("BBA"), list("AAA")], dtype=object)
    )
    domain2 = seqdef(
        np.array([list("yxx"), list("xxy"), list("xyy"), list("yyx"), list("xxy")], dtype=object)
    )
    return domain1, domain2


def _r_drawn_order(case: str, panel: int) -> list[int]:
    import csv

    from r_goldens_utils import R_GOLDENS_DIR

    with (R_GOLDENS_DIR / "seqplotmd_sort_order.csv").open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle)
            if row["case"] == case and int(row["panel"]) == panel
        ]
    return [int(row["row"]) - 1 for row in sorted(rows, key=lambda row: int(row["position"]))]


@pytest.mark.parametrize("sortv", ["from.start", "from.end"])
def test_seqplotMD_string_sortv_orders_every_panel_by_the_multidomain_sequences(sortv):
    domains = _sort_order_domains()
    fig, axes = seqplotMD(list(domains), type="i", sortv=sortv)
    for panel, (ax, domain) in enumerate(zip(axes, domains), start=1):
        # The index renderer draws its display order from the top row down.
        drawn = np.asarray(ax.images[0].get_array())
        expected = np.asarray(domain.data)[_r_drawn_order(sortv, panel)]
        np.testing.assert_array_equal(drawn, expected)
    plt.close(fig)


def test_seqplotMD_numeric_sortv_keeps_panels_aligned():
    domains = _sort_order_domains()
    order_key = np.array([3.0, 1.0, 4.0, 0.0, 2.0])
    fig, axes = seqplotMD(list(domains), type="i", sortv=order_key)
    for ax, domain in zip(axes, domains):
        drawn = np.asarray(ax.images[0].get_array())
        np.testing.assert_array_equal(drawn, np.asarray(domain.data)[[3, 1, 4, 0, 2]])
    plt.close(fig)
