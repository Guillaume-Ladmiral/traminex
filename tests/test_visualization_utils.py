from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from traminex.visualization.utils import (
    _reserve_legend_space,
    compute_even_ticks,
    determine_layout,
    figure_to_buffer,
    legend_figure,
    save_plot,
)


def test_compute_even_ticks_handles_empty_short_and_thinned_axes():
    assert compute_even_ticks(0) == []
    assert compute_even_ticks(4, max_labels=10) == [0, 1, 2, 3]
    assert compute_even_ticks(3, max_labels=None) == [0, 1, 2]
    assert compute_even_ticks(11, max_labels=3) == [0, 5, 10]
    assert compute_even_ticks(11, max_labels=1) == [0, 10]


def test_determine_layout_validates_explicit_and_automatic_grids():
    assert determine_layout(0) == (0, 0)
    assert determine_layout(5, nrows=2, ncols=3) == (2, 3)
    assert determine_layout(5, layout="column", default_cols=2) == (3, 2)
    assert determine_layout(5, layout="grid") == (2, 3)

    with pytest.raises(ValueError, match="too small"):
        determine_layout(5, nrows=2, ncols=2)
    with pytest.raises(ValueError, match="Unsupported layout style"):
        determine_layout(5, layout="diagonal")


def test_legend_reservation_supports_both_positions_and_rejects_unknown_position():
    fig, axes = plt.subplots(1, 2)
    try:
        assert _reserve_legend_space(fig, axes, 0.2, "right") == (0.8, 0.0, 0.2, 1.0)
        assert _reserve_legend_space(fig, axes, 0.25, "bottom") == (0.0, 0.0, 1.0, 0.25)
        with pytest.raises(ValueError, match="legend position"):
            _reserve_legend_space(fig, axes, 0.2, "left")
    finally:
        plt.close(fig)


def test_legend_figure_validates_parallel_inputs():
    with pytest.raises(ValueError, match="same length"):
        legend_figure(["A"], ["red", "blue"])

    fig, ax = legend_figure(["A", "B"], ["red", "blue"], ncol=10)
    try:
        assert not ax.axison
        assert len(ax.get_legend().get_patches()) == 2
    finally:
        plt.close(fig)


def test_figure_serialization_closes_figure_and_rewinds_buffer():
    fig, _ax = plt.subplots()
    figure_number = fig.number
    buffer = figure_to_buffer(fig, fmt="png", dpi=40)

    assert buffer.tell() == 0
    assert buffer.read(8) == b"\x89PNG\r\n\x1a\n"
    assert not plt.fignum_exists(figure_number)


def test_save_plot_accepts_figure_tuple_and_creates_parent(tmp_path):
    fig, ax = plt.subplots()
    output = tmp_path / "nested" / "plot.png"
    try:
        save_plot((fig, ax), str(output), dpi=40)
        assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        with pytest.raises(ValueError, match="filepath cannot be empty"):
            save_plot(fig, "")
    finally:
        plt.close(fig)
