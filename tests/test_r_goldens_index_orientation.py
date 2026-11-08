"""Index plots draw the first sequence at the bottom, as R does (TX-AUD-20261010-059).

R's index plot is barplot(horiz = TRUE) with one column per sequence
(plot.stslist.R:131,149-151): the first sequence is the lowest bar. traminex's
imshow kept Matplotlib's default origin and drew the first sequence at the top.
Fixture: ``scripts/r_goldens/make_goldens_index_orientation.R`` (barplot's bar
midpoints).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import traminex as tx  # noqa: E402

GOLDEN = pd.read_csv(Path(__file__).resolve().parent / "r_goldens" / "index_orientation.csv")
SEQS = tx.seqdef(
    np.array([list(s) for s in ["AABBC", "ABBBC", "BBAAC", "CCCBA"]], dtype=object),
    alphabet=["A", "B", "C"],
)


def _screen_heights(ax) -> list[float]:
    """Screen height of each drawn row's centre, in the order of the image's rows."""
    image = ax.get_images()[0]
    rows = np.asarray(image.get_array()).shape[0]
    _, _, bottom, top = image.get_extent()
    step = (top - bottom) / rows
    first = bottom if image.origin == "lower" else top
    sign = 1.0 if image.origin == "lower" else -1.0
    centres = [first + sign * (row + 0.5) * step for row in range(rows)]
    return [float(ax.transData.transform((0.0, y))[1]) for y in centres]


def test_r_draws_the_first_sequence_lowest() -> None:
    midpoints = GOLDEN.sort_values("sequence")["midpoint"].to_numpy(dtype=float)
    assert np.all(np.diff(midpoints) > 0)


def test_sequence_index_plot_draws_the_first_sequence_at_the_bottom() -> None:
    fig, ax = tx.sequence_index_plot(SEQS)
    fig.canvas.draw()
    heights = _screen_heights(ax)
    r_order = GOLDEN.sort_values("midpoint")["sequence"].tolist()
    assert [index + 1 for index in np.argsort(heights)] == r_order
    plt.close(fig)


def test_seqplot_type_i_draws_the_first_sequence_at_the_bottom() -> None:
    fig, ax = tx.seqplot(SEQS, type="i")
    fig.canvas.draw()
    assert np.all(np.diff(_screen_heights(ax)) > 0)
    plt.close(fig)


def test_grouped_index_plot_draws_its_first_row_at_the_bottom() -> None:
    fig, ax, order = tx.sequence_index_plot_grouped(SEQS, np.array(["b", "a", "b", "a"]))
    fig.canvas.draw()
    assert np.all(np.diff(_screen_heights(ax)) > 0)
    plt.close(fig)
