"""Small paths of what the other libraries use that no test here ran before.

From the supervision repository's measurement of 2026-11-02 (its
``docs/TRAMINEX_COMPONENT_GAPS.md``). py_seqplot's tests were the only ones to run
``resolve_sort_order``'s refusals and ``_exact_near_ties``' early return, and no test here
called ``apply_time_axis`` directly, though pyweightedcluster does and the function is
outside the public API. R's values come from
``scripts/r_goldens/make_goldens_consumed_items.R``; the ticks are the rows
``scripts/r_goldens/make_goldens_plot_items.R`` already records for R's own plots.
"""

from __future__ import annotations

import csv
import importlib

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR  # noqa: E402

from traminex import dissrep, seqdef  # noqa: E402
from traminex.cluster import seqrep  # noqa: E402
from traminex.visualization.series import resolve_sort_order  # noqa: E402
from traminex.visualization.ticks import apply_time_axis  # noqa: E402

diss_p1 = importlib.import_module("traminex.diss_p1")


def _read(name: str) -> dict[tuple[str, str], str]:
    with (R_GOLDENS_DIR / name).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(line for line in handle if not line.startswith("#")))
    return {(case, key): value for case, key, value in rows[1:]}


ITEMS = _read("consumed_items.csv")
PLOT_ITEMS = _read("plot_items.csv")
THREE = [list("ABB"), list("BBA"), list("AAA")]


def _numbers(text: str) -> np.ndarray:
    return np.array(text.split(" "), dtype=float)


def _assert_statistics(frame, case: str) -> None:
    # R's Statistics, written row by row: na, na(%), nb, nb(%), SD, MD, DC, V, Q; r1, then Total.
    expected = _numbers(ITEMS[(case, "statistics")]).reshape(2, 9)
    np.testing.assert_allclose(frame.to_numpy(dtype=float), expected, rtol=0, atol=R_GOLDEN_ATOL)


@pytest.mark.parametrize("criterion", ["density", "freq", "dist"])
def test_dissrep_of_a_single_object_matches_r(criterion):
    result = dissrep(np.zeros((1, 1)), criterion=criterion)
    case = f"dissrep_single_{criterion}"
    assert result["representatives"].tolist() == [
        int(value) - 1 for value in _numbers(ITEMS[(case, "representatives")])
    ]
    np.testing.assert_allclose(
        result["scores"], _numbers(ITEMS[(case, "scores")]), rtol=0, atol=R_GOLDEN_ATOL
    )
    _assert_statistics(result["statistics"], case)


def test_seqrep_of_a_single_sequence_matches_r():
    result = seqrep(seqdef([list("ABB")], alphabet=["A", "B"]), diss=np.zeros((1, 1)))
    _assert_statistics(result["statistics"], "seqrep_single")


def test_near_tie_resumming_leaves_a_single_score_alone():
    score = np.array([2.5])
    assert diss_p1._exact_near_ties(score, np.ones((1, 1)), np.ones(1)) is score


def test_an_unknown_sortv_string_is_refused_as_in_r():
    # R's message is the one for a wrong length, because it reads any other string as a
    # vector of length 1. With one sequence R would accept it, which this test leaves unpinned.
    assert "sortv must contain one value" in ITEMS[("sortv_unknown_string", "refusal")]
    x = seqdef(THREE, alphabet=["A", "B"])
    with pytest.raises(ValueError, match="sortv string must be 'from.start' or 'from.end'"):
        resolve_sort_order(x, "not-a-method")


@pytest.mark.parametrize("sortv", [[1, 2], np.array([1, 2]), np.zeros((3, 2))])
def test_a_sortv_of_the_wrong_shape_is_refused_as_in_r(sortv):
    assert "sortv must contain one value" in ITEMS[("sortv_wrong_length", "refusal")]
    x = seqdef(THREE, alphabet=["A", "B"])
    with pytest.raises(ValueError, match="sortv must contain one value for each sequence"):
        resolve_sort_order(x, sortv)


@pytest.mark.parametrize("weighted", [True, False])
def test_an_integer_sortv_array_orders_as_r(weighted):
    # The order does not depend on weighted: R's weights set bar widths, not the order.
    x = seqdef(THREE, alphabet=["A", "B"])
    order = resolve_sort_order(x, np.array([3, 1, 2]), weighted=weighted)
    assert order.tolist() == [
        int(value) - 1 for value in _numbers(ITEMS[("sortv_integers", "order")])
    ]


TICK_CASES = ["ticks_step6_72", "ticks_default_72", "ticks_step5_last_72", "ticks_step3_last_5"]


def _tick_inputs(case: str):
    length = int(PLOT_ITEMS[(case, "length")])
    cells = np.array(
        [[("A", "B")[(row + t) % 2] for t in range(length)] for row in range(3)], dtype=object
    )
    x = seqdef(
        cells,
        alphabet=["A", "B"],
        cnames=[f"c{t}" for t in range(1, length + 1)],
        xtstep=int(PLOT_ITEMS[(case, "xtstep")]),
        tick_last=PLOT_ITEMS[(case, "tick_last")] == "TRUE",
    )
    # R draws position t (1-based) at t - 0.5; apply_time_axis places position t (0-based) at t.
    positions = _numbers(PLOT_ITEMS[(case, "positions")]) - 0.5
    labels = PLOT_ITEMS[(case, "labels")].split(" ;; ")
    return x, length, positions, labels


@pytest.fixture
def axes():
    fig, ax = plt.subplots()
    yield ax
    plt.close(fig)


@pytest.mark.parametrize("case", TICK_CASES)
def test_apply_time_axis_places_r_s_ticks_and_labels(case, axes):
    x, length, positions, labels = _tick_inputs(case)
    apply_time_axis(axes, x, length)
    np.testing.assert_allclose(axes.get_xticks(), positions, rtol=0, atol=1e-12)
    formatter = axes.xaxis.get_major_formatter()
    assert [formatter(value) for value in axes.get_xticks()] == labels


@pytest.mark.parametrize("case", TICK_CASES)
def test_apply_time_axis_one_based_shifts_the_ticks_only(case, axes):
    x, length, positions, labels = _tick_inputs(case)
    apply_time_axis(axes, x, length, one_based=True)
    np.testing.assert_allclose(axes.get_xticks(), positions + 1.0, rtol=0, atol=1e-12)
    formatter = axes.xaxis.get_major_formatter()
    assert [formatter(value) for value in axes.get_xticks()] == labels


def test_apply_time_axis_max_labels_no_longer_limits_the_ticks(axes):
    x, length, positions, _ = _tick_inputs("ticks_default_72")
    apply_time_axis(axes, x, length, max_labels=2)
    assert axes.get_xticks().size == positions.size == 72


def test_apply_time_axis_styles_the_labels(axes):
    x, length, _, _ = _tick_inputs("ticks_step3_last_5")
    apply_time_axis(axes, x, length, label_color="red", label_fontsize=7.0)
    axes.figure.canvas.draw()
    labels = axes.xaxis.get_ticklabels()
    assert all(to_rgba(label.get_color()) == to_rgba("red") for label in labels)
    assert all(label.get_fontsize() == 7.0 for label in labels)


def test_apply_time_axis_leaves_the_default_style_alone(axes):
    x, length, _, _ = _tick_inputs("ticks_step3_last_5")
    before = axes.xaxis.get_tick_params()
    apply_time_axis(axes, x, length, label_color=None, label_fontsize=None)
    assert axes.xaxis.get_tick_params() == before


def test_apply_time_axis_numbers_from_start_without_column_labels(axes):
    # Without usable labels the positions are numbered from x.start.
    class Sequences:
        xtstep = 2
        start = 5
        tick_last = False
        xtlab = None

    apply_time_axis(axes, Sequences(), 5)
    np.testing.assert_array_equal(axes.get_xticks(), [0.0, 2.0, 4.0])
    formatter = axes.xaxis.get_major_formatter()
    assert [formatter(value) for value in axes.get_xticks()] == ["5", "7", "9"]
