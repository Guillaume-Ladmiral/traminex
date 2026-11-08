"""R4 sampling contract for TraMineR ``plot.stslist`` index plots."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import QuadMesh
from r_goldens_utils import R_GOLDENS_DIR

from traminex.sequences import STSList
from traminex.visualization import series
from traminex.visualization.plots.index import sequence_index_plot

_DEFAULT_WEIGHTS = object()


def _fixture(*, weights: np.ndarray | None | object = _DEFAULT_WEIGHTS) -> STSList:
    if weights is _DEFAULT_WEIGHTS:
        weights = np.array([1.0, 4.0, 0.0, 2.0, 3.0, 5.0])
    return STSList(
        data=np.array(
            [
                [1, 1, 0, 0],
                [0, 2, 2, 1],
                [0, 0, 1, 1],
                [2, 2, 2, 2],
                [0, 4, 4, 4],
                [1, 0, 0, 4],
            ],
            dtype=int,
        ),
        alphabet=["Alpha", "Beta", "Gamma"],
        labels=["Alpha", "Beta", "Gamma"],
        cpal=["#111111", "#555555", "#999999"],
        xtlab=["t1", "t2", "t3", "t4"],
        xtstep=2,
        tick_last=True,
        nr_code=3,
        void_code=4,
        weights=weights,
        start=7,
        ids=[f"row-{index}" for index in range(1, 7)],
    )


def _golden(case: str) -> pd.DataFrame:
    records = pd.read_csv(R_GOLDENS_DIR / "plot_stslist_sampling_r4_selections.csv")
    return records.loc[records["case"] == case].reset_index(drop=True)


def _condition(case: str) -> str:
    records = pd.read_csv(
        R_GOLDENS_DIR / "plot_stslist_sampling_r4_conditions.csv",
        dtype=str,
        keep_default_na=False,
    )
    return records.loc[records["case"] == case, "message"].item()


@pytest.mark.parametrize(
    ("case", "sortv", "sampling"),
    [
        ("input_systematic_3", None, 3),
        ("from_start_systematic_3", "from.start", 3),
        ("from_end_systematic_4", "from.end", 4),
        ("numeric_ties_systematic_3", np.array([2, 0, 2, 1, 0, 1]), 3),
        ("input_one", None, 1),
        ("input_equal", None, 6),
        ("input_above", None, 9),
        ("input_no_sampling", None, None),
    ],
)
def test_systematic_selection_matches_r_golden(case, sortv, sampling):
    x = _fixture()
    selected = series.select_index_rows(x, sortv=sortv, sampling=sampling)
    expected = _golden(case)

    np.testing.assert_array_equal(selected.source_rows, expected["source_index_zero"])
    np.testing.assert_array_equal(selected.display_rows, x.data[selected.source_rows])
    np.testing.assert_allclose(selected.weights, expected["selected_weight"], rtol=0, atol=0)
    assert selected.labels == tuple(expected["source_id"])


def test_factor_like_sort_uses_category_codes_like_r():
    sortv = pd.Categorical(["z", "x", "z", "y", "x", "y"], categories=["y", "z", "x"])
    selected = series.select_index_rows(_fixture(), sortv=sortv, sampling=4)
    expected = _golden("factor_ties_systematic_4")
    np.testing.assert_array_equal(selected.source_rows, expected["source_index_zero"])


def test_factor_like_sort_puts_missing_values_last_like_r():
    sortv = pd.Categorical(
        ["b", None, "a", "b", None, "a"],
        categories=["a", "b"],
        ordered=True,
    )
    selected = series.select_index_rows(_fixture(), sortv=sortv)
    expected = _golden("factor_missing_all")
    np.testing.assert_array_equal(selected.source_rows, expected["source_index_zero"])


def test_random_sampling_is_reproducible_sorted_and_weighted():
    first = series.select_index_rows(
        _fixture(), sampling=3, sample_meth="random", random_state=1729
    )
    second = series.select_index_rows(
        _fixture(), sampling=3, sample_meth="random", random_state=1729
    )
    np.testing.assert_array_equal(first.source_rows, second.source_rows)
    assert np.all(first.source_rows[:-1] < first.source_rows[1:])
    assert 2 not in first.source_rows  # zero-weight candidate cannot be drawn


def test_random_sampling_accepts_generator_and_does_not_touch_global_rng():
    np.random.seed(884)
    expected_global = np.random.random(4)
    np.random.seed(884)
    generator = np.random.default_rng(71)
    selected = series.select_index_rows(
        _fixture(), sampling=2, sample_meth="random", random_state=generator
    )
    actual_global = np.random.random(4)

    assert selected.source_rows.shape == (2,)
    np.testing.assert_array_equal(actual_global, expected_global)


def test_unweighted_random_sampling_returns_unit_weights():
    selected = series.select_index_rows(
        _fixture(), sampling=4, sample_meth="random", weighted=False, random_state=9
    )
    np.testing.assert_array_equal(selected.weights, np.ones(4))


@pytest.mark.parametrize("weights", [None, np.ones(6)])
def test_absent_and_unit_weights_select_as_unit_weights(weights):
    selected = series.select_index_rows(_fixture(weights=weights), sampling=3)
    np.testing.assert_array_equal(selected.weights, np.ones(3))


@pytest.mark.parametrize("sampling", [0, -1, 1.5, True])
def test_invalid_sampling_size_is_actionable(sampling):
    if sampling == 0:
        assert "differing number of rows" in _condition("zero_sampling")
    elif sampling == -1:
        assert "non-negative" in _condition("negative_sampling")
    with pytest.raises((TypeError, ValueError), match="sampling"):
        series.select_index_rows(_fixture(), sampling=sampling)


def test_invalid_sampling_method_is_actionable_even_without_subsampling():
    oracle_message = _condition("unknown_method")
    assert "systematic" in oracle_message and "random" in oracle_message
    with pytest.raises(ValueError, match="sample_meth.*systematic.*random"):
        series.select_index_rows(_fixture(), sample_meth="bogus")


@pytest.mark.parametrize(
    "weights",
    [
        np.array([1.0, 2.0]),
        np.array([1.0, -1.0, 1.0, 1.0, 1.0, 1.0]),
        np.array([1.0, np.nan, 1.0, 1.0, 1.0, 1.0]),
        np.zeros(6),
    ],
)
def test_invalid_weight_vectors_are_rejected_for_weighted_random_sampling(weights):
    with pytest.raises(ValueError, match="weights"):
        series.select_index_rows(
            _fixture(weights=weights), sampling=3, sample_meth="random", random_state=1
        )


def test_random_sampling_requires_enough_positive_weight_candidates():
    weights = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="positive weights"):
        series.select_index_rows(
            _fixture(weights=weights), sampling=3, sample_meth="random", random_state=1
        )


def test_selection_and_plot_do_not_mutate_the_caller():
    x = _fixture()
    original_data = x.data.copy()
    original_weights = x.weights.copy()
    original_ids = list(x.ids)

    fig, _ = sequence_index_plot(x, sortv="from.end", sampling=3)
    plt.close(fig)

    np.testing.assert_array_equal(x.data, original_data)
    np.testing.assert_array_equal(x.weights, original_weights)
    assert x.ids == original_ids


def _drawn_rows(ax) -> np.ndarray:
    """The state codes drawn: the image's, or the mesh's when weighted rows have their own heights."""
    if ax.images:
        return np.asarray(ax.images[0].get_array())
    mesh = next(c for c in ax.collections if isinstance(c, QuadMesh))
    rows, cols = (np.asarray(mesh.get_coordinates()).shape[i] - 1 for i in (0, 1))
    return np.asarray(mesh.get_array()).reshape(rows, cols)


def test_plot_exposes_sampled_rows_labels_and_post_sampling_highlight():
    x = _fixture()
    fig, ax = sequence_index_plot(x, sortv="from.start", sampling=3, highlight=[1])
    expected = _golden("from_start_systematic_3")
    source = expected["source_index_zero"].to_numpy()

    np.testing.assert_array_equal(_drawn_rows(ax), x.data[source])
    assert [tick.get_text() for tick in ax.get_yticklabels()] == list(expected["source_id"])
    # The highlighted second row lies between the cumulated weights of the rows drawn below and through it
    # (R's barplot(width = weights), TX-AUD-20261019-048).
    edges = np.concatenate(([0.0], np.cumsum(np.asarray(x.weights, dtype=float)[source])))
    assert [tuple(line.get_ydata()) for line in ax.lines] == [
        (edges[1], edges[1]),
        (edges[2], edges[2]),
    ]
    plt.close(fig)


def test_plot_accepts_supplied_axes_and_legacy_no_sampling_call():
    x = _fixture(weights=None)
    fig, supplied = plt.subplots()
    returned_fig, returned_ax = sequence_index_plot(x, ax=supplied)

    assert returned_fig is fig
    assert returned_ax is supplied
    np.testing.assert_array_equal(np.asarray(supplied.images[0].get_array()), x.data)
    plt.close(fig)


def test_empty_and_zero_width_inputs_are_safe():
    empty = STSList(data=np.empty((0, 4), dtype=int), alphabet=["A"])
    zero_width = STSList(data=np.empty((2, 0), dtype=int), alphabet=["A"], ids=["x", "y"])

    empty_selection = series.select_index_rows(empty, sampling=None)
    assert empty_selection.source_rows.size == 0
    empty_fig, empty_ax = sequence_index_plot(empty)
    zero_fig, zero_ax = sequence_index_plot(zero_width, sampling=1)

    assert not empty_ax.images
    assert not zero_ax.images
    plt.close(empty_fig)
    plt.close(zero_fig)


@pytest.mark.parametrize(
    "data",
    [
        np.array([[0, 3, 1, 4], [1, 1, 4, 4]], dtype=int),
        np.full((2, 4), 4, dtype=int),
        np.array([[0, 1, 1, 4]], dtype=int),
    ],
)
def test_ragged_internal_missing_all_void_and_single_rows_render(data):
    x = STSList(
        data=data,
        alphabet=["A", "B"],
        nr_code=3,
        void_code=4,
        ids=[f"special-{index}" for index in range(data.shape[0])],
    )
    fig, ax = sequence_index_plot(x, sampling=1)
    assert np.asarray(ax.images[0].get_array()).shape == (1, data.shape[1])
    plt.close(fig)


def test_oracle_metadata_records_the_two_deliberate_divergences():
    metadata = pd.read_csv(R_GOLDENS_DIR / "plot_stslist_sampling_r4_metadata.csv")
    divergence_keys = set(metadata.loc[metadata["section"] == "deliberate_divergence", "key"])
    assert divergence_keys == {"python_sampling_none", "python_random_stream"}
    assert metadata.loc[metadata["key"] == "sort_before_sampling", "value"].item() == "TRUE"
