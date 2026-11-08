"""Consumed-component gaps, part 4: dissimilarity analysis, event sequences, statistics, plots.

The supervision repository's ``docs/TRAMINEX_COMPONENT_GAPS.md`` lists the lines of the definitions the
three sibling libraries use that no suite ran (section 6; re-measured on 2026-11-08 at ``2ed89d77``).
Part 4 gives each a test or a recorded reason. The tests here cover the dissimilarity, event,
statistics and plotting rows: refusals are checked with their message, and the rare paths against an
equivalent call that other tests already compare with R, or against TraMineR's own definition.
"""

from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")

import importlib  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from matplotlib.colors import ListedColormap, to_rgba  # noqa: E402

import traminex as tx  # noqa: E402
from traminex import (  # noqa: E402
    _rlinpack,  # noqa: E402
    diss_mfa,
    repseq_rf,
)
from traminex.visualization import palette, series, ticks  # noqa: E402
from traminex.visualization.plots import dispatch, distribution  # noqa: E402

# ``traminex.disstree`` is the function; the module holds the C-order helpers.
disstree_module = importlib.import_module("traminex.disstree")

FRAME = pd.DataFrame(
    {
        "a": ["A", "B", "A", "C", "B", "A"],
        "b": ["B", "B", "A", "C", "A", "A"],
        "c": ["C", "A", "B", "A", "A", "B"],
        "d": ["C", "A", "B", "A", "C", "B"],
    }
)


@pytest.fixture(scope="module")
def seqs():
    return tx.seqdef(FRAME)


@pytest.fixture(scope="module")
def diss(seqs):
    return np.asarray(tx.seqdist(seqs, method="HAM"), dtype=float)


# --- diss.py -------------------------------------------------------------------------------------


def test_dissassoc_refuses_a_group_of_the_wrong_length(diss):
    with pytest.raises(ValueError, match="group must have one value per observation"):
        tx.dissassoc(diss, group=[0, 1, 0], R=1)


@pytest.mark.parametrize(
    ("weights", "message"),
    [
        ([1.0, -1.0, 1.0, 1.0, 1.0, 1.0], "random-sampling weights must be nonnegative"),
        ([0.0] * 6, "random-sampling weights must have a positive sum"),
    ],
)
def test_dissassoc_random_sampling_refuses_unusable_weights(diss, weights, message):
    with pytest.raises(ValueError, match=message):
        tx.dissassoc(
            diss,
            group=[0, 0, 0, 1, 1, 1],
            weights=weights,
            R=2,
            weight_permutation="random-sampling",
        )


def test_disscenter_reads_true_as_first_and_refuses_other_medoid_choices(diss):
    as_bool = tx.disscenter(diss, medoids_index=True)
    as_first = tx.disscenter(diss, medoids_index="first")
    assert np.array_equal(np.asarray(as_bool), np.asarray(as_first))
    with pytest.raises(ValueError, match='medoids_index must be one of "first", "all" or None'):
        tx.disscenter(diss, medoids_index="middle")
    with pytest.raises(ValueError, match="group must have one value per observation"):
        tx.disscenter(diss, group=[0, 1])


# --- diss_mfa.py ---------------------------------------------------------------------------------


def test_dissmfacw_helpers_refuse_misshapen_inputs(diss):
    with pytest.raises(ValueError, match="weights must have length equal to nrow"):
        diss_mfa._gower_matrix_weighted(diss, False, np.ones(3))
    with pytest.raises(ValueError, match="design and weights must have conformable rows"):
        diss_mfa._weighted_qr(np.ones((6, 2)), np.ones(5))
    with pytest.raises(ValueError, match="predictors must be one-dimensional"):
        diss_mfa._as_predictor_series(np.ones((2, 2)))
    with pytest.raises(ValueError, match="factors must contain at least one predictor"):
        diss_mfa._build_design({})
    with pytest.raises(ValueError, match="all predictors must have the same length"):
        diss_mfa._build_design({"g": ["a", "b", "a"], "h": ["x", "y"]})
    with pytest.raises(ValueError, match="factors must contain at least one predictor"):
        diss_mfa._complete_predictors({}, 6)


def test_dissmfacw_refuses_zero_weights_and_a_non_positive_nperm(diss):
    factors = {"g": ["a", "a", "b", "b", "a", "b"]}
    with pytest.raises(ValueError, match="weights must have a positive sum"):
        tx.dissmfacw(diss, factors, weights=[0.0] * 6, nperm=2)
    with pytest.raises(ValueError, match="nperm must be a positive integer"):
        tx.dissmfacw(diss, factors, nperm=0)


def test_dissmfacw_gower_reads_the_matrix_as_an_already_built_gower_matrix(diss):
    # dissmfacw.R:47-52: with gower = TRUE the matrix is used as given, else gower_matrix() builds it.
    factors = {"g": ["a", "a", "b", "b", "a", "b"]}
    plain = tx.dissmfacw(diss, factors, nperm=9, seed=1)
    gower = diss_mfa._gower_matrix_weighted(diss, False, np.ones(6))
    prebuilt = tx.dissmfacw(gower, factors, nperm=9, seed=1, gower=True)
    numbers = plain.mfac.select_dtypes("number").columns
    assert len(numbers) > 0
    assert np.allclose(
        plain.mfac[numbers].to_numpy(dtype=float), prebuilt.mfac[numbers].to_numpy(dtype=float)
    )


# --- diss_p1.py ----------------------------------------------------------------------------------


def test_dissrep_accepts_condensed_and_list_input_as_the_square_matrix(diss):
    from scipy.spatial.distance import squareform

    square = tx.dissrep(diss, criterion="density")
    condensed = tx.dissrep(squareform(diss), criterion="density")
    as_list = tx.dissrep(diss.tolist(), criterion="density")
    for other in (condensed, as_list):
        assert np.array_equal(other["representatives"], square["representatives"])
        assert np.array_equal(other["scores"], square["scores"])
    with pytest.raises(ValueError, match="weights must have length equal to number of objects"):
        tx.dissrep(diss, weights=[1.0, 2.0])


# --- disstree.py ---------------------------------------------------------------------------------


def test_disstree_refuses_predictors_of_another_size(diss):
    with pytest.raises(
        ValueError, match="dissimilarity matrix and predictors must have the same size"
    ):
        tx.disstree(diss, pd.DataFrame({"g": ["a", "b"]}), R=1)


def test_weighted_inter_inertia_of_an_empty_group_is_zero(diss):
    empty = np.array([], dtype=int)
    group = np.array([0, 1, 2])
    weights = np.ones(6)
    assert disstree_module._c_weighted_inter_inertia(diss, empty, group, weights) == 0.0
    assert disstree_module._c_weighted_inter_inertia(diss, group, empty, weights) == 0.0


# --- events --------------------------------------------------------------------------------------


def test_seqetm_refuses_an_unknown_method_and_warns_about_commas(seqs):
    with pytest.raises(ValueError, match="method must be one of 'transition', 'period', 'state'"):
        tx.seqetm(seqs, method="duration")
    commas = tx.seqdef(pd.DataFrame({"a": ["x,1", "y"], "b": ["y", "x,1"]}))
    with pytest.warns(UserWarning, match="should not contain commas"):
        tx.seqetm(commas, method="transition")


def test_seqecreate_takes_a_transition_table_and_refuses_misshapen_ones(seqs):
    table = tx.seqetm(seqs, method="transition")
    from_table = tx.seqecreate(seqs, tevent=table)
    from_name = tx.seqecreate(seqs, tevent="transition")
    assert [str(e) for e in from_table] == [str(e) for e in from_name]
    with pytest.raises(
        ValueError, match="tevent matrix must match the state transition dimensions"
    ):
        tx.seqecreate(seqs, tevent=np.zeros((2, 2), dtype=object))


def test_seqecreate_refuses_missing_timestamps_events_and_text_times():
    with pytest.raises(ValueError, match="Could not find a timestamp argument"):
        tx.seqecreate(id=[1, 1], event=["a", "b"])
    with pytest.raises(ValueError, match="Could not find an event argument"):
        tx.seqecreate(id=[1, 1], timestamp=[0.0, 1.0])
    with pytest.raises(ValueError, match="timestamp values must be numeric"):
        tx.seqecreate(id=[1, 1], timestamp=["early", "late"], event=["a", "b"])


# --- repseq_rf.py --------------------------------------------------------------------------------


def test_oriented_leaves_non_finite_coordinates_alone():
    coordinates = np.array([np.nan, np.inf, -np.inf])
    assert repseq_rf._oriented(coordinates) is coordinates


@pytest.mark.parametrize(
    ("matrix", "message"),
    [
        (np.zeros((0, 0)), "diss must contain at least one case"),
        (np.array([[0.0, np.inf], [np.inf, 0.0]]), "diss must contain only finite values"),
    ],
)
def test_dissrf_refuses_empty_and_non_finite_dissimilarities(matrix, message):
    with pytest.raises(ValueError, match=message):
        tx.dissrf(matrix, k=1)


@pytest.mark.parametrize(
    ("weights", "message"),
    [
        ([1.0, 1.0], "weights must contain one value for each case"),
        ([1.0, -1.0, 1.0, 1.0, 1.0, 1.0], "weights must be finite and non-negative"),
        ([0.0] * 6, "weights must have a positive sum"),
    ],
)
def test_dissrf_refuses_unusable_weights(diss, weights, message):
    with pytest.raises(ValueError, match=message):
        tx.dissrf(diss, k=2, weights=weights)


@pytest.mark.parametrize(
    ("sortv", "message"),
    [
        ("median", 'the only scalar string allowed for sortv is "mds"'),
        ([1.0, 2.0], "length of sortv must equal the number of cases"),
        ([1.0, 2.0, np.nan, 4.0, 5.0, 6.0], "sortv must contain only finite values"),
        ([b"x"] * 6, "sortv must contain numeric values"),
    ],
)
def test_dissrf_refuses_unusable_sort_keys(diss, sortv, message):
    with pytest.raises(ValueError, match=message):
        tx.dissrf(diss, k=2, sortv=sortv)


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        ({"k": "2"}, TypeError, "k must be an integer or None"),
        ({"k": 0}, ValueError, "k must be between 1 and the number of cases"),
        ({"k": 2, "pow": 0}, ValueError, "pow must be finite and positive"),
    ],
)
def test_dissrf_refuses_bad_k_and_pow(diss, kwargs, error, message):
    with pytest.raises(error, match=message):
        tx.dissrf(diss, **kwargs)


def test_seqrf_refuses_a_dissimilarity_of_another_size(seqs, diss):
    with pytest.raises(ValueError, match="diss size must equal the number of sequences"):
        tx.seqrf(seqs, diss[:3, :3], k=1)


# --- stats ---------------------------------------------------------------------------------------


def test_seqtab_refuses_ragged_ranks(seqs):
    with pytest.raises(TypeError, match="idxs must be scalar 0 or a nonempty vector"):
        tx.seqtab(seqs, idxs=[[1], [1, 2]])


def test_seqient_refuses_a_base_that_is_not_a_number(seqs):
    with pytest.raises(ValueError, match="base must be finite, positive, and not equal to 1"):
        tx.seqient(seqs, base="ten")


def test_seqtrate_refuses_a_non_integer_lag(seqs):
    with pytest.raises(TypeError, match="lag must be an integer"):
        tx.seqtrate(seqs, lag=1.5)


# --- visualization -------------------------------------------------------------------------------


def test_ensure_cmap_turns_one_colour_into_a_colormap_and_passes_others_through(seqs):
    cmap = palette.ensure_cmap(seqs, "red")
    assert isinstance(cmap, ListedColormap)
    assert tuple(cmap(0)) == to_rgba("red")
    sentinel = 5
    assert palette.ensure_cmap(seqs, sentinel) is sentinel


def test_palette_renderer_ignores_cpal_for_entropy_and_refuses_cpal_with_cmap():
    def renderer(*args, **kwargs):
        return None

    assert dispatch._palette_renderer("Ht", renderer, {"cpal": ["red"]}) is renderer
    with pytest.raises(TypeError, match="cpal and cmap cannot be given together"):
        dispatch._palette_renderer("d", renderer, {"cpal": ["red"], "cmap": "viridis"})


def test_seqplot_refuses_an_unknown_legend_choice(seqs):
    with pytest.raises(ValueError, match="with_legend must be True, False or 'auto'"):
        tx.seqplot(seqs, type="d", with_legend="maybe")


def test_entropy_mode_of_the_distribution_plot_labels_its_axes(seqs):
    fig, ax = plt.subplots()
    try:
        distribution.state_distribution_plot(seqs, mode="Ht", ax=ax)
        assert ax.get_title() == "Entropy over time (Ht)"
        assert ax.get_ylabel() == "Entropy"
    finally:
        plt.close(fig)


def test_time_labels_without_axes_or_texts_are_returned_as_they_are():
    formatter = ticks._RTimeLabels({1.0: "t1", 2.0: ""})
    assert formatter.format_ticks([]) == []
    assert formatter.format_ticks([1.0]) == ["t1"]


# --- _rlinpack.py --------------------------------------------------------------------------------


def test_dqrdc2_pivot_refuses_a_vector():
    with pytest.raises(ValueError, match="matrix must be two-dimensional"):
        _rlinpack.dqrdc2_pivot(np.ones(3), 1e-7)


def test_no_warning_escapes_the_quiet_paths(seqs):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        palette.ensure_cmap(seqs, "blue")


# --- rare paths, measured on 2026-11-08 ------------------------------------------------------------


def test_dissrep_with_no_representative_wanted_refuses_the_empty_set(diss):
    with pytest.raises(ValueError, match="No representatives selected"):
        tx.dissrep(diss, nrep=0)


def test_dissrf_lets_one_heavy_case_span_several_groups(diss):
    result = tx.dissrf(diss, k=5, sortv=[0, 1, 2, 3, 4, 5], weights=[10, 1, 1, 1, 1, 1])
    gsize = result["sizes"]["gsize"]
    assert [round(float(np.sum(w)), 9) for w in result["weights_list"]] == [round(gsize, 9)] * 5
    assert [list(map(int, cases)) for cases in result["index_list"]][:3] == [[0], [0], [0]]
    assert 0 in result["index_list"][3]


def test_seqtab_with_zero_weights_reports_counts_without_proportions(seqs):
    table = tx.seqtab(seqs, weights=np.zeros(6))
    assert (table["count"] == 0).all()
    assert table["prop"].isna().all() and table["cumprop"].isna().all()


def test_state_statistics_of_an_object_without_sequences_or_positions():
    from traminex.sequences import STSList

    empty = STSList(data=np.zeros((0, 2), dtype=int), alphabet=["A", "B"])
    assert tx.seqistatd(empty).shape == (0, 2)
    assert tx.seqistatd(empty).dtype.kind == "i"
    assert tx.seqistatd(empty, prop=True).dtype.kind == "f"
    no_positions = STSList(data=np.zeros((2, 0), dtype=int), alphabet=["A"])
    assert np.array_equal(series.resolve_sort_order(no_positions, sortv="from.start"), np.arange(2))
    assert np.isnan(tx.seqST(no_positions, norm=True)).all()


def test_entropy_plot_takes_the_groups_itself(seqs):
    result = tx.seqplot(seqs, type="Ht", group=[0, 0, 0, 1, 1, 1])
    try:
        figure = result[0]
        assert len(figure.axes) == 1
        assert len(figure.axes[0].lines) == 2
    finally:
        plt.close("all")


def test_time_labels_skip_empty_texts_when_an_axis_is_present():
    fig, ax = plt.subplots()
    try:
        ax.plot([1, 2, 3], [1, 2, 3])
        formatter = ticks._RTimeLabels({1.0: "t1", 2.0: "", 3.0: "t3"})
        ax.xaxis.set_major_formatter(formatter)
        texts = formatter.format_ticks([1.0, 2.0, 3.0])
        assert texts[0] == "t1" and texts[1] == "" and texts[2] == "t3"
    finally:
        plt.close(fig)


EVENTS = dict(
    id=[1, 1, 1, 3, 3, 3],
    timestamp=[0.0, 1.5, 3.25, 0.5, 1.0, 4.75],
    event=["a", "b", "c", "a", "c", "b"],
)


def test_an_event_sequence_with_only_its_end_event_counts_as_empty():
    eseq = tx.seqecreate(
        id=[1, 1, 1, 2, 3, 3, 3],
        timestamp=[0.0, 1.5, 3.25, 2.0, 0.5, 1.0, 4.75],
        event=["a", "b", "c", "END", "a", "c", "b"],
        end_event="END",
    )
    assert str(eseq[1]) == ""
    subsequences = tx.seqefsub(eseq, pmin_support=0.01)
    support = dict(zip((str(s) for s in subsequences.subseq), subsequences.data["Support"]))
    assert support["(a)"] == pytest.approx(2 / 3)
    ages = tx.seqeapplysub(subsequences, method="age", constraint=tx.seqeconstraint(age_max=1.0))
    # Sequence::count returns 0 for an empty sequence; a non-empty one without the subsequence gives -1.
    assert (ages.iloc[1] == 0).all()
    assert ages.loc[ages.index[0], "(b)"] == -1


@pytest.mark.parametrize("count_method", [4, 5])
def test_age_max_drops_subsequences_that_start_later(count_method):
    eseq = tx.seqecreate(**EVENTS)
    limited = tx.seqefsub(
        eseq,
        pmin_support=0.01,
        constraint=tx.seqeconstraint(count_method=count_method, age_max=1.0),
    )
    unlimited = tx.seqefsub(
        eseq, pmin_support=0.01, constraint=tx.seqeconstraint(count_method=count_method)
    )
    limited_names = {str(s) for s in limited.subseq}
    unlimited_names = {str(s) for s in unlimited.subseq}
    assert {"(b)", "(b)-(c)"} <= unlimited_names
    assert not {"(b)", "(b)-(c)"} & limited_names
    assert limited_names < unlimited_names
