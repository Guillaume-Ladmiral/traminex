"""R4 reference rendering for TraMineRextras ``plot.SAMM``."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDENS_DIR

import traminex
from traminex.samm import SAMMResult
from traminex.sequences import STSList


def _ragged_fixture() -> STSList:
    return STSList(
        data=np.array(
            [
                [0, 0, 1, 1, 2, 2, 2],
                [0, 1, 1, 2, 4, 4, 4],
                [2, 2, 2, 2, 2, 2, 2],
                [0, 0, 0, 3, 1, 1, 4],
            ],
            dtype=int,
        ),
        alphabet=["Alpha", "Beta", "Gamma"],
        labels=["Alpha", "Beta", "Gamma"],
        cpal=["#111111", "#222222", "#333333"],
        xtlab=["t1", "t2", "t3", "t4", "t5", "t6", "t7"],
        xtstep=2,
        tick_last=True,
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
        weights=np.array([1.0, 2.0, 3.0, 4.0]),
        start=5,
        ids=["one", "two", "three", "four"],
        missing_color="#999999",
    )


def _repeated_fixture() -> STSList:
    return STSList(
        data=np.array([[0, 1, 0, 0, 2, 2, 0], [0, 0, 1, 1, 0, 2, 2]], dtype=int),
        alphabet=["State A", "State B", "State C"],
        labels=["State A", "State B", "State C"],
        cpal=["#AA0000", "#00AA00", "#0000AA"],
        nr_code=3,
        void_code=4,
        start=4,
    )


def _constant_fixture() -> STSList:
    return STSList(
        data=np.array([[0, 0, 0, 0, 0], [1, 1, 1, 1, 1]], dtype=int),
        alphabet=["A", "B"],
        nr_code=2,
        void_code=3,
    )


def _plot(*args, **kwargs):
    return getattr(traminex, "plot_SAMM")(*args, **kwargs)


def _axes(result):
    return list(np.asarray(result[1], dtype=object).reshape(-1))


def _golden_selection(case: str) -> pd.DataFrame:
    records = pd.read_csv(R_GOLDENS_DIR / "plot_samm_r4_selection.csv")
    return records.loc[records["case"] == case].reset_index(drop=True)


def _golden_metadata(case: str, key: str) -> str:
    records = pd.read_csv(
        R_GOLDENS_DIR / "plot_samm_r4_metadata.csv", dtype=str, keep_default_na=False
    )
    return records.loc[(records["case"] == case) & (records["key"] == key), "value"].item()


def _golden_condition(case: str) -> str:
    records = pd.read_csv(
        R_GOLDENS_DIR / "plot_samm_r4_conditions.csv", dtype=str, keep_default_na=False
    )
    return records.loc[records["case"] == case, "value"].item()


_SUPPORTED_PLOT_TYPES = tuple(_golden_condition("supported_dispatch_types").split("\x1e"))


def _selected_distance(samm: SAMMResult) -> np.ndarray:
    selected = samm.sequences.data[samm.frame["transition"].to_numpy(dtype=bool)]
    return np.count_nonzero(selected[:, None, :] != selected[None, :, :], axis=2).astype(float)


def test_plot_samm_index_selection_rows_groups_and_titles_match_r_golden():
    samm = traminex.seqsamm(_ragged_fixture(), 3)
    fig, axes = _plot(samm, type="i")
    flat = _axes((fig, axes))
    golden = _golden_selection("ragged")
    group_order = _golden_metadata("ragged", "group_order").split("\x1e")

    assert [axis.get_title() for axis in flat] == [
        f"Transition out of {origin}" for origin in group_order
    ]
    for axis, origin in zip(flat, group_order):
        selected_indexes = golden.loc[golden["origin"] == origin, "selected_index_zero"].unique()
        expected = samm.sequences.data[samm.frame["transition"].to_numpy(dtype=bool)][
            selected_indexes
        ]
        np.testing.assert_array_equal(np.asarray(axis.images[0].get_array()), expected)
    plt.close(fig)


def test_plot_samm_selected_code_rows_and_source_positions_match_r_golden():
    samm = traminex.seqsamm(_ragged_fixture(), 3)
    golden = _golden_selection("ragged")
    transition_positions = np.flatnonzero(samm.frame["transition"].to_numpy(dtype=bool))
    selected = samm.sequences.data[transition_positions]

    np.testing.assert_array_equal(
        transition_positions,
        golden.groupby("selected_index_zero", sort=True)["source_frame_index_zero"].first(),
    )
    for selected_index, records in golden.groupby("selected_index_zero", sort=True):
        np.testing.assert_array_equal(selected[selected_index], records["state_code_zero"])


@pytest.mark.parametrize("plot_type", _SUPPORTED_PLOT_TYPES[:-2])
def test_plot_samm_supports_each_established_non_distance_dispatch_type(plot_type):
    samm = traminex.seqsamm(_repeated_fixture(), 2)
    fig, axes = _plot(samm, type=plot_type)
    assert len(_axes((fig, axes))) == 2
    plt.close(fig)


def test_plot_samm_supports_representative_route_with_selected_distance():
    samm = traminex.seqsamm(_repeated_fixture(), 2)
    fig, axes = _plot(samm, type="r", diss=_selected_distance(samm), nrep=1)
    assert len(_axes((fig, axes))) == 2
    assert all(axis.patches for axis in _axes((fig, axes)))
    plt.close(fig)


@pytest.mark.parametrize("which_plot", ["medoids", "diss.to.med", "both"])
def test_plot_samm_supports_relative_frequency_route_and_panel_modes(which_plot):
    samm = traminex.seqsamm(_repeated_fixture(), 2)
    fig, axes = _plot(
        samm,
        type="rf",
        diss=_selected_distance(samm),
        k=1,
        which_plot=which_plot,
    )
    expected_axes = 4 if which_plot == "both" else 2
    assert len(_axes((fig, axes))) == expected_axes
    plt.close(fig)


@pytest.mark.parametrize("plot_type", _SUPPORTED_PLOT_TYPES[-2:])
def test_plot_samm_distance_routes_require_selected_transition_distance(plot_type):
    oracle_message = _golden_condition(f"{plot_type}_requires_diss")
    assert any(token in oracle_message.lower() for token in ("diss", "distance"))
    with pytest.raises(ValueError, match="diss.*transition"):
        _plot(traminex.seqsamm(_repeated_fixture(), 2), type=plot_type)


def test_plot_samm_rejects_wrong_distance_size():
    samm = traminex.seqsamm(_repeated_fixture(), 2)
    with pytest.raises(ValueError, match="diss.*transition"):
        _plot(samm, type="r", diss=np.zeros((2, 2)))


def test_plot_samm_supplied_axes_and_style_forwarding():
    samm = traminex.seqsamm(_ragged_fixture(), 3)
    fig, supplied = plt.subplots(1, 3)
    returned_fig, axes = _plot(
        samm,
        type="i",
        axes=supplied,
        sampling=1,
        legend=False,
    )
    assert returned_fig is fig
    assert _axes((returned_fig, axes)) == list(supplied)
    assert all(np.asarray(axis.images[0].get_array()).shape[0] == 1 for axis in supplied)
    plt.close(fig)


def test_plot_samm_drops_unused_origin_levels_instead_of_drawing_zero_row_groups():
    samm = traminex.seqsamm(_ragged_fixture(), 3)
    fig, axes = _plot(samm)
    titles = [axis.get_title() for axis in _axes((fig, axes))]
    assert "Transition out of Gamma" not in titles
    assert "Transition out of %" not in titles
    plt.close(fig)


def test_plot_samm_no_transition_and_empty_results_raise_specific_errors():
    no_transitions = traminex.seqsamm(_constant_fixture(), 2)
    with pytest.raises(ValueError, match="no transition rows"):
        _plot(no_transitions)

    empty_frame = no_transitions.frame.iloc[:0].copy()
    empty = SAMMResult(
        frame=empty_frame,
        sequences=STSList(data=np.empty((0, 2), dtype=int), alphabet=["A", "B"]),
        spell=("A", "B"),
        typology=np.empty(0, dtype=object),
        state_columns=("s.1", "s.2"),
        sublength=2,
    )
    with pytest.raises(ValueError, match="no transition rows"):
        _plot(empty)


def test_plot_samm_validates_object_and_malformed_result_components():
    with pytest.raises(TypeError, match="SAMMResult"):
        _plot(object())

    valid = traminex.seqsamm(_repeated_fixture(), 2)
    malformed = [
        SAMMResult(
            frame=valid.frame.drop(columns="transition"),
            sequences=valid.sequences,
            spell=valid.spell,
            typology=valid.typology,
            state_columns=valid.state_columns,
            sublength=valid.sublength,
        ),
        SAMMResult(
            frame=valid.frame,
            sequences=valid.sequences,
            spell=valid.spell,
            typology=valid.typology[:-1],
            state_columns=valid.state_columns,
            sublength=valid.sublength,
        ),
        SAMMResult(
            frame=valid.frame,
            sequences=valid.sequences,
            spell=valid.spell,
            typology=valid.typology,
            state_columns=(),
            sublength=valid.sublength,
        ),
    ]
    for value in malformed:
        with pytest.raises(ValueError):
            _plot(value)


def test_plot_samm_invalid_type_is_actionable():
    with pytest.raises(ValueError, match="type.*i.*rf"):
        _plot(traminex.seqsamm(_repeated_fixture(), 2), type="unknown")


def test_plot_samm_is_nonmutating_and_top_level_identity_export():
    samm = traminex.seqsamm(_ragged_fixture(), 3)
    original_frame = samm.frame.copy(deep=True)
    original_data = samm.sequences.data.copy()
    original_typology = samm.typology.copy()
    fig, _ = _plot(samm)

    pd.testing.assert_frame_equal(samm.frame, original_frame)
    np.testing.assert_array_equal(samm.sequences.data, original_data)
    np.testing.assert_array_equal(samm.typology, original_typology)
    from traminex.visualization.plots import samm as plot_module

    assert traminex.plot_SAMM is plot_module.plot_SAMM
    assert "SAMMResult" not in traminex.__all__
    plt.close(fig)
