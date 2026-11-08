"""Executable public-series contract consumed by visualization siblings."""

from __future__ import annotations

import numpy as np
import pytest

from traminex import (
    STSList,
    seqdef,
    seqdss,
    seqmeant,
    seqmodst,
    seqrep,
    seqrf,
    seqstatd,
    seqtab,
    seqtrate,
)
from traminex.api_manifest import PUBLIC_API_SET
from traminex.visualization.series import resolve_sort_order


def _fixture() -> STSList:
    return seqdef(
        np.array(
            [
                ["A", "A", "B", "C"],
                ["A", "B", "B", "%"],
                ["B", "B", "*", "%"],
                ["C", "C", "C", "C"],
            ],
            dtype=object,
        ),
        alphabet=["A", "B", "C"],
        states=["Alpha", "Beta", "Gamma"],
        cpal=["#111111", "#222222", "#333333"],
        nr="*",
        void="%",
        right="DEL",
        weights=[1.0, 2.0, 1.0, 3.0],
        id=["a", "b", "c", "d"],
        cnames=["t1", "t2", "t3", "t4"],
        start=7,
        xtstep=2,
        tick_last=True,
    )


def _distance() -> np.ndarray:
    return np.array(
        [
            [0.0, 1.0, 2.0, 3.0],
            [1.0, 0.0, 1.5, 2.5],
            [2.0, 1.5, 0.0, 3.0],
            [3.0, 2.5, 3.0, 0.0],
        ]
    )


def test_stslist_shape_order_codebook_and_metadata_contract():
    x = _fixture()

    assert x.data.shape == (4, 4)
    assert x.data.dtype.kind in "iu"
    assert x.alphabet == ["Alpha", "Beta", "Gamma", "*", "%"]
    assert x.labels == x.alphabet
    assert x.cpal == ["#111111", "#222222", "#333333", "darkgrey", "white"]
    assert (x.nr_code, x.void_code) == (3, 4)
    np.testing.assert_array_equal(x.weights, [1.0, 2.0, 1.0, 3.0])
    assert x.ids == ["a", "b", "c", "d"]
    assert (x.start, x.xtlab, x.xtstep, x.tick_last) == (
        7,
        ["t1", "t2", "t3", "t4"],
        2,
        True,
    )
    np.testing.assert_array_equal(resolve_sort_order(x, None), np.arange(4))


def test_distinct_successive_state_output_preserves_rows_states_and_case_metadata():
    x = _fixture()
    dss = seqdss(x)

    assert isinstance(dss, STSList)
    assert dss.data.shape == (4, 3)
    np.testing.assert_array_equal(
        dss.data,
        np.array([[0, 1, 2], [0, 1, 4], [1, 4, 4], [2, 4, 4]]),
    )
    assert dss.alphabet == x.alphabet
    assert dss.labels == x.labels
    assert dss.cpal == x.cpal
    assert dss.ids == x.ids
    np.testing.assert_array_equal(dss.weights, x.weights)
    assert (dss.nr_code, dss.void_code, dss.start) == (x.nr_code, x.void_code, x.start)


def test_seqstatd_orientation_denominators_missing_policy_and_entropy_contract():
    x = _fixture()
    result = seqstatd(x)

    assert list(result) == ["Frequencies", "ValidStates", "Entropy"]
    assert result["Frequencies"].shape == (x.seqlen, len(x.alphabet))
    assert result["ValidStates"].shape == (x.seqlen,)
    assert result["Entropy"].shape == (x.seqlen,)
    np.testing.assert_allclose(result["ValidStates"], [7.0, 7.0, 6.0, 4.0])
    np.testing.assert_allclose(result["Frequencies"].sum(axis=1), np.ones(4))
    np.testing.assert_array_equal(result["Frequencies"][:, [x.nr_code, x.void_code]], 0)

    with_missing = seqstatd(x, include_missing=True)
    assert with_missing["ValidStates"][2] == 7.0
    assert with_missing["Frequencies"][2, x.nr_code] == pytest.approx(1 / 7)
    unweighted = seqstatd(x, weighted=False)
    np.testing.assert_allclose(unweighted["ValidStates"], [4.0, 4.0, 3.0, 2.0])


def test_transition_rate_orientation_is_origin_by_destination_in_state_order():
    rates = seqtrate(_fixture(), as_dataframe=False)
    assert rates.shape == (3, 3)
    np.testing.assert_allclose(
        rates,
        np.array([[0.25, 0.75, 0.0], [0.0, 0.75, 0.25], [0.0, 0.0, 1.0]]),
    )
    np.testing.assert_allclose(rates.sum(axis=1), 1.0)


def test_mean_modal_and_tabulation_outputs_keep_state_time_and_frequency_order():
    x = _fixture()
    mean = seqmeant(x, serr=True)
    modal = seqmodst(x)
    tab = seqtab(x)

    assert list(mean.columns) == ["Mean", "Var", "Stdev", "SE"]
    assert list(mean.index) == x.alphabet
    assert isinstance(modal, STSList)
    assert modal.data.shape == (1, x.seqlen)
    np.testing.assert_array_equal(modal.data[0], [0, 1, 1, 2])
    assert modal.frequencies_.shape == (x.seqlen,)
    assert modal.nbseq_ == pytest.approx(x.weights.sum())
    assert list(tab.columns) == ["sequence", "count", "prop", "cumprop"]
    assert tab.iloc[0]["sequence"] == "Gamma/4"
    assert tab.iloc[0]["count"] == 3.0
    assert tab.iloc[-1]["cumprop"] == pytest.approx(1.0)
    assert tab.attrs["format"] == "SPS"
    assert tab.attrs["source_indices"][0] == 3
    np.testing.assert_array_equal(tab.attrs["representatives"].data[0], x.data[3])


def test_representative_and_rf_orientations_and_weighted_box_inputs():
    x = _fixture()
    distance = _distance()
    representatives = seqrep(x, diss=distance, nrep=2)
    relative = seqrf(x, distance, k=2, sortv=np.arange(x.nseq), weighted=True)

    assert representatives["representatives"].shape == (2, x.seqlen)
    assert representatives["rep_indices"].shape == (2,)
    assert len(representatives["statistics"]) == 3  # representatives + Total
    assert relative["seqtoplot"].shape == (2, x.seqlen)
    assert relative["medoids"].shape == (2,)
    np.testing.assert_allclose(relative["heights"].sum(), 1.0)

    rf = relative["rf"]
    for indexes, distances, weights in zip(
        rf["index_list"], rf["dist_list"], rf["weights_list"], strict=True
    ):
        assert len(indexes) == len(distances) == len(weights)
        assert np.isfinite(distances).all()
        assert np.isfinite(weights).all()
        assert (weights >= 0).all()


def test_public_series_calls_are_nonmutating():
    x = _fixture()
    original = x.copy()
    distance = _distance()
    original_distance = distance.copy()

    seqdss(x)
    seqstatd(x)
    seqtrate(x, as_dataframe=False)
    seqmeant(x, serr=True)
    seqmodst(x)
    seqtab(x)
    seqrep(x, diss=distance, nrep=2)
    seqrf(x, distance, k=2, sortv=np.arange(x.nseq))

    np.testing.assert_array_equal(x.data, original.data)
    np.testing.assert_array_equal(x.weights, original.weights)
    assert x.alphabet == original.alphabet
    assert x.ids == original.ids
    np.testing.assert_array_equal(distance, original_distance)


def test_contract_refuses_malformed_core_inputs_and_keeps_private_helpers_private():
    with pytest.raises(ValueError, match="2D"):
        STSList(data=np.array([0, 1]), alphabet=["A", "B"])
    with pytest.raises(ValueError, match="ids"):
        STSList(data=np.zeros((2, 2), dtype=int), alphabet=["A"], ids=["one"])
    with pytest.raises(ValueError, match="weights"):
        seqstatd(_fixture(), weights=np.ones(2))
    with pytest.raises(ValueError, match="square"):
        seqrf(_fixture(), np.ones((3, 4)), k=2)
    with pytest.raises(TypeError, match="state sequence"):
        seqrf(object(), _distance(), k=2)

    assert "_wtd_boxplot_stats" not in PUBLIC_API_SET
    assert "py_ggplotseq" not in PUBLIC_API_SET
