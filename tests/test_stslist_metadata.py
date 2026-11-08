"""Identity and missing-color metadata contracts for ``STSList``."""

from __future__ import annotations

from importlib import import_module

import numpy as np
import pandas as pd
import pytest
from matplotlib import pyplot as plt

import traminex.distances as distances
import traminex.stats_utils as stats_utils
from traminex._sequence_builders import build_sts_from_state_matrix
from traminex.events._helpers import EventSequence, events_to_state_sequences, frequent_ngrams
from traminex.format import STS_to_SRS, seqdef
from traminex.seqformat import STS_to_TSE, seqformat, seqstart, to_DSS, to_SPS, toPersonPeriod
from traminex.sequences import STSList
from traminex.stats import seqistatd, seqstatd
from traminex.stats_complexity import seqdss
from traminex.stats_p1 import seqmodst, seqstatf, seqstatl
from traminex.stats_p2 import seqconc, seqgranularity, seqimplic, seqpolyads, seqrecode, seqsha
from traminex.visualization.plots.summary import _subset_stslist

chronogram_plots = import_module("traminex.visualization.plots.chronogram")


def _rich_stslist() -> STSList:
    return STSList(
        data=np.array(
            [
                [0, 0, 1, 1],
                [1, 1, 0, 0],
                [0, 1, 0, 1],
            ],
            dtype=int,
        ),
        alphabet=["A", "B", "%"],
        labels=["Alpha", "Beta", "Void"],
        cpal=["#111111", "#eeeeee", "#ffffff"],
        xtlab=["2000", "2001", "2002", "2003"],
        xtstep=2,
        tick_last=True,
        nr="?",
        void="%",
        nr_code=None,
        void_code=2,
        weights=np.array([0.5, 1.5, 2.0]),
        start=2000,
        ids=["case-a", 7, ("case", "c")],
        missing_color="#aabbcc",
        issues={
            "truncated": [1, 3],
            "empty.seq": [2],
            "start.before.birth": [],
            "spell.before.birth": [3],
        },
    )


def _sparse_external_marker_stslist() -> STSList:
    return STSList(
        np.array([[0, 3, 1, 4], [1, 1, 4, 4]], dtype=int),
        alphabet=["A", "B"],
        labels=["Alpha", "Beta"],
        cpal=["red", "blue"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
        ids=["p1", "p2"],
        missing_color="gray",
    )


def test_stslist_metadata_defaults_and_validation():
    default = STSList(np.array([[0, 0]], dtype=int), alphabet=["A"])
    assert default.ids is None
    assert default.missing_color == "darkgrey"
    assert default.issues is None

    source_ids = ["one", 2]
    identified = STSList(
        np.array([[0], [0]], dtype=int),
        alphabet=["A"],
        ids=source_ids,
    )
    source_ids[0] = "changed"
    assert identified.ids == ["one", 2]

    with pytest.raises(ValueError, match="ids must have length 2"):
        STSList(
            np.array([[0], [0]], dtype=int),
            alphabet=["A"],
            ids=["only-one"],
        )
    with pytest.raises(TypeError, match="missing_color must be a string"):
        STSList(
            np.array([[0]], dtype=int),
            alphabet=["A"],
            missing_color=123,  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("dtype", [np.int8, np.int64, np.uint8, np.uint64])
def test_stslist_accepts_integral_code_dtypes(dtype):
    sequences = STSList(np.array([[0, 1]], dtype=dtype), alphabet=["A", "B"])

    assert sequences.to_strings() == ["A-B"]


@pytest.mark.parametrize(
    "data",
    [
        np.array([[False, True]], dtype=bool),
        np.array([[0.0, 1.0]], dtype=float),
        np.array([[0.0, 1.5]], dtype=float),
        np.array([[0, 1]], dtype=object),
        np.array([[0, "1"]], dtype=object),
    ],
)
def test_stslist_refuses_non_integral_code_dtypes(data):
    with pytest.raises(TypeError, match="data must have an integral dtype"):
        STSList(data, alphabet=["A", "B"])


@pytest.mark.parametrize(
    ("data", "message"),
    [
        (np.array([[-1]], dtype=int), "non-negative"),
        (np.array([[2]], dtype=int), "undeclared state code 2"),
        (np.array([[99]], dtype=np.uint64), "undeclared state code 99"),
    ],
)
def test_stslist_refuses_values_outside_the_declared_code_domain(data, message):
    with pytest.raises(ValueError, match=message):
        STSList(data, alphabet=["A", "B"])


@pytest.mark.parametrize("attribute", ["nr_code", "void_code"])
@pytest.mark.parametrize("value", [True, 1.0, "1"])
def test_stslist_refuses_non_integral_marker_codes(attribute, value):
    kwargs = {attribute: value}

    with pytest.raises(TypeError, match=rf"{attribute} must be an integer"):
        STSList(np.array([[0]], dtype=int), alphabet=["A"], **kwargs)


@pytest.mark.parametrize("attribute", ["nr_code", "void_code"])
def test_stslist_refuses_negative_marker_codes(attribute):
    kwargs = {attribute: -1}

    with pytest.raises(ValueError, match=rf"{attribute} must be non-negative"):
        STSList(np.array([[0]], dtype=int), alphabet=["A"], **kwargs)


def test_stslist_refuses_colliding_missing_and_void_codes():
    with pytest.raises(ValueError, match="nr_code and void_code must be distinct"):
        STSList(
            np.array([[0, 1]], dtype=int),
            alphabet=["A", "marker"],
            nr_code=1,
            void_code=1,
        )


def test_stslist_accepts_sparse_out_of_alphabet_marker_codes():
    sequences = STSList(
        np.array([[0, 3, 1, 4]], dtype=int),
        alphabet=["A", "B"],
        nr_code=3,
        void_code=4,
    )

    assert sequences.alphabet == ["A", "B"]
    assert sequences.to_strings() == ["A-*-B-%"]


def test_stslist_preserves_valid_out_of_alphabet_markers_without_recoding_data():
    data = np.array([[0, 2, 3], [1, 1, 3]], dtype=int)
    sequences = STSList(
        data=data,
        alphabet=["A", "B"],
        labels=["Alpha", "Beta"],
        nr="?",
        void="!",
        nr_code=2,
        void_code=3,
    )

    assert sequences.alphabet == ["A", "B"]
    assert sequences.labels == ["Alpha", "Beta"]
    assert sequences.nr_code == 2 and sequences.void_code == 3
    np.testing.assert_array_equal(sequences.data, data)
    assert sequences.to_strings() == ["A-?-!", "B-B-!"]


def test_stslist_preserves_declared_but_unobserved_marker_codes():
    sequences = STSList(
        np.array([[0, 1]], dtype=int),
        alphabet=["A", "B"],
        nr_code=2,
        void_code=3,
    )

    assert sequences.alphabet == ["A", "B"]
    assert sequences.to_strings() == ["A-B"]


@pytest.mark.parametrize("dtype", [np.int8, np.uint8])
def test_unobserved_external_marker_codes_do_not_depend_on_matrix_dtype_range(dtype):
    sequences = STSList(
        np.array([[0]], dtype=dtype),
        alphabet=["A"],
        nr_code=300,
        void_code=301,
    )

    assert sequences.nr_code == 300
    assert sequences.void_code == 301
    assert sequences.to_strings() == ["A"]


def test_stslist_decoder_uses_declared_marker_roles_in_loose_codebooks():
    sequences = STSList(
        np.array([[0, 1, 2]], dtype=int),
        alphabet=["A", "encoded missing", "encoded void"],
        nr="*",
        void="%",
        nr_code=1,
        void_code=2,
    )

    assert sequences.to_strings() == ["A-*-%"]


@pytest.mark.parametrize(
    "data",
    [
        np.empty((0, 3), dtype=int),
        np.empty((3, 0), dtype=int),
        np.empty((0, 0), dtype=int),
    ],
)
def test_stslist_preserves_valid_empty_shapes(data):
    sequences = STSList(data, alphabet=[])

    assert sequences.data.shape == data.shape
    assert sequences.to_strings() == [""] * data.shape[0]


def test_valid_out_of_alphabet_markers_are_safe_across_real_consumers():
    sequences = STSList(
        np.array([[0, 3, 1, 4], [1, 1, 4, 4]], dtype=int),
        alphabet=["A", "B"],
        nr_code=3,
        void_code=4,
        ids=["p1", "p2"],
    )

    assert sequences.to_strings() == ["A-*-B-%", "B-B-%-%"]
    # Marker-padded like the in-alphabet twin ``A, B, *, %`` (void column always
    # zero): per-state shape must not depend on the marker representation
    # (2026-09-22 audit F-04; this used to pin the unpadded ``A, B, *``).
    np.testing.assert_array_equal(
        seqistatd(sequences, include_missing=True),
        np.array([[1, 1, 1, 0], [0, 2, 0, 0]]),
    )
    assert seqstatf(sequences, include_missing=True)["Freq"].tolist() == [1.0, 3.0, 1.0]
    assert frequent_ngrams(sequences, n=1, min_support=1) == [
        (("%",), 3, 3 / 8),
        (("B",), 3, 3 / 8),
        (("*",), 1, 1 / 8),
        (("A",), 1, 1 / 8),
    ]

    shifted = STS_to_SRS(sequences)
    assert shifted.shape[0] == 8
    assert set(shifted.iloc[:, 2:].stack().dropna()) == {"A", "B"}

    person_period = toPersonPeriod(sequences)
    assert person_period["id"].tolist() == ["p1", "p1", "p1", "p2", "p2"]
    assert person_period["state"].tolist()[0] == "A"
    assert pd.isna(person_period["state"].tolist()[1])
    assert person_period["state"].tolist()[2:] == ["B", "B", "B"]

    copied = sequences.copy()
    subset = _subset_stslist(sequences, np.array([1]))
    rebuilt = seqdef(sequences)
    assert copied.alphabet == subset.alphabet == sequences.alphabet == ["A", "B"]
    assert copied.nr_code == subset.nr_code == 3
    assert copied.void_code == subset.void_code == 4
    assert rebuilt.to_strings() == sequences.to_strings()


def test_sparse_external_markers_are_decoded_by_public_text_converters():
    sequences = _sparse_external_marker_stslist()

    # Void cells get no record (TX-AUD-20261019-082); the missing cell keeps its marker.
    assert STS_to_TSE(sequences, state_changes_only=False)["state"].tolist() == [
        "A",
        "*",
        "B",
        "B",
        "B",
    ]
    # Void cells are left out, as in R's SPS and DSS forms
    # (tests/test_r_goldens_sts_sps_dss_void.py, TX-AUD-20261010-033).
    assert to_SPS(sequences) == ["A/1-*/1-B/1", "B/2"]
    assert to_DSS(sequences) == ["A-*-B", "B"]
    assert seqstatl(sequences) == ["*", "%", "A", "B"]
    assert seqconc(sequences).tolist() == ["A-*-B", "B-B"]


def test_sparse_external_markers_are_decoded_by_seqformat_routes():
    sequences = _sparse_external_marker_stslist()
    states = ["A", "B", "*", "%"]
    transition_events = pd.DataFrame(None, index=states, columns=states, dtype=object)
    for source in states:
        for target in states:
            transition_events.loc[source, target] = (
                f"start-{source}" if source == target else f"{source}>{target}"
            )

    with pytest.warns(UserWarning, match="'missing' set"):
        events = seqformat(
            sequences,
            from_="STS",
            to="TSE",
            tevent=transition_events,
            id=["p1", "p2"],
        )
    assert isinstance(events, pd.DataFrame)
    assert events.to_dict("list") == {
        "id": ["p1", "p1", "p1", "p1", "p2", "p2"],
        "time": [0, 1, 2, 3, 0, 2],
        "event": ["start-A", "A>*", "*>B", "B>%", "start-B", "B>%"],
    }

    with pytest.warns(UserWarning, match="'missing' set"):
        shifted = seqformat(sequences, from_="STS", to="SRS", covar=[2])
    assert isinstance(shifted, pd.DataFrame)
    assert shifted["T2"].astype(object).tolist() == ["*", "*", "*", "*", "B", "B", "B", "B"]


def test_sparse_external_markers_are_preserved_by_transform_consumers():
    sequences = _sparse_external_marker_stslist()

    recoded = seqrecode(sequences, {"AB": ["A", "B"]})
    assert recoded.to_strings() == ["AB-*-AB-%", "AB-AB-%-%"]
    assert recoded.alphabet == ["AB", "*", "%"]
    assert recoded.nr_code == 1
    assert recoded.void_code == 2

    history = seqsha(sequences, [4, 4], [True, False])
    assert history.loc[3, ["T1", "T2", "T3", "T4"]].tolist() == ["A", "*", "B", "%"]
    assert history.loc[7, ["T1", "T2", "T3", "T4"]].tolist() == ["B", "B", "%", "%"]

    implication = seqimplic(sequences, ["g1", "g2"], with_missing=True)
    assert implication.alphabet == ["A", "B", "*"]
    assert implication.labels == ["Alpha", "Beta", "Missing"]
    assert implication.cpal == ["red", "blue", "gray"]
    assert implication.indices.shape == (2, 3, 4)


def test_seqstatd_counts_sparse_external_missing_state():
    sequences = _sparse_external_marker_stslist()

    result = seqstatd(sequences, include_missing=True, norm=False)

    # Marker-padded like the in-alphabet twin ``A, B, *, %``: the trailing void
    # column is zero, including at the all-void position (2026-09-22 audit F-04;
    # this used to pin the unpadded three columns).
    np.testing.assert_allclose(
        result["Frequencies"],
        np.array(
            [
                [0.5, 0.5, 0.0, 0.0],
                [0.0, 0.5, 0.5, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [np.nan, np.nan, np.nan, 0.0],
            ]
        ),
        rtol=0,
        atol=0,
        equal_nan=True,
    )
    np.testing.assert_array_equal(result["ValidStates"], [2.0, 2.0, 1.0, 0.0])
    assert result.cpal == ["red", "blue", "gray"]


def test_seqmodst_maps_distribution_columns_back_to_sparse_marker_codes():
    sequences = STSList(
        np.array([[3], [3], [0]], dtype=int),
        alphabet=["A"],
        nr="*",
        nr_code=3,
    )

    modal = seqmodst(sequences, include_missing=True)

    assert modal.to_strings() == ["*"]
    assert modal.nr_code is not None
    np.testing.assert_allclose(modal.frequencies_, [2 / 3], rtol=0, atol=1e-15)
    assert modal.occurrences_ == 2


def test_stslist_issues_validation_and_defensive_copy():
    source = {
        "truncated": [1, 2],
        "empty.seq": [],
    }
    sequences = STSList(
        np.array([[0], [0]], dtype=int),
        alphabet=["A"],
        issues=source,
    )
    source["truncated"][0] = 2
    source["empty.seq"].append(1)
    source["new"] = [1]
    assert sequences.issues == {"truncated": [1, 2], "empty.seq": []}
    assert sequences.issues is not source
    assert sequences.issues["truncated"] is not source["truncated"]

    with pytest.raises(TypeError, match="issues must be a dictionary"):
        STSList(
            np.array([[0]], dtype=int),
            alphabet=["A"],
            issues=[],  # type: ignore[arg-type]
        )
    with pytest.raises(TypeError, match="issues keys must be strings"):
        STSList(
            np.array([[0]], dtype=int),
            alphabet=["A"],
            issues={1: [1]},  # type: ignore[dict-item]
        )
    with pytest.raises(TypeError, match="issues values must be lists of sequence indexes"):
        STSList(
            np.array([[0]], dtype=int),
            alphabet=["A"],
            issues={"truncated": (1,)},  # type: ignore[dict-item]
        )
    with pytest.raises(TypeError, match="issue indexes must be integers"):
        STSList(
            np.array([[0]], dtype=int),
            alphabet=["A"],
            issues={"truncated": [True]},
        )
    with pytest.raises(ValueError, match="issue indexes must be between 1 and 2"):
        STSList(
            np.array([[0], [0]], dtype=int),
            alphabet=["A"],
            issues={"truncated": [0, 3]},
        )


def test_stslist_copy_preserves_all_metadata_with_isolation():
    original = _rich_stslist()
    copied = original.copy()

    assert copied.ids == original.ids
    assert copied.ids is not original.ids
    assert copied.issues == original.issues
    assert copied.issues is not original.issues
    assert copied.issues["truncated"] is not original.issues["truncated"]  # type: ignore[index]
    assert copied.missing_color == "#aabbcc"
    assert copied.labels == original.labels
    assert copied.cpal == original.cpal
    assert copied.xtlab == original.xtlab
    assert copied.xtstep == 2
    assert copied.tick_last is True
    assert copied.start == 2000
    assert copied.nr == "?" and copied.void == "%"
    assert copied.nr_code is None and copied.void_code == 2
    np.testing.assert_array_equal(copied.weights, original.weights)
    assert copied.weights is not original.weights

    copied.data[0, 0] = 1
    copied.alphabet[0] = "changed"
    copied.labels[0] = "changed"  # type: ignore[index]
    copied.cpal[0] = "changed"  # type: ignore[index]
    copied.xtlab[0] = "changed"  # type: ignore[index]
    copied.weights[0] = 99  # type: ignore[index]
    copied.ids[0] = "changed"  # type: ignore[index]
    copied.issues["truncated"][0] = 2  # type: ignore[index]
    assert original.data[0, 0] == 0
    assert original.alphabet[0] == "A"
    assert original.labels[0] == "Alpha"  # type: ignore[index]
    assert original.cpal[0] == "#111111"  # type: ignore[index]
    assert original.xtlab[0] == "2000"  # type: ignore[index]
    assert original.weights[0] == 0.5  # type: ignore[index]
    assert original.ids[0] == "case-a"  # type: ignore[index]
    assert original.issues["truncated"] == [1, 3]  # type: ignore[index]


@pytest.mark.parametrize(
    "transform",
    [
        lambda x: seqdss(x),
        lambda x: seqgranularity(x, tspan=2, method="last"),
        lambda x: seqstart(x, data_start=2000, new_start=2000),
    ],
)
def test_row_preserving_transforms_keep_identity_metadata(transform):
    original = _rich_stslist()
    transformed = transform(original)
    assert transformed.ids == original.ids
    assert transformed.ids is not original.ids
    assert transformed.missing_color == original.missing_color


@pytest.mark.parametrize(
    "transform",
    [
        lambda x: seqgranularity(x, tspan=1, method="last"),
        lambda x: seqrecode(x, {"Combined": ["A", "B"]}),
    ],
)
def test_row_preserving_transforms_keep_axis_metadata(transform):
    original = _rich_stslist()

    transformed = transform(original)

    assert transformed.xtlab == original.xtlab
    assert transformed.xtlab is not original.xtlab
    assert transformed.xtstep == original.xtstep
    assert transformed.tick_last is original.tick_last


def test_seqdss_uses_width_aligned_spell_axis_metadata():
    original = _rich_stslist()
    original.data = np.zeros_like(original.data)

    transformed = seqdss(original)

    assert transformed.seqlen == 1
    assert transformed.xtlab == ["ST1"]
    assert len(transformed.xtlab) == transformed.seqlen
    assert transformed.xtstep == original.xtstep
    assert transformed.tick_last is original.tick_last


def test_seqgranularity_coarsens_xtstep_and_preserves_tick_last():
    original = _rich_stslist()
    original.xtstep = 3

    transformed = seqgranularity(original, tspan=2, method="last")

    assert transformed.xtstep == 2
    assert transformed.tick_last is True


def test_seqstart_rebuilds_axis_labels_from_the_aligned_grid():
    original = STSList(
        data=np.array([[0, 1], [1, 0]], dtype=int),
        alphabet=["A", "B"],
        xtlab=["2000", "2001"],
        ids=["case-a", "case-b"],
    )

    transformed = seqstart(
        original,
        data_start=[2000, 2000],
        new_start=[2000, 1999],
    )

    assert transformed.data.shape == (2, 3)
    assert transformed.xtlab == ["0", "1", "2"]
    assert len(transformed.xtlab) == transformed.seqlen


def test_seqrecode_preserves_ids_missing_color_and_issues():
    original = _rich_stslist()

    transformed = seqrecode(original, {"Combined": ["A", "B"]})

    assert transformed.ids == original.ids
    assert transformed.ids is not original.ids
    assert transformed.missing_color == original.missing_color
    assert transformed.issues == original.issues
    assert transformed.issues is not original.issues
    assert transformed.issues["truncated"] is not original.issues["truncated"]  # type: ignore[index]


@pytest.mark.parametrize(
    "transform",
    [
        lambda x: seqdss(x),
        lambda x: seqgranularity(x, tspan=2, method="last"),
        lambda x: seqstart(x, data_start=2000, new_start=2000),
        lambda x: seqrecode(x, {"Combined": ["A", "B"]}),
    ],
)
def test_row_preserving_transforms_keep_issues_with_isolation(transform):
    original = _rich_stslist()
    transformed = transform(original)
    assert transformed.issues == original.issues
    assert transformed.issues is not original.issues
    assert transformed.issues["truncated"] is not original.issues["truncated"]  # type: ignore[index]


def test_row_transforms_align_state_metadata_when_special_states_are_generated():
    no_void_metadata = STSList(
        data=np.array([[0, 0, 1], [1, 1, 1]], dtype=int),
        alphabet=["A", "B"],
        labels=["Alpha", "Beta"],
        cpal=["red", "blue"],
        ids=["left", "right"],
        missing_color="silver",
    )
    dss = seqdss(no_void_metadata)
    assert dss.alphabet == ["A", "B", "%"]
    assert dss.labels == ["Alpha", "Beta", "%"]
    assert dss.cpal == ["red", "blue", "white"]
    assert dss.ids == ["left", "right"]

    missing_block = STSList(
        data=np.array([[2, 2, 0, 0], [0, 0, 2, 2]], dtype=int),
        alphabet=["A", "B", "%"],
        labels=["Alpha", "Beta", "Void"],
        cpal=["red", "blue", "white"],
        nr="?",
        void="%",
        nr_code=None,
        void_code=2,
        ids=["left", "right"],
        missing_color="silver",
    )
    granular = seqgranularity(missing_block, tspan=2, method="first.valid")
    assert granular.alphabet == ["A", "B", "%", "?"]
    assert granular.labels == ["Alpha", "Beta", "Void", "?"]
    assert granular.cpal == ["red", "blue", "white", "silver"]
    assert granular.ids == ["left", "right"]


def test_row_subset_helper_subsets_ids_weights_and_color():
    original = _rich_stslist()
    subset = _subset_stslist(original, np.array([2, 0]))
    assert subset.ids == [("case", "c"), "case-a"]
    assert subset.missing_color == "#aabbcc"
    assert subset.issues == {
        "truncated": [1, 2],
        "empty.seq": [],
        "start.before.birth": [],
        "spell.before.birth": [1],
    }
    np.testing.assert_allclose(subset.weights, [2.0, 0.5])

    masked = _subset_stslist(original, np.array([True, False, True]))
    assert masked.ids == ["case-a", ("case", "c")]
    assert masked.issues == {
        "truncated": [1, 2],
        "empty.seq": [],
        "start.before.birth": [],
        "spell.before.birth": [2],
    }
    np.testing.assert_allclose(masked.weights, [0.5, 2.0])


def test_chronogram_subsets_follow_cluster_row_order(monkeypatch):
    original = _rich_stslist()
    captured: list[STSList] = []

    def capture_stslist(**kwargs):
        subset = STSList(**kwargs)
        captured.append(subset)
        return subset

    monkeypatch.setattr(chronogram_plots, "STSList", capture_stslist)
    labels = np.array([2, 1, 2])

    fig, _ = chronogram_plots.cluster_chronograms(original, labels)
    assert [subset.ids for subset in captured] == [[7], ["case-a", ("case", "c")]]
    assert [subset.issues for subset in captured] == [
        {
            "truncated": [],
            "empty.seq": [1],
            "start.before.birth": [],
            "spell.before.birth": [],
        },
        {
            "truncated": [1, 2],
            "empty.seq": [],
            "start.before.birth": [],
            "spell.before.birth": [2],
        },
    ]
    assert all(subset.missing_color == "#aabbcc" for subset in captured)
    plt.close(fig)

    captured.clear()
    fig, _ = chronogram_plots.cluster_chronograms_with_representatives(
        original,
        labels,
        np.zeros((original.nseq, original.nseq)),
    )
    assert [subset.ids for subset in captured] == [[7], ["case-a", ("case", "c")]]
    assert [subset.issues for subset in captured] == [
        {
            "truncated": [],
            "empty.seq": [1],
            "start.before.birth": [],
            "spell.before.birth": [],
        },
        {
            "truncated": [1, 2],
            "empty.seq": [],
            "start.before.birth": [],
            "spell.before.birth": [2],
        },
    ]
    assert all(subset.missing_color == "#aabbcc" for subset in captured)
    plt.close(fig)


def test_seqpolyads_stack_has_no_ambiguous_row_ids(monkeypatch):
    first = _rich_stslist()
    second = _rich_stslist()
    second.ids = ["role-2-a", "role-2-b", "role-2-c"]
    second.missing_color = "orange"
    captured: list[STSList] = []

    def fake_seqdist(x, **_kwargs):
        captured.append(x)
        return np.zeros((x.nseq, x.nseq))

    monkeypatch.setattr(distances, "seqdist", fake_seqdist)
    seqpolyads([first, second], T=2)

    assert len(captured) == 1
    assert captured[0].ids is None
    assert captured[0].issues is None
    assert captured[0].missing_color == first.missing_color


def test_builders_and_seqdef_attach_explicit_ids_and_missing_color():
    matrix = np.array([["A", "B"], ["B", "A"]], dtype=object)
    kwargs = {"ids": ["u1", "u2"], "missing_color": "orange"}
    from_states = build_sts_from_state_matrix(matrix, kwargs=kwargs)
    assert from_states.ids == ["u1", "u2"]
    assert from_states.missing_color == "orange"

    defined = seqdef(matrix, id=["person-1", "person-2"])
    assert defined.ids == ["person-1", "person-2"]


@pytest.mark.parametrize("align_missing", [False, True])
def test_events_to_state_sequences_uses_event_ids_as_fresh_row_ids(align_missing):
    events = EventSequence(
        id=["z", "a", "z"],
        time=np.array([1, 1, 2]),
        label=["A", "B", "C"],
    )
    result = events_to_state_sequences(
        events,
        time_grid=np.array([1, 2]),
        align_missing=align_missing,
    )
    assert result.ids == ["a", "z"]
    assert result.missing_color == "darkgrey"


def test_seqindic_dyn_windows_preserve_identity_metadata(monkeypatch):
    original = _rich_stslist()
    seen: list[STSList] = []

    def fake_seqindic(x, **_kwargs):
        seen.append(x)
        return pd.DataFrame({"Cplx": np.zeros(x.nseq)})

    monkeypatch.setattr(stats_utils, "seqindic", fake_seqindic)
    result = stats_utils.seqindic_dyn(original, window_size=2, wstep=1)
    assert result.shape == (3, 3)
    assert len(seen) == 3
    expected_xtlab = [
        ["2000", "2001"],
        ["2001", "2002"],
        ["2002", "2003"],
    ]
    for window, labels in zip(seen, expected_xtlab, strict=True):
        assert window.ids == original.ids
        assert window.ids is not original.ids
        assert window.missing_color == "#aabbcc"
        assert window.issues == original.issues
        assert window.issues is not original.issues
        assert window.xtlab == labels
        assert window.xtstep == 2
        assert window.tick_last is True
