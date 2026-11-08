"""TX-AUD-20261005-011: marker positions must not change analysis.

The factory case and literal unweighted expectations come from the pinned
TraMineR 2.2-14 / TraMineRextras 0.6.9 oracle retained in audit-only-20261005.
EMLT's label comparisons exclude missing and void; seqlogp includes missing
when requested and indexes probability models in ordinary-state-first order.
"""

import numpy as np
import pytest

import traminex as tx

_DISJUNCTIVE = np.array(
    [
        [0, 0, 1, 0, 1, 0, 0, 1],
        [0, 0, 0, 1, 1, 0, 0, 1],
        [0, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 0, 1, 0, 1, 1, 0],
    ],
    dtype=float,
)
# Only nonzero components are needed for the sign/rotation-invariant Gram
# matrix. These are literal coordinates from emlt-R-coordinates.csv.
_R_COORDINATES = np.array(
    [
        [-3.9556382062758795, 1.9894387727998963, -0.0053603514969431387],
        [-3.9556382062758804, -1.9894387727998957, -0.0053603514969433608],
        [3.9556382062758804, 1.9894387727998957, 0.0053603514969435828],
        [3.9556382062758795, -1.9894387727998963, 0.0053603514969433608],
    ]
)
_WEIGHTS = np.array([0.5, 1.0, 1.5, 2.0])
_LOGICAL_ROWS = [
    [None, "a", "a", "b"],
    [None, "b", "a", "b"],
    [None, "a", "b", "a"],
    [None, "b", "b", "a"],
]


def _factory(states=None, *, missing=True, fill="_"):
    events = tx.EventSequence(
        ["p1", "p1", "p2", "p2", "p2", "p3", "p3", "p3", "p4", "p4"],
        np.array([1, 3, 1, 2, 3, 1, 2, 3, 1, 3]),
        ["a", "b", "b", "a", "b", "a", "b", "a", "b", "a"],
    )
    return tx.events_to_state_sequences(
        events,
        states=states,
        time_grid=np.arange(4) if missing else np.arange(1, 4),
        fill=fill,
        align_missing=True,
    )


_LAYOUTS = [
    pytest.param(None, id="factory-missing-first"),
    pytest.param(["a", "_", "b"], id="factory-missing-interleaved"),
    pytest.param(["a", "b", "_"], id="factory-missing-last"),
]


@pytest.mark.parametrize("states", _LAYOUTS)
def test_emlt_factory_markers_do_not_count_as_ordinary_situations(states):
    result = tx.seqemlt(_factory(states))
    assert result.states == ["a", "b"]
    np.testing.assert_array_equal(result.sit_freq, [0, 0, 2, 2, 2, 2, 2, 2])
    np.testing.assert_array_equal(result.disjunctive, _DISJUNCTIVE)
    np.testing.assert_allclose(
        result.coordinates @ result.coordinates.T,
        _R_COORDINATES @ _R_COORDINATES.T,
        rtol=0,
        atol=1e-10,
    )


@pytest.mark.parametrize("states", _LAYOUTS)
@pytest.mark.parametrize("weighted", [False, True])
def test_emlt_factory_fractional_weights_preserve_decoded_geometry(states, weighted):
    factory = _factory(states)
    factory.weights = _WEIGHTS.copy()
    canonical = tx.seqdef(_LOGICAL_ROWS, alphabet=["a", "b"], weights=_WEIGHTS)
    actual = tx.seqemlt(factory, weighted=weighted)
    expected = tx.seqemlt(canonical, weighted=weighted)
    np.testing.assert_array_equal(actual.disjunctive, _DISJUNCTIVE)
    np.testing.assert_array_equal(actual.sit_freq, [0, 0, 2, 2, 2, 2, 2, 2])
    # The pinned R implementation uses crossproducts of weighted indicators.
    np.testing.assert_allclose(
        actual.sit_transrate[2, [2, 4, 5]],
        [1.25, 0.125, 1.125] if weighted else [1.0, 0.5, 0.5],
        rtol=0,
        atol=1e-14,
    )
    np.testing.assert_allclose(
        actual.coordinates @ actual.coordinates.T,
        expected.coordinates @ expected.coordinates.T,
        rtol=0,
        atol=1e-10,
    )


@pytest.mark.parametrize("states", _LAYOUTS)
@pytest.mark.parametrize("prob,expected", [("trate", np.log(4)), ("freq", np.log(8))])
def test_logp_factory_layout_matches_pinned_r_models(states, prob, expected):
    actual = tx.seqlogp(_factory(states), include_missing=True, prob=prob)
    np.testing.assert_allclose(actual, np.full(4, expected), rtol=0, atol=1e-14)


@pytest.mark.parametrize("states", _LAYOUTS)
@pytest.mark.parametrize("weighted", [False, True])
def test_logp_factory_fractional_weights_preserve_path_probabilities(states, weighted):
    factory = _factory(states)
    factory.weights = _WEIGHTS.copy()
    # All four observed paths branch into unique trajectories. Their fitted
    # time-varying Markov probabilities are their relative case weights.
    expected = -np.log([0.1, 0.2, 0.3, 0.4]) if weighted else np.full(4, np.log(4))
    actual = tx.seqlogp(factory, include_missing=True, weighted=weighted)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)


def test_factory_without_missing_keeps_ordinary_emlt_and_likelihood():
    factory = _factory(missing=False)
    result = tx.seqemlt(factory)
    np.testing.assert_array_equal(result.sit_freq, [2, 2, 2, 2, 2, 2])
    np.testing.assert_array_equal(result.disjunctive, _DISJUNCTIVE[:, 2:])
    np.testing.assert_allclose(
        result.coordinates @ result.coordinates.T,
        _R_COORDINATES @ _R_COORDINATES.T,
        rtol=0,
        atol=1e-10,
    )
    np.testing.assert_allclose(tx.seqlogp(factory), np.full(4, np.log(4)), atol=1e-14)
    np.testing.assert_allclose(tx.seqlogp(factory, prob="freq"), np.full(4, np.log(8)), atol=1e-14)


def test_emlt_interleaved_void_excludes_padding_and_preserves_geometry():
    padded = tx.STSList(
        np.array([[0, 0, 2, 1], [2, 0, 2, 1], [0, 2, 0, 1], [2, 2, 0, 1]]),
        ["a", "%", "b"],
        void_code=1,
    )
    actual = tx.seqemlt(padded)
    np.testing.assert_array_equal(actual.sit_freq, [2, 2, 2, 2, 2, 2, 0, 0])
    np.testing.assert_array_equal(
        actual.disjunctive, np.column_stack([_DISJUNCTIVE[:, 2:], np.zeros((4, 2))])
    )
    np.testing.assert_allclose(
        actual.coordinates @ actual.coordinates.T,
        _R_COORDINATES @ _R_COORDINATES.T,
        rtol=0,
        atol=1e-10,
    )


def _exploratory_marker_pair(marker_kind):
    """Supported STSList objects with identical decoded states and markers."""
    if marker_kind == "both":
        data = np.array([[0, 2, 1, 3], [1, 1, 2, 3], [0, 0, 0, 3], [1, 0, 1, 3]])
        markers = {"nr_code": 2, "void_code": 3}
        canonical_alphabet = ["A", "B", "*", "%"]
    elif marker_kind == "missing":
        data = np.array([[0, 2, 1], [1, 1, 2], [0, 0, 0], [1, 0, 1]])
        markers = {"nr_code": 2}
        canonical_alphabet = ["A", "B", "*"]
    else:
        data = np.array([[0, 0, 1, 2], [1, 1, 1, 2], [0, 0, 0, 2], [1, 0, 1, 2]])
        markers = {"void_code": 2}
        canonical_alphabet = ["A", "B", "%"]
    external = tx.STSList(data.copy(), ["A", "B"], **markers)
    canonical = tx.STSList(data.copy(), canonical_alphabet, **markers)
    assert external.to_strings() == canonical.to_strings()
    return external, canonical


_TREE_FEATURE_KINDS = ["durations", "transitions", "summary", "durations+transitions+summary"]


@pytest.mark.parametrize("marker_kind", ["missing", "void", "both"])
@pytest.mark.parametrize("kind", _TREE_FEATURE_KINDS)
def test_exploratory_tree_features_match_for_declared_external_markers(marker_kind, kind):
    # This pins representation consistency for the Python exploratory feature
    # contract (the features leave the markers out; TX-AUD-20261010-081).
    from traminex.seqtree_explore import _features_from_sequences

    external, canonical = _exploratory_marker_pair(marker_kind)
    actual, actual_names = _features_from_sequences(external, kind=kind)
    expected, expected_names = _features_from_sequences(canonical, kind=kind)
    assert actual.shape[1] == len(actual_names)
    assert actual_names == expected_names
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)


@pytest.mark.parametrize("kind", _TREE_FEATURE_KINDS)
def test_exploratory_tree_external_markers_preserve_fit_and_cross_layout_predict(kind):
    pytest.importorskip("sklearn")
    external, canonical = _exploratory_marker_pair("both")
    # B B * % and A A A % have equal summary features (one state, no transition between two
    # states; batch-12 review M-2), so they share a label every feature kind can learn.
    labels = np.array([0, 1, 1, 0])
    options = {"features": kind, "min_samples_leaf": 1, "random_state": 0}
    external_result = tx.seqtree_fit(external, labels, **options)
    canonical_result = tx.seqtree_fit(canonical, labels, **options)
    assert external_result.feature_names == canonical_result.feature_names
    np.testing.assert_array_equal(external_result.leaf_index, canonical_result.leaf_index)
    for fitted in [external_result, canonical_result]:
        for new_data in [external, canonical]:
            np.testing.assert_array_equal(
                tx.seqtree_predict(fitted, new_data, features=kind), labels
            )


def test_exploratory_tree_default_features_accept_declared_external_markers():
    pytest.importorskip("sklearn")
    external, canonical = _exploratory_marker_pair("both")
    labels = np.array([0, 1, 0, 1])
    actual = tx.seqtree_fit(external, labels, min_samples_leaf=1)
    expected = tx.seqtree_fit(canonical, labels, min_samples_leaf=1)
    assert actual.feature_names == expected.feature_names
    np.testing.assert_array_equal(tx.seqtree_predict(actual, canonical), labels)
    np.testing.assert_array_equal(tx.seqtree_predict(expected, external), labels)


def test_exploratory_tree_canonical_marker_features_leave_the_markers_out():
    # The features are over the alphabet's states (TX-AUD-20261010-081): the
    # rows are A * B %, B B * %, A A A %, B A B %, and a marker cell has no
    # column, joins no transition and is outside the entropy.
    from traminex.seqtree_explore import _features_from_sequences

    _, canonical = _exploratory_marker_pair("both")
    durations, duration_names = _features_from_sequences(canonical, kind="durations")
    assert duration_names == ["dur_A", "dur_B"]
    np.testing.assert_array_equal(durations, [[1, 1], [0, 2], [3, 0], [1, 2]])
    transitions, transition_names = _features_from_sequences(canonical, kind="transitions")
    assert transition_names == ["tr_A_A", "tr_A_B", "tr_B_A", "tr_B_B"]
    np.testing.assert_array_equal(
        transitions, [[0, 0, 0, 0], [0, 0, 0, 1], [2, 0, 0, 0], [0, 1, 1, 0]]
    )
    summary, summary_names = _features_from_sequences(canonical, kind="summary")
    assert summary_names == ["entropy_norm", "transition_count_norm"]
    # Entropy of the state shares over log(2): {A, B} 1, {B} 0, {A} 0, {A, B, B} H(1/3, 2/3);
    # changes between two adjacent state cells over seqlen - 1: 0, 0, 0, 2 (batch-12 review M-2).
    np.testing.assert_allclose(
        summary,
        [[1, 0], [0, 0], [0, 0], [0.9182958340544896, 2 / 3]],
        rtol=0,
        atol=1e-10,
    )


def test_exploratory_tree_ordinary_no_marker_feature_matrix_and_names():
    from traminex.seqtree_explore import _features_from_sequences

    ordinary = tx.seqdef([list("AAB"), list("BBA")], alphabet=["A", "B"])
    actual, names = _features_from_sequences(ordinary)
    assert names == [
        "dur_A",
        "dur_B",
        "tr_A_A",
        "tr_A_B",
        "tr_B_A",
        "tr_B_B",
        "entropy_norm",
        "transition_count_norm",
    ]
    # Durations 2/1 and 1/2; each row has one self-transition and one change;
    # entropy of (2/3,1/3) in bits, and one change over two adjacent pairs.
    np.testing.assert_allclose(
        actual,
        [[2, 1, 1, 1, 0, 0, 0.9182958340544896, 0.5], [1, 2, 0, 0, 1, 1, 0.9182958340544896, 0.5]],
        rtol=0,
        atol=1e-10,
    )


_NUMERIC_REAL = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 0], [1, 1, 0]], dtype=float)


@pytest.mark.parametrize("states", _LAYOUTS)
@pytest.mark.parametrize("with_missing", [False, True])
def test_seqasnum_factory_maps_real_alphabet_first_and_missing_last(states, with_missing):
    # TraMineR:::seqasnum compares labels against the real alphabet, then nr.
    factory = _factory(states)
    original_data = factory.data.copy()
    original_alphabet = list(factory.alphabet)
    missing = np.full((4, 1), 2.0 if with_missing else np.nan)
    expected = np.hstack([missing, _NUMERIC_REAL])
    np.testing.assert_allclose(
        tx.seqasnum(factory, with_missing=with_missing), expected, rtol=0, atol=0
    )
    np.testing.assert_array_equal(factory.data, original_data)
    assert factory.alphabet == original_alphabet
    assert factory.nr == "_"


@pytest.mark.parametrize("states", _LAYOUTS)
@pytest.mark.parametrize("with_missing", [False, True])
def test_seqnum_factory_recodes_labels_and_missing_metadata_with_numeric_roundtrip(
    states, with_missing
):
    factory = _factory(states)
    original_data = factory.data.copy()
    original_alphabet = list(factory.alphabet)
    original_missing_code = factory.nr_code
    numeric = tx.seqnum(factory, with_missing=with_missing)
    assert tx.alphabet(numeric) == ["0", "1"]
    assert numeric.nr == ("2" if with_missing else "_")
    assert numeric.void == "%"
    assert numeric.nr_code is not None
    assert numeric.alphabet[numeric.nr_code] == numeric.nr
    prefix = "2" if with_missing else "_"
    assert numeric.to_strings() == [
        f"{prefix}-0-0-1",
        f"{prefix}-1-0-1",
        f"{prefix}-0-1-0",
        f"{prefix}-1-1-0",
    ]
    expected = np.hstack([np.full((4, 1), 2.0 if with_missing else np.nan), _NUMERIC_REAL])
    np.testing.assert_allclose(
        tx.seqasnum(numeric, with_missing=with_missing), expected, rtol=0, atol=0
    )
    np.testing.assert_array_equal(factory.data, original_data)
    assert factory.alphabet == original_alphabet
    assert factory.nr == "_"
    assert factory.nr_code == original_missing_code


@pytest.mark.parametrize("states", _LAYOUTS)
@pytest.mark.parametrize("with_missing", [False, True])
def test_seqtransn_pweight_factory_matches_pinned_r_after_numeric_recoding(states, with_missing):
    # Global real-state rates are AA=BB=.25, AB=BA=.75. Each ordinary change
    # adds .25; a transition involving missing contributes NA, as in R.
    expected = [np.nan] * 4 if with_missing else [0.25, 0.5, 0.5, 0.25]
    np.testing.assert_allclose(
        tx.seqtransn(_factory(states), pweight=True, with_missing=with_missing),
        expected,
        rtol=0,
        atol=1e-14,
    )


@pytest.mark.parametrize("marker_kind", ["missing", "void", "both"])
@pytest.mark.parametrize("with_missing", [False, True])
def test_numeric_recoding_preserves_external_marker_twins_and_void_semantics(
    marker_kind, with_missing
):
    external, canonical = _exploratory_marker_pair(marker_kind)
    original_external = external.copy()
    if marker_kind == "both":
        expected = np.array(
            [[0, 2, 1, np.nan], [1, 1, 2, np.nan], [0, 0, 0, np.nan], [1, 0, 1, np.nan]]
        )
    elif marker_kind == "missing":
        expected = np.array([[0, 2, 1], [1, 1, 2], [0, 0, 0], [1, 0, 1]], dtype=float)
    else:
        expected = np.array(
            [[0, 0, 1, np.nan], [1, 1, 1, np.nan], [0, 0, 0, np.nan], [1, 0, 1, np.nan]]
        )
    if marker_kind != "void" and not with_missing:
        expected[expected == 2] = np.nan
    converted = []
    for source in [external, canonical]:
        np.testing.assert_allclose(tx.seqasnum(source, with_missing=with_missing), expected)
        numeric = tx.seqnum(source, with_missing=with_missing)
        assert tx.alphabet(numeric) == ["0", "1"]
        if marker_kind != "void":
            assert numeric.nr == ("2" if with_missing else "*")
        assert numeric.void == "%"
        np.testing.assert_allclose(tx.seqasnum(numeric, with_missing=with_missing), expected)
        converted.append(numeric.to_strings())
    assert converted[0] == converted[1]
    np.testing.assert_array_equal(external.data, original_external.data)
    assert external.alphabet == original_external.alphabet
    assert external.nr == original_external.nr
    assert external.void == original_external.void


def test_numeric_recoding_and_pweight_ordinary_no_marker_controls():
    ordinary = _factory(missing=False)
    np.testing.assert_array_equal(tx.seqasnum(ordinary), _NUMERIC_REAL)
    numeric = tx.seqnum(ordinary)
    assert tx.alphabet(numeric) == ["0", "1"]
    assert numeric.to_strings() == ["0-0-1", "1-0-1", "0-1-0", "1-1-0"]
    np.testing.assert_array_equal(tx.seqasnum(numeric), _NUMERIC_REAL)
    np.testing.assert_allclose(tx.seqtransn(ordinary, pweight=True), [0.25, 0.5, 0.5, 0.25])


def test_seqnum_with_missing_updates_nr_symbol_even_without_observed_missing():
    ordinary = tx.seqdef([list("AAB"), list("BBA")], alphabet=["A", "B"])
    assert ordinary.nr_code is None
    numeric = tx.seqnum(ordinary, with_missing=True)
    assert numeric.nr == "2"
    assert tx.alphabet(numeric) == ["0", "1"]
    assert numeric.to_strings() == ["0-0-1", "1-1-0"]
    assert ordinary.nr == "*"


@pytest.mark.parametrize("states", _LAYOUTS)
@pytest.mark.parametrize("include_missing", [False, True])
@pytest.mark.parametrize("function", ["seqintegr", "seqibad"])
def test_integrative_columns_and_badness_follow_real_then_missing_order(
    states, include_missing, function
):
    factory = _factory(states)
    original = factory.data.copy()
    expected_integration = (
        [[0.5, 0.4, 0.1], [0.3, 0.6, 0.1], [0.6, 0.3, 0.1], [0.4, 0.5, 0.1]]
        if include_missing
        else [[5 / 9, 4 / 9], [3 / 9, 6 / 9], [6 / 9, 3 / 9], [4 / 9, 5 / 9]]
    )
    expected_badness = [0.3, 0.4, 0.25, 0.35] if include_missing else [4 / 9, 6 / 9, 3 / 9, 5 / 9]
    expected = expected_integration if function == "seqintegr" else expected_badness
    np.testing.assert_allclose(
        getattr(tx, function)(factory, include_missing=include_missing),
        expected,
        rtol=0,
        atol=1e-14,
    )
    np.testing.assert_array_equal(factory.data, original)


@pytest.mark.parametrize("states", _LAYOUTS)
@pytest.mark.parametrize("include_missing", [False, True])
@pytest.mark.parametrize("function", ["seqprecarity", "seqinsecurity"])
def test_integrative_first_state_precedence_uses_label_rank_not_marker_position(
    states, include_missing, function
):
    # Use the default R nr='*' to isolate state ordering: R seqdss itself
    # rewrites custom nr symbols, which is a separate source behavior.
    star_states = None if states is None else ["*" if label == "_" else label for label in states]
    factory = _factory(star_states, fill="*")
    expected_precarity = [1, 1, 1, 1] if include_missing else [0, 1, 0, 1]
    expected_insecurity = (
        [0.8943108670863426, 1.0728281607546042, 1.0728281607546042, 0.8943108670863426]
        if include_missing
        else [0.5532617927812262, 0.9490969975607296, 0.782430330894063, 1.0532617927812262]
    )
    # correction=0 isolates the start-cost term from degradation estimation.
    options = (
        {"include_missing": include_missing, "otto": 1}
        if function == "seqprecarity"
        else {"with_missing": include_missing}
    )
    expected = expected_precarity if function == "seqprecarity" else expected_insecurity
    np.testing.assert_allclose(
        getattr(tx, function)(factory, correction=0, **options),
        expected,
        rtol=0,
        atol=1e-13,
    )


@pytest.mark.parametrize("with_void", [False, True])
@pytest.mark.parametrize("include_missing", [False, True])
@pytest.mark.parametrize("threshold", [False, True])
def test_main_tokens_selects_dominant_declared_external_missing_state(
    with_void, include_missing, threshold
):
    external_options = {"nr_code": 9}
    canonical_options = {"nr_code": 2}
    canonical_alphabet = ["a", "b", "*"]
    if with_void:
        external_options["void_code"] = 12
        canonical_options["void_code"] = 3
        canonical_alphabet.append("%")
    external = tx.STSList(np.array([[9, 9, 9, 0], [9, 9, 9, 1]]), ["a", "b"], **external_options)
    canonical = tx.STSList(
        np.array([[2, 2, 2, 0], [2, 2, 2, 1]]), canonical_alphabet, **canonical_options
    )
    assert external.to_strings() == canonical.to_strings()
    expected_top = [2] if include_missing else [0]
    expected_threshold = [2] if include_missing else []
    original = external.data.copy()
    options = {"k": 8, "mint": 1} if threshold else {"k": 1}
    expected = expected_threshold if threshold else expected_top
    for source in [external, canonical]:
        np.testing.assert_array_equal(
            tx.seqmaintokens(source, include_missing=include_missing, **options), expected
        )
    np.testing.assert_array_equal(external.data, original)
    assert external.nr_code == 9


@pytest.mark.parametrize("include_missing", [False, True])
def test_integrative_external_marker_twins_preserve_label_order(include_missing):
    external, canonical = _exploratory_marker_pair("both")
    for name, kwargs in [
        ("seqintegr", {"include_missing": include_missing}),
        ("seqibad", {"include_missing": include_missing}),
        ("seqprecarity", {"include_missing": include_missing, "correction": 0, "otto": 1}),
        ("seqinsecurity", {"with_missing": include_missing, "correction": 0}),
    ]:
        np.testing.assert_allclose(
            getattr(tx, name)(external, **kwargs),
            getattr(tx, name)(canonical, **kwargs),
            rtol=0,
            atol=1e-14,
        )


def test_integrative_ordinary_no_missing_controls_match_pinned_r():
    ordinary = _factory(missing=False)
    np.testing.assert_allclose(
        tx.seqintegr(ordinary),
        [[0.5, 0.5], [1 / 3, 2 / 3], [2 / 3, 1 / 3], [0.5, 0.5]],
        rtol=0,
        atol=1e-14,
    )
    np.testing.assert_allclose(tx.seqibad(ordinary), [0.5, 2 / 3, 1 / 3, 0.5], rtol=0, atol=1e-14)
    np.testing.assert_array_equal(tx.seqprecarity(ordinary, correction=0, otto=1), [0, 1, 0, 1])
    np.testing.assert_allclose(
        tx.seqinsecurity(ordinary, correction=0),
        [0.6776045432457228, 1.1249442016503943, 0.9582775349837277, 1.1776045432457227],
        rtol=0,
        atol=1e-13,
    )
