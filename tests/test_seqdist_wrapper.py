import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import r_anchor_values

import traminex as tx
from traminex.distances import pairwise, pairwise_between_sets


def make_simple_sts():
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "B"],
            "t2": ["B", "A", "B"],
            "t3": ["B", "B", "A"],
        }
    )
    x = tx.seqdef(df)
    return x


def test_seqdist_refseq_zero_breaks_frequency_ties_by_r_numeric_sps_order():
    x = tx.seqdef(
        np.asarray([["B", "B"], ["A", "A"], ["C", "C"]], dtype=object),
        alphabet=["A", "B", "C"],
    )

    actual = tx.seqdist(x, method="LCS", norm="none", refseq=0)

    np.testing.assert_array_equal(actual, np.asarray([4.0, 0.0, 4.0]))


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Empty\ sequence\(s\)\ \ \d+,\ \d+\ \ may\ cause\ inconsistent\ results!:UserWarning"
)
def test_seqdist_refseq_zero_keeps_structurally_empty_frequency_rows():
    x = tx.seqdef(
        np.asarray(
            [
                ["%", "%", "%"],
                ["%", "%", "%"],
                ["A", "*", "A"],
                ["A", "B", "%"],
            ],
            dtype=object,
        ),
        alphabet=["A", "B"],
        nr="*",
    )

    actual = tx.seqdist(x, method="LCS", with_missing=True, norm="none", refseq=0)

    np.testing.assert_array_equal(actual, np.asarray([0.0, 0.0, 3.0, 2.0]))


def test_seqdist_matches_pairwise_om():
    x = make_simple_sts()
    D1 = tx.seqdist(x, method="OM", indel=1.0, subm=2.0, norm="length")
    D2 = pairwise(x.data, method="OM", indel=1.0, subm=2.0, norm="length")
    assert D1.shape == (x.nseq, x.nseq)
    assert np.allclose(D1, D2)


def test_seqdist_om_method_default_matches_explicit_om():
    x = make_simple_sts()

    default = tx.seqdist(x, indel=1.0, subm=2.0, norm="length")
    explicit = tx.seqdist(x, method="OM", indel=1.0, subm=2.0, norm="length")

    np.testing.assert_array_equal(default, explicit)


@pytest.mark.parametrize(
    ("method", "method_options"),
    [
        ("OM", {}),
        ("OMLOC", {}),
        ("OMSLEN", {}),
        ("OMSPELL", {}),
        ("OMSTRAN", {"otto": 0.5}),
        ("HAM", {}),
        ("LCS", {}),
        ("TWED", {"nu": 0.5, "h": 1.0}),
    ],
)
@pytest.mark.parametrize("cost_method", ["TRATE", "CONSTANT", "INDELS", "INDELSLOG"])
def test_seqdist_resolves_every_r_cost_method_string(method, method_options, cost_method):
    x = make_simple_sts()
    cval = 1.0 if method == "HAM" and cost_method == "CONSTANT" else 2.0
    costs = tx.seqcost(x, method=cost_method, cval=cval, miss_cost=cval)
    explicit_options = dict(method_options)
    if method != "LCS":
        explicit_options["subm"] = np.asarray(costs["sm"], dtype=float)
    if method not in {"OMLOC", "HAM", "LCS", "TWED"}:
        explicit_options["indel"] = costs["indel"]

    expected = tx.seqdist(x, method=method, **explicit_options)
    actual = tx.seqdist(x, method=method, subm=cost_method, **method_options)

    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)


def test_seqdist_resolves_cost_method_string_for_raw_integer_matrix():
    sequences = np.asarray([[0, 1, 1], [0, 0, 1]], dtype=int)

    actual = tx.seqdist(sequences, method="OM", subm="CONSTANT")

    assert actual[1, 0] == pytest.approx(2.0)


@pytest.mark.parametrize("method", ["OM", "OMLOC", "OMSLEN", "OMSPELL", "HAM", "LCS", "TWED"])
def test_seqdist_refuses_non_source_active_cost_method_strings(method):
    options = {"nu": 0.5, "h": 1.0} if method == "TWED" else {}
    with pytest.raises(ValueError, match="TRATE.*CONSTANT.*INDELS.*INDELSLOG"):
        tx.seqdist(make_simple_sts(), method=method, subm="FUTURE", **options)


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        ("OM", 4.05),
        ("OMLOC", 4.05),
        ("OMSPELL", 1.0),
        ("OMSLEN", 2.863782464),
        ("HAM", 4.05),
    ],
)
def test_seqdist_accepts_unreachable_infinite_indels_cost_cells(method, expected):
    x = tx.seqdef(
        [["A", "A", "B"], ["A", "B", "B"], ["B", "B", "A"]],
        alphabet=["A", "B", "C", "D"],
    )

    actual = tx.seqdist(x, method=method, subm="INDELS")

    assert actual[1, 0] == pytest.approx(expected, abs=1e-9)


def test_seqdist_omstran_preserves_r_failure_for_generated_nonfinite_costs():
    x = tx.seqdef(
        [["A", "A", "B"], ["A", "B", "B"], ["B", "B", "A"]],
        alphabet=["A", "B", "C", "D"],
    )

    # R stops too (in checkcost); the refusal names the non-finite substitution costs, not the
    # indel (batch-9 review I-1, TX-AUD-20261019-079).
    with pytest.raises(ValueError, match=r"substitution costs\b.*not finite"):
        tx.seqdist(x, method="OMSTRAN", subm="INDELS", otto=0.5)


def test_seqdist_twed_unobserved_states_change_no_distance():
    # R gives Inf here: the INDELS stand-in of C and D, which no sequence contains, enters
    # TWED's costs. traminex corrects it (TX-AUD-20261019-002, seqdist DIVERGENCE row): the
    # distances equal those without the unobserved states; tests/test_indels_unobserved_state_norm.py
    # compares that rule with R's values on the restricted alphabet.
    rows = [["A", "A", "B"], ["A", "B", "B"], ["B", "B", "A"]]
    x = tx.seqdef(rows, alphabet=["A", "B", "C", "D"])
    observed_only = tx.seqdef(rows, alphabet=["A", "B"])

    actual = tx.seqdist(x, method="TWED", subm="INDELS", nu=0.5, h=1.0)
    expected = tx.seqdist(observed_only, method="TWED", subm="INDELS", nu=0.5, h=1.0)

    assert np.all(np.isfinite(actual))
    np.testing.assert_allclose(actual, expected, rtol=1e-12)


@pytest.mark.parametrize("method", ["OM", "HAM"])
def test_pairwise_accepts_nonfinite_costs_only_for_unobserved_states(method):
    sequences = np.asarray([[0, 0, 1], [0, 1, 1]], dtype=int)
    costs = np.asarray(
        [
            [0.0, 4.0, np.inf],
            [4.0, 0.0, np.inf],
            [np.inf, np.inf, 0.0],
        ]
    )
    options = {"subm": costs}
    if method == "OM":
        options["indel"] = 2.0

    actual = pairwise(sequences, method=method, **options)

    assert actual[1, 0] == pytest.approx(4.0)


def test_negative_scalar_substitution_cost_names_the_supplied_argument():
    with pytest.raises(ValueError, match="OM substitution costs cannot be negative"):
        tx.seqdist(make_simple_sts(), method="OM", subm=-1.0)


def test_seqdist_hamming_shape():
    x = make_simple_sts()
    D = tx.seqdist(x, method="HAM")
    assert D.shape == (x.nseq, x.nseq)
    assert np.all(D >= 0)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ Empty\ sequence\(s\)\ \ \d+,\ \d+,\ \d+\ \ may\ cause\ inconsistent\ results!:UserWarning"
)
def test_seqdist_hamming_zero_length_sequences_returns_zero():
    x = tx.seqdef(np.empty((3, 0), dtype=int))
    D = tx.seqdist(x, method="HAM")
    assert D.shape == (3, 3)
    assert np.all(np.isfinite(D))
    assert np.allclose(D, 0.0)

    D_pairwise = pairwise(x.data, method="HAM")
    assert D_pairwise.shape == (3, 3)
    assert np.all(np.isfinite(D_pairwise))
    assert np.allclose(D_pairwise, 0.0)


@pytest.mark.parametrize("indel", [-1.0, 0.0, np.nan, -np.inf])
def test_seqdist_om_rejects_invalid_indel(indel):
    sequences = np.array([[0, 0], [1, 1]], dtype=int)

    with pytest.raises(ValueError, match="indel cost must be strictly positive"):
        tx.seqdist(sequences, method="OM", subm=2.0, indel=indel)


@pytest.mark.parametrize("indel", [-1.0, 0.0, np.nan, -np.inf])
def test_public_pairwise_om_rejects_invalid_indel(indel):
    sequences = np.array([[0, 0], [1, 1]], dtype=int)

    with pytest.raises(ValueError, match="indel cost must be strictly positive"):
        pairwise(sequences, method="OM", indel=indel)


@pytest.mark.parametrize(
    ("method", "subm"),
    [
        ("OM", -1.0),
        ("OM", np.array([[0.0, -1.0], [-1.0, 0.0]])),
        ("OM", np.array([[1.0, 2.0], [2.0, 0.0]])),
        ("OM", np.ones((2, 3))),
        ("HAM", np.array([[0.0, -1.0], [-1.0, 0.0]])),
        ("HAM", np.array([[1.0, 2.0], [2.0, 0.0]])),
        ("HAM", np.ones((2, 3))),
    ],
)
def test_public_distance_dispatchers_reject_invalid_substitution_costs(method, subm):
    sequences = np.array([[0, 0], [1, 1]], dtype=int)
    options = {"subm": subm}
    if method == "OM":
        options["indel"] = 1.0

    with pytest.raises(
        ValueError, match="substitution costs|substitution-cost|subm must be square"
    ):
        pairwise(sequences, method=method, **options)
    with pytest.raises(
        ValueError, match="substitution costs|substitution-cost|subm must be square"
    ):
        tx.seqdist(sequences, method=method, **options)


@pytest.mark.parametrize(
    "indel_by_pos",
    [
        np.array([1.0, -1.0]),
        np.array([1.0, np.nan]),
        np.array([1.0]),
        np.ones((2, 1)),
    ],
)
def test_om_rejects_invalid_position_indel_costs(indel_by_pos):
    sequences = np.array([[0, 0], [1, 1]], dtype=int)

    with pytest.raises(ValueError, match="indel_by_pos"):
        pairwise(
            sequences,
            method="OM",
            indel=1.0,
            subm=2.0,
            indel_by_pos=indel_by_pos,
        )


@pytest.mark.parametrize(
    "subm_by_pos",
    [
        np.array([[[0.0, -1.0], [-1.0, 0.0]]] * 2),
        np.array([[[0.0, np.nan], [1.0, 0.0]]] * 2),
        np.array([[[1.0, 2.0], [2.0, 0.0]]] * 2),
        np.zeros((1, 2, 2)),
        np.zeros((2, 2, 3)),
    ],
)
def test_om_rejects_invalid_position_substitution_costs(subm_by_pos):
    sequences = np.array([[0, 0], [1, 1]], dtype=int)

    with pytest.raises(ValueError, match="subm_by_pos"):
        pairwise(
            sequences,
            method="OM",
            indel=1.0,
            subm=2.0,
            subm_by_pos=subm_by_pos,
        )


@pytest.mark.parametrize("method", ["OM", "DHD", "CHI2"])
@pytest.mark.parametrize(
    ("bad", "message"),
    [
        (np.array([0, 1]), "2D"),
        (np.array([[0.0, 1.0], [1.0, 0.0]]), "integer state codes"),
        (np.array([[False, True], [True, False]]), "integer state codes"),
        (np.array([[0, -1], [1, 0]]), "non-negative"),
    ],
)
def test_seqdist_rejects_invalid_dense_state_matrices(method, bad, message):
    with pytest.raises((TypeError, ValueError), match=message):
        tx.seqdist(bad, method=method)


@pytest.mark.parametrize(
    ("bad", "message"),
    [
        (np.array([0, 1]), "2D"),
        (np.array([[0.0, 1.0], [1.0, 0.0]]), "integer state codes"),
        (np.array([[False, True], [True, False]]), "integer state codes"),
        (np.array([[0, -1], [1, 0]]), "non-negative"),
    ],
)
def test_pairwise_rejects_invalid_dense_state_matrices(bad, message):
    with pytest.raises((TypeError, ValueError), match=message):
        pairwise(bad, method="OM")


def test_seqdist_rejects_refseq_with_ref_sets():
    x = make_simple_sts()
    with pytest.raises(ValueError, match="refseq and ref_sets"):
        tx.seqdist(x, method="LCS", refseq=1, ref_sets=([0], [1]))


def test_seqdist_external_reference_refusal_contracts():
    x = make_simple_sts()
    two_rows = tx.seqdef([["A", "B"], ["B", "A"]], alphabet=x.alphabet)
    different_alphabet = tx.seqdef([["A", "C"]], alphabet=["A", "C"])
    too_long = tx.seqdef([["A", "B", "A", "B"]], alphabet=x.alphabet)

    with pytest.raises(ValueError, match="exactly one sequence"):
        tx.seqdist(x, method="LCS", refseq=two_rows)
    with pytest.raises(ValueError, match="same alphabet"):
        tx.seqdist(x, method="LCS", refseq=different_alphabet)
    with pytest.raises(ValueError, match="cannot be longer"):
        tx.seqdist(x, method="LCS", refseq=too_long)


def test_seqdist_refseq_sets_require_integral_1_based_indices():
    x = make_simple_sts()
    with pytest.raises(ValueError, match="positive integer"):
        tx.seqdist(x, method="LCS", refseq=([1.5], [2]))


def test_pairwise_between_sets_accepts_canonical_ragged_rows():
    left = [np.array([0, 1, 2]), np.array([0, 1])]
    right = [np.array([0, 2]), np.array([2])]
    actual = pairwise_between_sets(left, right, method="LCS")
    expected = np.array([[1.0, 2.0], [2.0, 3.0]])
    np.testing.assert_array_equal(actual, expected)


def test_pairwise_between_sets_deduplicates_canonical_ragged_rows():
    left = [np.array([0, 1]), np.array([0, 1]), np.array([2])]
    right = [np.array([2]), np.array([0, 1])]

    deduplicated = pairwise_between_sets(left, right, method="LCS", deduplicate=True)
    baseline = pairwise_between_sets(left, right, method="LCS", deduplicate=False)

    np.testing.assert_array_equal(deduplicated, baseline)


@pytest.mark.parametrize(
    ("method", "kwargs"),
    [
        ("CHI2", {}),
        ("DHD", {}),
        (
            "OMSTRAN",
            {"otto": 0.5, "transindel": "prob", "subm": 2.0, "indel": 1.0},
        ),
        ("OM", {"subm": "TRATE"}),
    ],
)
def test_deduplicate_preserves_dataset_fitted_distances(method, kwargs):
    rows = np.asarray(
        [[0, 0, 1], [0, 0, 1], [0, 1, 1], [1, 1, 0]],
        dtype=np.int32,
    )

    baseline = pairwise(rows, method=method, deduplicate=False, **kwargs)
    deduplicated = pairwise(rows, method=method, deduplicate=True, **kwargs)
    # The cross-block route for a dataset-fitted method is ``pairwise(ref_sets=...)``:
    # it carries the whole corpus, so the fit is the corpus fit.  A standalone
    # ``pairwise_between_sets(A, B)`` cannot know that corpus and refuses (F-001),
    # so the deduplication-invariance contract is asserted on the route that
    # supports it rather than dropped.
    sets = (np.array([0, 1]), np.array([2, 3]))
    between_baseline = pairwise(rows, method=method, ref_sets=sets, deduplicate=False, **kwargs)
    between_deduplicated = pairwise(rows, method=method, ref_sets=sets, deduplicate=True, **kwargs)

    np.testing.assert_allclose(deduplicated, baseline, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(between_deduplicated, between_baseline, rtol=0.0, atol=1e-12)
    # ...and the standalone two-set call refuses rather than silently re-fitting.
    with pytest.raises(ValueError, match="estimated from the whole corpus"):
        pairwise_between_sets(rows[:2], rows[2:], method=method, **kwargs)


def test_seqdist_accepts_generic_data_wrapper_without_alphabet_metadata():
    class DataWrapper:
        def __init__(self, data):
            self.data = data

    data = np.asarray([[0, 1, 1], [0, 0, 1]], dtype=np.int32)

    actual = tx.seqdist(DataWrapper(data), method="LCS", norm="none")
    expected = tx.seqdist(data, method="LCS", norm="none")

    np.testing.assert_array_equal(actual, expected)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ seqdist.\ 'with\.missing'\ set\ as\ FALSE\ as\ 'seqdata'\ has\ no\ non\-void\ missing\ values:UserWarning"
)
def test_dhd_ignores_declared_missing_state_when_no_row_contains_it():
    class DeclaredMissingWrapper:
        def __init__(self) -> None:
            self.data = np.asarray([[0, 0, 1], [0, 1, 1]], dtype=np.int32)
            self.alphabet = ["A", "B", "*"]
            self.nr = "*"
            self.void = "%"
            self.nr_code = 2
            self.void_code = None
            self.weights = None

    costs = np.broadcast_to(
        np.asarray([[0.0, 2.0], [2.0, 0.0]]),
        (3, 2, 2),
    ).copy()
    actual = tx.seqdist(
        DeclaredMissingWrapper(),
        method="DHD",
        with_missing=True,
        norm="none",
        subm=costs,
    )
    np.testing.assert_allclose(actual, np.asarray([[0.0, 2.0], [2.0, 0.0]]))


def test_external_refseq_missing_requires_opt_in_even_when_source_is_complete():
    source = tx.seqdef(
        np.array([["A", "A"], ["A", "B"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )
    reference = tx.seqdef(
        np.array([["A", "*"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )
    with pytest.raises(ValueError, match="with_missing must be true"):
        tx.seqdist(source, method="LCS", refseq=reference, with_missing=False)


def test_external_refseq_can_introduce_missing_state():
    """TraMineR accepts a missing reference for a complete compatible source."""
    source = tx.seqdef(
        np.array([["A", "A", "A"], ["A", "B", "A"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )
    reference = tx.seqdef(
        np.array([["A", "*", "A"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )

    actual = tx.seqdist(
        source,
        method="NMS",
        refseq=reference,
        with_missing=True,
    )

    # TraMineR 2.2-14: c(sqrt(10), sqrt(8)), recorded by
    # scripts/r_goldens/make_goldens_audit20260922_s10.R.
    r_values = r_anchor_values("nms_refseq_with_missing")
    np.testing.assert_allclose(r_values, [np.sqrt(10.0), np.sqrt(8.0)], rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(actual, r_values, rtol=0.0, atol=1e-10)


def test_dhd_generated_costs_safely_support_ref_only_missing_state():
    """Avoid TraMineR 2.2-14's undersized-cost undefined-memory defect."""
    source = tx.seqdef(
        np.array([["A", "A", "A"], ["A", "B", "A"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )
    reference = tx.seqdef(
        np.array([["A", "*", "A"]], dtype=object),
        alphabet=["A", "B"],
        nr="*",
    )

    actual = tx.seqdist(
        source,
        method="DHD",
        refseq=reference,
        with_missing=True,
        norm="none",
    )

    np.testing.assert_allclose(actual, [4.0, 4.0], rtol=0.0, atol=1e-12)
    assert np.isfinite(actual).all()


@pytest.mark.parametrize("method", ["LCS", "OM"])
def test_generic_alphabet_wrapper_without_marker_attributes(method):
    class DataWrapper:
        def __init__(self, data):
            self.data = data
            self.alphabet = ["*", "%"]

    data = np.asarray([[0, 1], [1, 0]], dtype=np.int32)
    kwargs = {"subm": np.array([[0.0, 2.0], [2.0, 0.0]]), "indel": 1.0} if method == "OM" else {}
    actual = tx.seqdist(DataWrapper(data), method=method, norm="none", **kwargs)
    np.testing.assert_array_equal(actual, [[0.0, 2.0], [2.0, 0.0]])
