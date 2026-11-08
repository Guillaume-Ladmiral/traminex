import numpy as np
import pytest

from traminex import seqdef
from traminex.stats import entropy, seqindic, seqnum, seqtab, seqtabstocc, seqtransn, seqtrate


def test_seqtab_counts_and_props():
    mat = np.array(
        [
            ["A", "B", "C"],
            ["A", "B", "C"],
            ["B", "A", "C"],
        ]
    )
    x = seqdef(mat)
    tab = seqtab(x, idxs=0, format="STS", sep="-")
    # Expect two occurrences of A-B-C and one of B-A-C
    assert tab["count"].sum() == 3
    assert tab.iloc[0]["count"] >= tab.iloc[-1]["count"]  # sorted
    # Weighted variant
    w = np.array([1.0, 2.0, 1.0])
    tabw = seqtab(x, idxs=0, format="STS", weights=w)
    # Weighted total should be 4.0
    assert pytest.approx(tabw["count"].sum(), rel=1e-8) == 4.0
    assert pytest.approx(tabw.loc[tabw["sequence"] == "A-B-C", "count"].item(), rel=1e-8) == 3.0


def test_seqtransn_per_sequence_counts():
    # seqtransn is the per-sequence transition count = length(seqdss) - 1.
    mat = np.array(
        [
            [0, 1, 1, 2],  # dss 0,1,2 -> 2 transitions
            [1, 0, 1, 2],  # dss 1,0,1,2 -> 3 transitions
        ]
    )
    x = seqdef(mat)
    np.testing.assert_array_equal(seqtransn(x), np.array([2.0, 3.0]))
    # norm divides by (length - 1); both sequences have length 4.
    np.testing.assert_allclose(seqtransn(x, norm=True), np.array([2.0 / 3.0, 1.0]))


def test_transition_count_matrix_private():
    from traminex.stats_basic import _transition_count_matrix

    mat = np.array([[0, 1, 1, 2], [1, 0, 1, 2]])
    x = seqdef(mat)
    M = _transition_count_matrix(x)
    K = len(x.alphabet)
    assert M.shape == (K, K)
    assert M[0, 1] == 2 and M[1, 1] == 1 and M[1, 2] == 2 and M[1, 0] == 1
    P = _transition_count_matrix(x, normalize=True)
    for i in range(K):
        if M[i].sum() > 0:
            assert pytest.approx(P[i].sum(), rel=1e-8) == 1.0


def test_entropy_normalized_uniform():
    # 3 states uniformly distributed across all positions
    mat = np.array(
        [
            [0, 1, 2],
            [2, 0, 1],
            [1, 2, 0],
        ]
    )
    x = seqdef(mat)
    Hn = entropy(x, normalized=True)
    assert 0.99 <= Hn <= 1.01


def test_seqtransn_short_sequences_zero():
    mat = np.array([["A"], ["B"]], dtype=object)
    x = seqdef(mat)
    # Single-position sequences have a one-state DSS -> 0 transitions.
    np.testing.assert_array_equal(seqtransn(x), np.array([0.0, 0.0]))


def test_seqtabstocc_include_missing_flag():
    mat = np.array([["A", "*"], ["A", "B"]], dtype=object)
    x = seqdef(mat, nr="*")

    out_ex = seqtabstocc(x, include_missing=False)
    out_in = seqtabstocc(x, include_missing=True)

    assert all("*" not in p for p in out_ex["pattern"])
    assert any("*" in p for p in out_in["pattern"])


def test_seqtabstocc_accepts_r_equivalent_with_missing_kwarg():
    mat = np.array([["A", "*", "B", "*"]], dtype=object)
    x = seqdef(mat, missing="*")

    result = seqtabstocc(x, with_missing=True)
    assert any("*" in pattern for pattern in result["pattern"])


def test_seqindic_default_indicator_table():
    """seqindic returns the composite indicator table with R column labels."""
    mat = np.array([["A", "A", "B", "C"], ["A", "B", "B", "C"], ["C", "C", "B", "A"]])
    x = seqdef(mat)
    df = seqindic(x)  # default indic = visited/trans/entr/cplx/turb2n
    assert list(df.columns) == ["Visited", "Trans", "Entr", "Cplx", "Turb2n"]
    assert len(df) == x.nseq
    # all three sequences visit {A,B,C} with two transitions
    assert df["Visited"].tolist() == [3.0, 3.0, 3.0]
    assert df["Trans"].tolist() == [2.0, 2.0, 2.0]


def test_seqindic_rejects_invalid_indicator():
    x = seqdef(np.array([["A", "B"], ["B", "A"]]))
    with pytest.raises(ValueError, match="invalid values in indic"):
        seqindic(x, indic=["not-an-indicator"])


def test_seqnum_seqasnum_numeric_recoding():
    """seqasnum is the numeric coding; seqnum recodes the alphabet to 0..K-1."""
    from traminex.stats import seqasnum

    mat = np.array([["A", "*", "B"], ["*", "A", "*"]], dtype=object)
    x = seqdef(mat, alphabet=["A", "B"], nr="*")  # A=0, B=1, nr=2 (no void)

    a = seqasnum(x)  # nr -> NaN by default
    assert a[0, 0] == 0.0 and a[0, 2] == 1.0 and np.isnan(a[0, 1])
    am = seqasnum(x, with_missing=True)  # nr keeps its code (2)
    assert am[0, 1] == 2.0 and am[1, 0] == 2.0

    num = seqnum(x)
    assert num.alphabet[0] == "0" and num.alphabet[1] == "1"
    assert num.alphabet[2] == "*"  # nr symbol preserved (with_missing=False)


def test_seqtrate_include_missing_flag():
    mat = np.array([["A", "*", "B"], ["A", "B", "B"]], dtype=object)
    x = seqdef(mat, nr="*")

    out_ex = seqtrate(x, include_missing=False, as_dataframe=True)
    out_in = seqtrate(x, include_missing=True, as_dataframe=True)

    assert "*" not in out_ex.index
    assert "*" not in out_ex.columns
    assert "*" in out_in.index
    assert "*" in out_in.columns


@pytest.mark.parametrize("dtype", [np.int8, np.uint8, np.int32, np.int64])
@pytest.mark.parametrize("normalize", [False, True])
def test_transition_matrix_widens_index_arithmetic_before_multiplication(dtype, normalize):
    from traminex import STSList
    from traminex.stats_basic import _transition_count_matrix

    data = np.array([[19, 19], [19, 0]], dtype=dtype)
    data.setflags(write=False)
    source = STSList(data, [f"S{i}" for i in range(20)], weights=np.array([0.25, 1.5]))
    expected = np.zeros((20, 20))
    expected[19, 19] = 1.0 / 7.0 if normalize else 0.25
    expected[19, 0] = 6.0 / 7.0 if normalize else 1.5
    np.testing.assert_array_equal(_transition_count_matrix(source, normalize=normalize), expected)
    np.testing.assert_array_equal(source.data, data)


@pytest.mark.parametrize("dtype", [np.int8, np.uint8])
def test_transition_heatmap_draws_narrow_state_counts(dtype):
    plt = pytest.importorskip("matplotlib.pyplot")
    from traminex import STSList, transition_heatmap

    source = STSList(np.array([[19, 19], [19, 0]], dtype=dtype), [f"S{i}" for i in range(20)])
    expected = np.zeros((20, 20))
    expected[19, 19] = expected[19, 0] = 1.0
    fig, ax = transition_heatmap(source, normalize=False)
    try:
        np.testing.assert_array_equal(ax.images[0].get_array(), expected)
        assert [label.get_text() for label in ax.get_xticklabels()] == source.alphabet
    finally:
        plt.close(fig)
