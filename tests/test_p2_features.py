import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.cluster import seqrep
from traminex.diss import dissdomassoc
from traminex.sequences import SeqMD, seqMD_distances


def _make_sts_small():
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "B", "B"],
            "t2": ["B", "A", "B", "A"],
            "t3": ["B", "B", "A", "A"],
            "t4": ["C", "B", "A", "B"],
        }
    )
    return tx.seqdef(df)


def test_chi2_windowed_basic_properties():
    x = _make_sts_small()
    D = tx.seqdist(x, method="CHI2", step=2, overlap=False)
    assert D.shape == (x.nseq, x.nseq)
    assert np.allclose(np.diag(D), 0.0)
    assert np.all(D >= 0)
    assert np.allclose(D, D.T)
    # identical rows zero
    X2 = tx.seqdef(pd.DataFrame({"t1": ["A", "A"], "t2": ["B", "B"], "t3": ["C", "C"]}))
    D2 = tx.seqdist(X2, method="CHI2", step=2, overlap=True)
    assert D2[0, 1] == pytest.approx(0.0, abs=1e-12)


def test_seqdiff_shapes_bounds_and_no_effect_when_groups_identical():
    # RP4 T16: seqdiff is the R-faithful sliding-window dissassoc -> a dict of
    # {stat, discrepancy} tables, not the old per-position residuals.
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "B", "B"],
            "t2": ["B", "B", "A", "A"],
            "t3": ["A", "A", "A", "A"],
        }
    )
    X = tx.seqdef(df)
    group = np.array([0, 1, 0, 1])
    from traminex import seqdiff as _seqdiff

    out = _seqdiff(X, group)
    assert set(out.keys()) == {"stat", "discrepancy"}
    assert list(out["stat"].columns) == [
        "Pseudo F",
        "Pseudo Fbf",
        "Pseudo R2",
        "Bartlett",
        "Levene",
    ]
    assert list(out["discrepancy"].columns) == ["0", "1", "Total"]
    # cmprange (0, 1) over 3 columns -> 2 window positions.
    assert out["stat"].shape[0] == 2
    # Both groups hold the same sequences (rows 0==1, 2==3) -> no between-group
    # variance -> Pseudo R2 == 0.
    np.testing.assert_allclose(out["stat"]["Pseudo R2"].to_numpy(), 0.0, atol=1e-12)


def test_dissdomassoc_basic():
    # RP4 T14: R-faithful inter-domain association -- a weighted Pearson
    # correlation matrix (+ optional joint column), not the old variance summary.
    D1 = np.array(
        [[0.0, 1.0, 2.0, 1.5], [1.0, 0.0, 3.0, 2.5], [2.0, 3.0, 0.0, 1.0], [1.5, 2.5, 1.0, 0.0]],
        float,
    )
    D2 = np.array(
        [[0.0, 0.5, 1.5, 2.0], [0.5, 0.0, 1.0, 0.5], [1.5, 1.0, 0.0, 2.5], [2.0, 0.5, 2.5, 0.0]],
        float,
    )
    out = dissdomassoc([D1, D2], what=("pearson",))
    corr = out["Pearson"]
    assert corr.shape == (2, 2)
    np.testing.assert_allclose(np.diag(corr), 1.0)
    np.testing.assert_allclose(corr, corr.T)
    assert "p.pearson" in out

    # jointdiss appends an extra column/row to the correlation matrix.
    out_j = dissdomassoc([D1, D2], jointdiss=D1 + D2, what=("pearson",))
    assert out_j["Pearson"].shape == (3, 3)


def test_seqMD_distances_additivity_and_errors():
    # Two channels with 3 sequences each
    A = tx.seqdef(pd.DataFrame({"t1": ["A", "A", "B"], "t2": ["B", "A", "B"]}))
    B = tx.seqdef(pd.DataFrame({"t1": ["X", "Y", "X"], "t2": ["Y", "Y", "X"]}))
    obj = SeqMD([A, B], names=["A", "B"])
    D1 = tx.seqdist(A, method="HAM")
    D2 = tx.seqdist(B, method="HAM")
    w = [0.3, 1.7]
    Dag = seqMD_distances(obj, method="HAM", weights=w)
    assert np.allclose(Dag, w[0] * D1 + w[1] * D2)
    # Weight length mismatch raises
    with pytest.raises(ValueError):
        _ = seqMD_distances(obj, method="HAM", weights=[1.0])


def test_seqMD_distances_refuses_an_empty_channel_list():
    """cycle-6 F-010: the declared `-> np.ndarray` must not be a None.

    A zero-channel SeqMD used to fall out of the summation loop with `Dsum`
    still None, so the caller met an unrelated
    `AttributeError: 'NoneType' object has no attribute 'shape'` at whatever it
    did next, and a `# type: ignore[return-value]` kept mypy quiet about it.
    """
    with pytest.raises(ValueError, match="channels"):
        seqMD_distances(SeqMD(channels=[]), method="HAM")


def test_seqMD_distances_refusal_names_the_argument_not_a_downstream_symptom():
    with pytest.raises(ValueError) as excinfo:
        seqMD_distances(SeqMD(channels=[]), method="HAM")
    message = str(excinfo.value)
    assert "channels" in message
    assert "NoneType" not in message


def test_seqmd_constructor_refuses_empty_and_unaligned_channels():
    two_rows = tx.seqdef([["A", "B"], ["B", "A"]])
    one_row = tx.seqdef([["A", "B"]])

    with pytest.raises(ValueError, match="channels must contain at least one channel"):
        SeqMD([])
    with pytest.raises(ValueError, match="channels must contain the same number of sequences"):
        SeqMD([two_rows, one_row])


def test_seqmd_constructor_refuses_non_sequence_channels():
    sequences = tx.seqdef([["A", "B"], ["B", "A"]])

    with pytest.raises(TypeError, match="channels must contain only STSList"):
        SeqMD([sequences, np.array([[0, 1]])])  # type: ignore[list-item]


@pytest.mark.parametrize(
    ("weights", "error", "message"),
    [
        (1.0, ValueError, "one-dimensional"),
        (np.array([[1.0, 1.0]]), ValueError, "one-dimensional"),
        ([1.0], ValueError, "2 elements"),
        ([1.0, "heavy"], TypeError, "numeric"),
        ([1.0, "2.0"], TypeError, "numeric"),
        ([True, 1.0], TypeError, "numeric"),
        ([1.0, -0.1], ValueError, "non-negative"),
        ([1.0, np.nan], ValueError, "finite"),
        ([1.0, np.inf], ValueError, "finite"),
        ([1.0, -np.inf], ValueError, "finite"),
        ([0.0, 0.0], ValueError, "positive sum"),
    ],
)
def test_seqmd_distances_rejects_malformed_analytical_weights(weights, error, message):
    sequences = tx.seqdef([["A", "B"], ["B", "A"]])
    obj = SeqMD([sequences, sequences.copy()])

    with pytest.raises(error, match=message):
        seqMD_distances(obj, method="HAM", weights=weights)


def test_seqmd_distances_accepts_zero_component_with_positive_total():
    first = tx.seqdef([["A", "B"], ["B", "A"]])
    second = tx.seqdef([["X", "X"], ["Y", "Y"]])
    obj = SeqMD([first, second])

    actual = seqMD_distances(obj, method="HAM", weights=[1.5, 0.0])

    np.testing.assert_allclose(actual, 1.5 * tx.seqdist(first, method="HAM"))


def test_seqmd_distances_revalidates_channels_after_mutation():
    two_rows = tx.seqdef([["A", "B"], ["B", "A"]])
    obj = SeqMD([two_rows, two_rows.copy()])
    obj.channels[1] = tx.seqdef([["A", "B"]])

    with pytest.raises(ValueError, match="channels must contain the same number of sequences"):
        seqMD_distances(obj, method="HAM")


def test_seqrep_frequency_and_distance_modes():
    # Build dataset with repeated sequences
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "A", "B", "B"],
            "t2": ["B", "B", "B", "A", "A"],
            "t3": ["A", "A", "A", "A", "A"],
        }
    )
    X = tx.seqdef(df)

    # RP4 T13c: seqrep wraps dissrep -> rep_indices + statistics/quality, no
    # coverage key. Density criterion, diss computed internally.
    out = seqrep(X, criterion="density", coverage=0.5, method="OM", indel=1.0, subm=2.0)
    reps = out["rep_indices"]
    assert reps.ndim == 1 and reps.size >= 1
    assert np.all((reps >= 0) & (reps < X.nseq))
    assert "statistics" in out and "quality" in out
    assert out["representatives"].shape[0] == reps.size

    # Distance-based selection of a fixed number of representatives.
    D = tx.seqdist(X, method="OM", indel=1.0, subm=2.0)
    out2 = seqrep(X, diss=D, nrep=2, criterion="dist")
    reps2 = out2["rep_indices"]
    assert reps2.size == 2
    assert len(set(reps2.tolist())) == 2
