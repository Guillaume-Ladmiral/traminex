import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex import seqformat as sf
from traminex.events import (
    events_to_state_sequences,
    events_transition_matrix,
    frequent_event_pairs,
    frequent_event_pairs_timegap,
    seqecreate,
)
from traminex.stats import seqlogp, seqtabstocc

matplotlib.use("Agg")


def _decode_states(x):
    return np.asarray(x.alphabet, dtype=object)[x.data]


def test_seqdef_from_long_keeps_a_star_state_and_voids_trailing_missing():
    # A state value spelled like the marker is a real state in the long builders
    # (TX-AUD-20261019-005; seqdef's STS route alone reads such a spelling as the
    # marker). Missing data in long input is a missing state or an absent row.
    long = pd.DataFrame(
        {
            "id": ["p1", "p1", "p2", "p2"],
            "time": [1, 2, 1, 2],
            "state": ["A", "*", "A", None],
        }
    )

    x = tx.seqdef_from_long(long, nr="*")

    assert tx.alphabet(x) == ["*", "A"]
    # p2's missing state is the last cell: void by seqdef's default right = "DEL",
    # so no cell needs the nr marker and the literal "*" keeps its own code.
    assert x.nr_code is None
    assert x.void_code is not None
    assert np.asarray(tx.seqlength(x, with_missing=False)).ravel().tolist() == [2, 1]
    out = seqtabstocc(x, include_missing=False)
    assert any("*" in pattern for pattern in out["pattern"])


def test_hspell_to_sts_propagates_fill_missing_code_metadata():
    df = pd.DataFrame(
        {
            "id": [1, 1, 2, 2],
            "start": [1, 3, 1, 3],
            "end": [1, 3, 1, 3],
            "state": ["A", "B", "A", "B"],
        }
    )

    x = tx.HSPELL_to_STS(
        df,
        id_col="id",
        start_col="start",
        end_col="end",
        state_col="state",
        fillblanks="*",
        nr="*",
    )

    assert x.nr == "*"
    assert "*" in x.alphabet
    assert x.nr_code == x.alphabet.index("*")
    out = seqtabstocc(x, include_missing=False)
    assert all("*" not in pattern for pattern in out["pattern"])


def test_seqformat_roundtrip_preserves_identifier_column():
    wide = pd.DataFrame(
        {
            "id": ["p1", "p2"],
            1: ["A", "B"],
            2: ["B", "B"],
            3: ["C", "B"],
        }
    )

    long = sf(wide, from_="wide", to="long", id_col="id")
    back = sf(long, from_="long", to="wide", id_col="id", time_col="time", state_col="state")

    assert list(back.columns) == ["id", 1, 2, 3]
    assert back["id"].tolist() == ["p1", "p2"]
    pd.testing.assert_frame_equal(back[[1, 2, 3]], wide[[1, 2, 3]])


def test_tied_event_timestamps_preserve_input_order_for_pairs():
    # The non-R convenience helper frequent_event_pairs keeps INPUT tie order (B,A).
    # (The faithful seqefsub, by contrast, comma-groups simultaneous events by event
    # code -> "(A,B)"; that behaviour is golden-pinned in test_r_goldens_seqefsub.py.)
    df = pd.DataFrame({"id": ["p1", "p1", "p1"], "timestamp": [1, 1, 2], "event": ["B", "A", "C"]})
    events = seqecreate(df, id_col="id", time_col="timestamp", event_col="event")

    pair_counts = {pair: count for pair, count, _ in frequent_event_pairs(events, min_support=1)}
    assert pair_counts[("B", "A")] == 1
    assert ("A", "B") not in pair_counts


def test_tied_event_timestamps_preserve_input_order_for_transition_counts_and_timegap_pairs():
    df = pd.DataFrame({"id": ["p1", "p1", "p1"], "timestamp": [1, 1, 2], "event": ["B", "A", "C"]})
    events = seqecreate(df, id_col="id", time_col="timestamp", event_col="event")

    labs, M = events_transition_matrix(events)
    idx = {lab: i for i, lab in enumerate(labs)}
    assert M[idx["B"], idx["A"]] == 1.0
    assert M[idx["A"], idx["C"]] == 1.0
    assert M[idx["A"], idx["B"]] == 0.0

    pair_counts = {
        pair: count
        for pair, count, _ in frequent_event_pairs_timegap(events, max_gap=1, min_support=1)
    }
    assert pair_counts[("B", "A")] == 1
    assert ("A", "B") not in pair_counts


def test_events_to_state_sequences_accepts_explicit_states_when_fill_appears_on_grid():
    df = pd.DataFrame({"id": ["p1", "p1"], "timestamp": [2, 4], "event": ["A", "B"]})
    events = seqecreate(df, id_col="id", time_col="timestamp", event_col="event")

    sts = events_to_state_sequences(
        events,
        states=["A", "B"],
        fill="_",
        align_missing=False,
        time_grid=np.arange(1, 5),
    )

    assert "_" in sts.alphabet
    assert _decode_states(sts).tolist() == [["_", "A", "A", "B"]]


@pytest.mark.parametrize("with_representatives", [False, True])
def test_cluster_chronogram_helpers_accept_state_labels(with_representatives):
    mat = np.array(
        [
            [0, 1, 2, 2],
            [0, 1, 1, 2],
            [2, 2, 1, 0],
            [2, 2, 2, 0],
        ],
        dtype=int,
    )
    x = tx.seqdef(mat, labels=["State 0", "State 1", "State 2"])
    groups = np.array([0, 0, 1, 1])

    if with_representatives:
        D = tx.seqdist(x, method="HAM")
        fig, axes = tx.cluster_chronograms_with_representatives(
            x, groups, D, include_missing=False, cols=2
        )
    else:
        fig, axes = tx.cluster_chronograms(x, groups, include_missing=False, cols=2)

    assert axes is not None
    plt.close(fig)


def test_seqlogp_raises_on_missing_without_include_missing():
    x = tx.seqdef(np.array([["A", "*", "*"], ["A", "A", "A"]], dtype=object), nr="*")
    with pytest.raises(ValueError):
        seqlogp(x, include_missing=False)


@pytest.mark.parametrize("variant", ["traminer", "tmean"])
def test_seqdist_rejects_unimplemented_dhd_variants(variant):
    x = tx.seqdef(np.array([[0, 1, 2, 2, 1], [0, 1, 1, 2, 2], [0, 2, 2, 2, 2]], dtype=int))

    with pytest.raises(ValueError, match="variant"):
        tx.seqdist(x, method="DHD", variant=variant)


def test_seqdist_rejects_unused_dur_eps_parameter():
    x = tx.seqdef(np.array([[0, 1, 2, 2, 1], [0, 1, 1, 2, 2], [0, 2, 2, 2, 2]], dtype=int))

    with pytest.raises(ValueError, match="eps"):
        tx.seqdist(x, method="DUR", metric="chi2", eps=1e-6)
