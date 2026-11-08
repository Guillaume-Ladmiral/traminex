import numpy as np
import pandas as pd

import traminex as tx
from traminex import seqdef_from_long
from traminex import seqformat as sf


def test_seqformat_roundtrip():
    wide = pd.DataFrame(
        {
            "id": ["p1", "p2"],
            1: ["A", "B"],
            2: ["B", "B"],
            3: ["C", "B"],
        }
    )
    long = sf(wide, from_="wide", to="long", id_col="id")
    assert set(long.columns) == {"id", "time", "state"}
    back = sf(long, from_="long", to="wide", id_col="id", time_col="time", state_col="state")
    assert back.shape == (2, 4)
    assert list(back.columns) == ["id", 1, 2, 3]
    assert back["id"].tolist() == ["p1", "p2"]


def test_seqdef_from_long():
    long = pd.DataFrame(
        {
            "id": ["p1", "p1", "p2", "p2"],
            "time": [1, 2, 1, 2],
            "state": ["A", "B", "B", "B"],
        }
    )
    x = seqdef_from_long(long)
    assert x.nseq == 2 and x.seqlen == 2


def test_seqdef_from_long_keeps_ids_aligned_in_first_appearance_order():
    # R's route keeps the persons in first-appearance order (TX-AUD-20261019-039).
    long = pd.DataFrame(
        {
            "id": ["case-z", "case-a", "case-z", "case-a"],
            "time": [1, 1, 2, 2],
            "state": ["A", "B", "C", "A"],
        }
    )

    x = seqdef_from_long(long)
    decoded = np.asarray(x.alphabet, dtype=object)[x.data]

    assert x.ids == ["case-z", "case-a"]
    assert dict(zip(x.ids, decoded.tolist(), strict=True)) == {
        "case-a": ["B", "A"],
        "case-z": ["A", "C"],
    }


def test_seqformat_wide_to_long_auto_id_one_based():
    wide = pd.DataFrame({"t1": ["A", "B"], "t2": ["B", "C"]})
    long = sf(wide, from_="wide", to="long")

    assert long["id"].min() == 1
    assert set(long["id"].unique()) == {1, 2}


def test_seqstart_preserves_labels_when_fill_state_added():
    x = tx.seqdef(
        np.array([["A", "B"], ["A", "B"]], dtype=object),
        labels=["State A", "State B"],
        nr="*",
    )

    y = tx.seqstart(x, data_start=[1, 1], new_start=[1, 2])

    assert y.labels is not None
    assert len(y.labels) == len(y.alphabet)
    state_to_label = dict(zip(y.alphabet, y.labels))
    assert state_to_label["A"] == "State A"
    assert state_to_label["B"] == "State B"
    assert "*" in y.alphabet
