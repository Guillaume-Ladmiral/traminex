"""STS_to_TSE's id column carries the object's ids (TX-AUD-20261010-033).

traminex's STS_to_TSE is its own state-at-time table (no tevent; ledger row
STS_to_TSE). It numbered the sequences 1..n even when the object had ids. R's
STS_to_TSE does that only when no id is passed, and seqformat(to="TSE") warns
then that it replaces the original ids; traminex has no id argument, so it uses
the object's ids, and 1..n (R's default, from the golden) when there are none.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import traminex as tx

DATA = np.array([list("AAB"), list("ABB")], dtype=object)
GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "sts_sps_dss_void.csv", keep_default_na=False
)


def test_sts_to_tse_uses_the_object_ids() -> None:
    table = tx.STS_to_TSE(tx.seqdef(DATA, id=["p1", "p2"]))
    assert table["id"].tolist() == ["p1", "p1", "p2", "p2"]


def test_sts_to_tse_without_ids_numbers_the_sequences_as_r() -> None:
    table = tx.STS_to_TSE(tx.seqdef(DATA))
    r_ids = GOLDEN[GOLDEN["form"] == "TSE_id_default"].sort_values("position")["token"]
    assert (
        [str(value) for value in table["id"]]
        == [str(value) for value in r_ids]
        == ["1", "1", "2", "2"]
    )


def test_sts_to_tse_keeps_the_input_order_of_the_sequences() -> None:
    # R's STS_to_TSE loops over the sequences in order; sorting by id would put
    # "x" before "z" now that the ids are the object's.
    table = tx.STS_to_TSE(tx.seqdef(DATA, id=["z", "x"]))
    assert table["id"].tolist() == ["z", "z", "x", "x"]
    assert table["time"].tolist() == [1, 3, 1, 2]
