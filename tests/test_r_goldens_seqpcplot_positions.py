"""seqpcplot places events at R's time ranks and draws state sequences (TX-AUD-20261010-045, 11-027).

R's seqpcplot (order.align = "first") places each event at the rank of its
distinct time within its sequence, so simultaneous events share an x
(seqpcplot.R:209-211, 985-990), and draws a state sequence object through its
DSS (:131-140); seqplot(type = "pc") routes there. seqdss drops missing cells
before it collapses runs, and the y levels are the state labels (:142, :149-157;
batch-4 review I-3). traminex placed the events at 1..n and refused a state
sequence object. Fixture:
``scripts/r_goldens/make_goldens_seqpcplot_positions.R`` (seqpcplot's points).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402

HERE = Path(__file__).resolve().parent / "r_goldens"
GOLDEN = pd.read_csv(HERE / "seqpcplot_positions.csv")
TSE = pd.read_csv(HERE / "seqpcplot_positions_tse.csv")
STS = tx.seqdef(
    np.array([list(s) for s in ["AABBC", "CCAAB", "BBBBB", "ABABA"]], dtype=object),
    alphabet=["A", "B", "C"],
)
STS_MISSING = tx.seqdef(
    np.array([list(s) for s in ["A*ABC", "*CAAB", "B*B*B", "ABAB*"]], dtype=object),
    alphabet=["A", "B", "C"],
    missing="*",
    labels=["lab-A", "lab-B", "lab-C"],
)


def _r_paths(case: str) -> list[list[tuple[float, float]]]:
    rows = GOLDEN[GOLDEN["case"] == case].sort_values(["id", "position"])
    return [
        list(zip(group["x"].astype(float), group["y"].astype(float)))
        for _, group in rows.groupby("id")
    ]


def _paths(ax) -> list[list[tuple[float, float]]]:
    return [
        list(zip(map(float, line.get_xdata()), map(float, line.get_ydata())))
        for line in ax.get_lines()
    ]


def test_simultaneous_events_share_their_time_rank() -> None:
    eseq = tx.seqecreate(id=TSE["id"], timestamp=TSE["time"], event=TSE["event"])
    fig, ax = tx.seqpcplot(eseq)
    assert _paths(ax) == _r_paths("events")
    assert [label.get_text() for label in ax.get_yticklabels()] == ["a", "b", "c"]
    plt.close(fig)


def test_a_state_sequence_object_is_drawn_through_its_dss() -> None:
    fig, ax = tx.seqpcplot(STS)
    assert _paths(ax) == _r_paths("states")
    assert [label.get_text() for label in ax.get_yticklabels()] == ["A", "B", "C"]
    plt.close(fig)


def test_missing_cells_are_dropped_before_the_dss() -> None:
    fig, ax = tx.seqpcplot(STS_MISSING)
    assert _paths(ax) == _r_paths("states_missing")
    rows = GOLDEN[GOLDEN["case"] == "states_missing"].drop_duplicates("y").sort_values("y")
    assert [label.get_text() for label in ax.get_yticklabels()] == rows["ylevel"].tolist()
    plt.close(fig)


def test_seqplot_routes_type_pc() -> None:
    fig, ax = tx.seqplot(STS, type="pc")
    assert _paths(ax) == _r_paths("states")
    plt.close(fig)


def test_other_inputs_are_refused() -> None:
    with pytest.raises(TypeError, match="event-sequence .* or state-sequence"):
        tx.seqpcplot([1, 2, 3])
