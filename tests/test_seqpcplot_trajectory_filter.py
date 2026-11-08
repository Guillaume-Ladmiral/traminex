"""seqpcplot applies a seqpcfilter to distinct trajectories, as R does.

What was wrong (2026-09-22 audit S6-02): ``seqpcplot(filt=...)`` evaluated the
filter on a per-*sequence* vector in which every sequence carried its
trajectory's share. ``cumfreq`` then kept only as many *sequences* as were
needed to reach the level: with trajectory A>B at a share of 0.6, level 0.5
drew one of its three sequences and level 0.7 drew two A>B lines and no A>C.

Source of truth: TraMineR 2.2-14 ``seqpcplot.R``. Lines 350-364 build one
trajectory per distinct event pattern (events grouped by timestamp rank,
``order.align = "first"``), lines 515-517 compute one weighted share per
trajectory plus a final "no observation" share, and lines 633-639 evaluate the
filter on that vector. The expected colouring is the live R result in
``tests/r_goldens/seqpcplot_filter_trajectories.csv``
(``scripts/r_goldens/make_goldens_audit20260922_s6.R``).
"""

from __future__ import annotations

import csv
from itertools import groupby
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from traminex import seqecreate, seqpcfilter
from traminex.visualization.plots.events_viz import seqpcplot

FIXTURE = Path(__file__).with_name("r_goldens") / "seqpcplot_filter_trajectories.csv"

# The same inputs as make_goldens_audit20260922_s6.R (pc_cases).
CASES = {
    "shares_602020": (
        [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
        ["A", "B", "A", "B", "A", "B", "A", "C", "B", "C"],
        [1, 1, 1, 1, 1],
    ),
    "first_appearance_tie": (
        [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
        ["B", "C", "A", "B", "A", "B", "A", "B", "A", "C"],
        [1, 1, 1, 1, 1],
    ),
    "weighted": (
        [1, 1, 2, 2, 3, 3, 4, 4],
        [0, 1, 0, 1, 0, 1, 0, 1],
        ["A", "B", "A", "B", "B", "C", "A", "C"],
        [1, 1, 3, 0.5],
    ),
    "simultaneous": (
        [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 1],
        ["A", "B", "A", "B", "A", "B", "A", "B", "A", "B"],
        [1, 1, 1, 1, 1],
    ),
}


def _event_sequences(case: str):
    ids, timestamps, events, weights = CASES[case]
    sequences = seqecreate(pd.DataFrame({"id": ids, "timestamp": timestamps, "event": events}))
    sequences.weights = weights
    return sequences


def _r_filter_runs():
    with FIXTURE.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    runs = []
    for (case, method, level), members in groupby(
        rows, key=lambda row: (row["case"], row["method"], row["level"])
    ):
        members = sorted(members, key=lambda row: int(row["sequence"]))
        drawn = [row["drawn"] == "TRUE" for row in members]
        runs.append((case, method, None if level == "none" else float(level), drawn))
    return runs


R_RUNS = _r_filter_runs()


def _filter(method: str, level):
    return seqpcfilter(method) if level is None else seqpcfilter(method, level=level)


@pytest.mark.parametrize(
    ("case", "method", "level", "drawn"),
    R_RUNS,
    ids=[f"{case}-{method}-{level}" for case, method, level, _ in R_RUNS],
)
def test_filter_keeps_the_sequences_whose_trajectory_r_colours(case, method, level, drawn):
    from traminex.visualization.plots.events_viz import _filtered_sequences

    sequences = _event_sequences(case)

    keep = _filtered_sequences(list(sequences.sequences), _filter(method, level))

    assert keep.tolist() == drawn


@pytest.mark.parametrize(
    ("case", "method", "level", "drawn"),
    R_RUNS,
    ids=[f"{case}-{method}-{level}" for case, method, level, _ in R_RUNS],
)
def test_plot_draws_exactly_the_kept_sequences(case, method, level, drawn):
    sequences = _event_sequences(case)

    fig, ax = seqpcplot(sequences, filt=_filter(method, level))

    expected = [
        np.asarray(sequence.codes, dtype=float).tolist()
        for sequence, keep in zip(sequences.sequences, drawn)
        if keep
    ]
    assert [line.get_ydata().tolist() for line in ax.lines] == expected
    plt.close(fig)


def test_audit_reproduction_draws_every_sequence_of_the_kept_trajectory():
    # The audit's probe: A>B has a share of 0.6, so level 0.5 keeps all three
    # A>B sequences (it drew one) and level 0.7 adds A>C (it drew two A>B).
    sequences = _event_sequences("shares_602020")
    for level, count in ((0.5, 3), (0.7, 4)):
        fig, ax = seqpcplot(sequences, filt=seqpcfilter("cumfreq", level=level))
        assert len(ax.lines) == count
        plt.close(fig)
