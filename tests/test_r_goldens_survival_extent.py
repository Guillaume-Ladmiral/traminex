"""Survival curves run to each stratum's last follow-up time (TX-AUD-20261010-042, 11-026).

R draws seqsplot's and seqedplot's survival curves with plot.survfit: each
starts at (0, 1), steps at its event times and runs on to its last time,
censored or not, so a state whose spells are all censored is a flat line at 1.
traminex drew only the event rows: curves stopped at the last event and a state
without events had no curve. Fixture: ``scripts/r_goldens/make_goldens_survival_extent.R``
(R's drawn coordinates, traced from lines()).
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
from traminex.events import seqecreate  # noqa: E402

HERE = Path(__file__).resolve().parent / "r_goldens"
GOLDEN = pd.read_csv(HERE / "survival_extent.csv", dtype=str, keep_default_na=False)
MVAD = pd.read_csv(HERE / "survival_extent_mvad.csv", dtype=str, keep_default_na=False)
STATES = MVAD.drop(columns=["male"]).to_numpy(dtype=object)
# The generator's alphabet: R's sort(unique(...)) of the mvad states.
MV = tx.seqdef(STATES, alphabet=["employment", "FE", "HE", "joblessness", "school", "training"])
TOY = tx.seqdef(
    np.array([list(s) for s in ["AABBC", "ABBBC", "BBAAC", "AABBB", "BBBAA"]], dtype=object),
    alphabet=["A", "B", "C"],
)


def _eseq():
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s2", "s3", "s3", "s4", "s4", "s5"],
            "timestamp": [1, 3, 2, 1, 4, 1, 3, 2],
            "event": ["A", "B", "A", "B", "A", "A", "A", "B"],
        }
    )
    eseq = seqecreate(tse)
    eseq.lengths = [5, 5, 6, 4, 5]
    return eseq


def _jumps(xs, ys) -> list[tuple[float, float]]:
    """(x, value from x on) for a drawn staircase, consecutive repeats merged."""
    points: list[tuple[float, float]] = []
    for x, y in zip(xs, ys):
        if points and np.isclose(points[-1][0], x, rtol=0, atol=1e-12):
            points[-1] = (points[-1][0], float(y))
        else:
            points.append((float(x), float(y)))
    merged = [points[0]]
    for x, y in points[1:]:
        if np.isclose(merged[-1][1], y, rtol=0, atol=1e-12) and x != merged[-1][0]:
            merged.append((x, y))  # a flat run to a later time is kept: it is the extent
        elif not np.isclose(merged[-1][1], y, rtol=0, atol=1e-12):
            merged.append((x, y))
    return merged


def _r_curves(case: str) -> dict[int, list[list[tuple[float, float]]]]:
    panels: dict[int, list] = {}
    for _, row in GOLDEN[GOLDEN["case"] == case].iterrows():
        xs = [float(v) for v in row["x"].split("|")]
        ys = [float(v) for v in row["y"].split("|")]
        panels.setdefault(int(row["panel"]), []).append(_jumps(xs, ys))
    return panels


def _ours(ax) -> list[list[tuple[float, float]]]:
    return [_jumps(line.get_xdata(), line.get_ydata()) for line in ax.get_lines()]


def _same(ours, r) -> None:
    assert len(ours) == len(r)
    for got, want in zip(ours, r):
        assert len(got) == len(want)
        np.testing.assert_allclose(np.array(got), np.array(want), rtol=0, atol=1e-12)


@pytest.mark.parametrize("case, seqs", [("toy", TOY), ("mvad150", MV)])
def test_seqsurvplot_curves_run_to_the_last_follow_up(case: str, seqs) -> None:
    fig, ax = tx.seqsurvplot(seqs)
    _same(_ours(ax), _r_curves(case)[1])
    plt.close(fig)


def test_a_state_with_censored_spells_only_is_flat_at_one() -> None:
    fig, ax = tx.seqsurvplot(TOY)
    flat = _ours(ax)[2]
    assert flat == [(0.0, 1.0), (1.0, 1.0)]
    plt.close(fig)


def test_seqsplot_group_panels_run_to_the_last_follow_up() -> None:
    fig, axes = tx.seqsplot(MV, group=MVAD["male"].tolist())
    r = _r_curves("mvad150_male")
    for panel, ax in enumerate(axes[: len(r)], start=1):
        _same(_ours(ax), r[panel])
    plt.close(fig)


def test_seqedplot_survival_runs_to_the_last_follow_up() -> None:
    fig, axes = tx.seqedplot(_eseq(), type="survival")
    ax = axes[0] if isinstance(axes, (list, tuple, np.ndarray)) else axes
    _same(_ours(ax), _r_curves("seqedplot")[1])
    plt.close(fig)
