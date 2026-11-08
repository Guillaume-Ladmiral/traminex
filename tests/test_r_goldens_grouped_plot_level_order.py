"""Grouped plots order text groups as R does (batch-2 review M-1, TX-REP-20261017-002).

R's seqedplot (seqedplot.R:53) and seqplot.tentrop (seqplot-tentrop.R:10) take
their panels and lines from ``factor(group)``; seqplot, behind seqHtplot,
from the internal ``group()`` (seqplot.R:163-184), which gives text the same
order with a space pasted before each level. Text sorts in R's collation
(ICU). traminex's shared panel helper sorted text by code point, and
seqplot_tentrop used pandas' categorical order, so ``_c``, ``A>%`` or ``é``
came out elsewhere than in R.

``group()`` also orders numeric groups as text (``1, 10, 2``) and titles them
with that space; traminex keeps ``factor()``'s numeric order and plain titles,
an R quirk corrected (ledger row seqHtplot, TX-REP-20261017-002).

Fixtures: ``scripts/r_goldens/make_goldens_factor_level_order.R`` (twelve
labels) and ``scripts/r_goldens/make_goldens_audit20260922_s6.R``
(``plot_group_levels.tsv``, numeric groups).
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import traminex as tx
from traminex.plotting import seqHtplot, seqplot_tentrop
from traminex.visualization.series import _factor_panels

GOLDENS = Path(__file__).resolve().parent / "r_goldens"
# The labels make_goldens_factor_level_order.R orders.
LABELS = ["b", "a", "B", "_c", "é", "e", "A>%", "A>*", "10", "9", "faible", "élevé"]


def _levels() -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    with (GOLDENS / "factor_level_order.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["case"] == "labels":
                out[row["kind"]].append(row["value"])
    return out


def _plot_group_levels() -> dict[tuple[str, str], list[str]]:
    with (GOLDENS / "plot_group_levels.tsv").open(newline="", encoding="utf-8") as handle:
        return {
            (row["case"], row["builder"]): row["levels"].split("|")
            for row in csv.DictReader(handle, delimiter="\t")
        }


R = _levels()
R_PLOT = _plot_group_levels()


def _seqs(n: int) -> tx.STSList:
    rows = [["A", "B", "A"] if k % 2 else ["B", "B", "A"] for k in range(n)]
    return tx.seqdef(rows, alphabet=["A", "B"])


def _legend(ax) -> list[str]:
    return [text.get_text() for text in ax.get_legend().get_texts()]


def test_r_seqplot_orders_text_groups_as_factor_does() -> None:
    # The premise of the seqHtplot test below, checked on R's own output.
    assert [level.removeprefix(" ") for level in R["seqplot_group_levels"]] == R["levels"]


def test_seqedplot_panels_order_text_groups_as_r_factor() -> None:
    rows, titles = _factor_panels(np.asarray(LABELS, dtype=object))
    assert titles == R["levels"]
    assert [LABELS[members[0]] for members in rows] == R["levels"]


def test_seqhtplot_orders_text_groups_as_r() -> None:
    fig, ax = seqHtplot(_seqs(len(LABELS)), group=LABELS)
    try:
        assert _legend(ax) == [level.removeprefix(" ") for level in R["seqplot_group_levels"]]
    finally:
        plt.close(fig)


def test_seqplot_tentrop_orders_text_groups_as_r_factor() -> None:
    fig, ax = seqplot_tentrop(_seqs(len(LABELS)), LABELS)
    try:
        assert _legend(ax) == R["levels"]
    finally:
        plt.close(fig)


def test_seqhtplot_orders_numeric_groups_by_value_not_as_r_seqplot_text() -> None:
    # make_goldens_audit20260922_s6.R section 4, case "numeric".
    group = [2, 10, 2, 10, 1, 1]
    fig, ax = seqHtplot(_seqs(len(group)), group=group)
    try:
        assert _legend(ax) == R_PLOT[("numeric", "factor")] == ["1", "2", "10"]
    finally:
        plt.close(fig)
    # R's seqplot orders and titles the same groups as text (the quirk corrected).
    assert R_PLOT[("numeric", "TraMineR:::group")] == [" 1", " 10", " 2"]


@pytest.mark.parametrize("plot", ["seqHtplot", "seqplot_tentrop"])
def test_mixed_text_and_number_groups_are_refused(plot: str) -> None:
    seqs = _seqs(3)
    with pytest.raises(TypeError, match="mix text and numbers"):
        if plot == "seqHtplot":
            seqHtplot(seqs, group=["a", 1, "b"])
        else:
            seqplot_tentrop(seqs, ["a", 1, "b"])
