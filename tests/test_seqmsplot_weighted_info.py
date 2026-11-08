"""seqmsplot's info text divides by the weighted number of sequences.

What was wrong (2026-09-22 audit S6-05): the "Modal state sequence (k
occurrences, freq=...%)" annotation divided the occurrence count by the number
of rows, so weighted data showed ``freq=25.0%`` where TraMineR shows 17.8%. The
unweighted text also printed ``25.0`` where R prints ``25``.

Source of truth: TraMineR 2.2-14 ``plot.stslist.modst.R`` computes
``ctfreq <- round((nbrep / n) * 100, 1)`` with ``n <- attr(x, "nbseq")``, which
``seqmodst.R:65`` sets to the sum of the weights, and pastes it as text. The
expected strings are live R output in ``tests/r_goldens/seqmsplot_info_text.tsv``
(``scripts/r_goldens/make_goldens_audit20260922_s6.R``).
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from traminex import seqdef
from traminex.plotting import seqmsplot

FIXTURE = Path(__file__).with_name("r_goldens") / "seqmsplot_info_text.tsv"

# The same data as make_goldens_audit20260922_s6.R (section 3).
ROWS = [
    ["A", "A", "B", "B", "C", "C"],
    ["A", "B", "B", "C", "C", "C"],
    ["B", "B", None, "C", "D", "D"],
    ["A", "A", "A", "A", "B", "B"],
    ["C", "C", "C", "D", "D", None],
    ["A", None, "B", "B", "C", "C"],
    ["D", "D", "D", "D", "D", "D"],
    ["A", "A", "B", "B", "C", "C"],
]
WEIGHTS = {"weighted": np.array([1.0, 2.0, 0.5, 3.0, 1.5, 1.0, 0.25, 2.0]), "unweighted": None}


def _r_texts():
    with FIXTURE.open(newline="") as handle:
        return {row["case"]: row["text"] for row in csv.DictReader(handle, delimiter="\t")}


R_TEXTS = _r_texts()


@pytest.mark.parametrize("case", ["weighted", "unweighted"])
def test_info_text_matches_r(case):
    x = seqdef(ROWS, alphabet=["A", "B", "C", "D"], weights=WEIGHTS[case])

    fig, ax = seqmsplot(x)

    assert [text.get_text() for text in ax.texts] == [R_TEXTS[case]]
    plt.close(fig)
