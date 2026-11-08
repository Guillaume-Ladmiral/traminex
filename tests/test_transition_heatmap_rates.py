"""transition_heatmap draws seqtrate's rates (TX-AUD-20261011-028).

transition_heatmap is a traminex helper whose docstring said its rates were
seqtrate's; it counted pairs over every code, the missing and void markers
included, always weighted (row B 0.2449/0.449/0.1224/0/0.1837 where seqtrate,
equal to R, gives 0.3/0.55/0.15). It now draws seqtrate's rates (counts with
``normalize=False``), and takes ``weighted=``. The expected values are R's
seqtrate on the same input (``scripts/r_goldens/make_goldens_plot_items.R``, rows
``heatmap``; TX-REP-20261019-003), not traminex's own.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402

SEQS = tx.seqdef(
    np.array(
        [list("AABBC"), ["A", "B", None, "B", "C"], list("CCBAA"), ["B", "B", "A", None, None]],
        dtype=object,
    ),
    weights=[1.0, 2.0, 1.0, 0.5],
)


with open(
    Path(__file__).resolve().parent / "r_goldens" / "plot_items.csv", encoding="utf-8", newline=""
) as handle:
    GOLDEN = {
        row["field"]: row["value"] for row in csv.DictReader(handle) if row["case"] == "heatmap"
    }
STATES = GOLDEN["states"].split(" ;; ")


def _r_matrix(field: str) -> np.ndarray:
    return np.array([float(v) for v in GOLDEN[field].split()]).reshape(len(STATES), len(STATES))


def _drawn(ax) -> np.ndarray:
    return np.asarray(ax.get_images()[0].get_array(), dtype=float)


def test_the_input_is_the_golden_s() -> None:
    rows = [row.split() for row in GOLDEN["cells"].split(" ;; ")]
    decoded = [[SEQS.alphabet[int(code)] for code in row] for row in SEQS.data]
    markers = {SEQS.alphabet[SEQS.nr_code], SEQS.alphabet[SEQS.void_code]}
    assert len(decoded) == len(rows)
    for ours, theirs in zip(decoded, rows):
        assert len(ours) == len(theirs)
        assert all(o in markers if t == "NA" else o == t for o, t in zip(ours, theirs))
    assert [float(w) for w in GOLDEN["weights"].split()] == list(SEQS.weights)


@pytest.mark.parametrize("weighted", [True, False])
def test_rates_are_r_seqtrate_rates(weighted: bool) -> None:
    fig, ax = tx.transition_heatmap(SEQS, weighted=weighted)
    expected = _r_matrix("rates_weighted" if weighted else "rates_unweighted")
    np.testing.assert_allclose(_drawn(ax), expected, rtol=0, atol=1e-12)
    assert [t.get_text() for t in ax.get_yticklabels()] == STATES
    assert [t.get_text() for t in ax.get_xticklabels()] == STATES
    plt.close(fig)


def test_counts_are_r_seqtrate_counts() -> None:
    fig, ax = tx.transition_heatmap(SEQS, normalize=False)
    np.testing.assert_allclose(_drawn(ax), _r_matrix("counts_weighted"), rtol=0, atol=1e-12)
    plt.close(fig)
