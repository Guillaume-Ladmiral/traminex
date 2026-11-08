"""seqplot_rf prints R's squared-distance statistic (TX-AUD-20261010-041).

TraMineRextras' seqplot.rf computes its own quality statistic from squared
distances to the group medoids and to the global medoid (seqplot-rf.R:22-24,
139-146), prints it with format() and writes it into the subtitle. traminex
printed seqrf's statistic under the same message. Fixture:
``scripts/r_goldens/make_goldens_seqplot_rf_stat.R`` (distinct MDS scores, so R's
random tie-breaking does not enter).
"""

from __future__ import annotations

import contextlib
import io
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402

# k = n leaves no residual degrees of freedom; seqrf repeats R's dissrf warning.
pytestmark = pytest.mark.filterwarnings(r"ignore:\s*\[!!\] Negative  or zero df:UserWarning")

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqplot_rf_stat.csv",
    dtype=str,
    keep_default_na=False,
)
STRINGS = ["AAABBC", "AABBCC", "ABBBCC", "AAAABB", "BBBCCC", "ABCABC", "AACCBB", "CCBBAA",
           "ABBBBB", "AAACCC", "BCCCCC", "AAAAAC", "BBAACC", "CABBBA"]  # fmt: skip
SEQS = tx.seqdef(np.array([list(s) for s in STRINGS], dtype=object), alphabet=["A", "B", "C"])
DISS = np.asarray(tx.seqdist(SEQS, method="LCS"), dtype=float)


def _sortv(row: dict) -> np.ndarray:
    # R's cmdscale(diss, k = 1); the sign of its axis is TX-REP-20261017-006.
    return np.array([float(value) for value in row["sortv"].split("|")])


def _run(row: dict):
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        options = {} if row["k"] == "default" else {"k": int(row["k"])}
        fig, axes = tx.seqplot_rf(SEQS, DISS, sortv=_sortv(row), **options)
    return fig, stderr.getvalue()


@pytest.mark.parametrize("row", GOLDEN.to_dict("records"), ids=lambda r: f"k{r['k']}")
def test_printed_statistic_is_r_seqplot_rf(row: dict) -> None:
    fig, text = _run(row)
    assert re.search(r"Pseudo/median-based-R2: (\S+)", text).group(1) == row["R2"]
    assert re.search(r"Pseudo/median-based-F statistic: (\S+)", text).group(1) == row["Fstat"]
    plt.close(fig)


@pytest.mark.parametrize("row", GOLDEN.to_dict("records"), ids=lambda r: f"k{r['k']}")
def test_subtitle_reports_r_statistic(row: dict) -> None:
    fig, _ = _run(row)
    r2, fstat = (
        row[name] if row[name] == "NaN" else f"{float(row[name]):.2f}" for name in ("R2", "Fstat")
    )
    quality = f"Representation quality: R2={r2} and F={fstat}"  # R's sprintf("%0.2f")
    assert fig._suptitle is not None and quality in fig._suptitle.get_text()
    plt.close(fig)


@pytest.mark.parametrize(
    "row",
    [r for r in GOLDEN.to_dict("records") if r["k"] not in {"1", "default"}],
    ids=lambda r: f"k{r['k']}",
)
def test_groups_and_medoids_are_r_seqplot_rf(row: dict) -> None:
    fitted = tx.seqrf(
        SEQS, DISS, k=int(row["k"]), weighted=False, grp_meth="first", sortv=_sortv(row)
    )
    medoid_index = [int(i) + 1 for i in fitted["rf"]["kmedoid_index"]]
    assert medoid_index == [int(i) for i in row["kmedoid_index"].split("|")]
    assert int(tx.disscenter(DISS, medoids_index="first")) + 1 == int(row["global_medoid"])
