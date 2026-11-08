"""seqmtplot draws R's bars, without error bars by default (TX-AUD-20261010-043, -058).

R's seqmtplot is seqplot(type="mt"), which calls seqmeant with serr = FALSE by
default (seqplot.R:345-346; seqmeant.R:5), so no error bars are drawn unless
asked, and draws one bar per seqmeant row: the alphabet, plus the missing state
with with.missing = TRUE. traminex drew 95% error bars by default and one bar
per marker-padded seqmeant row (A, B, *, % for ragged data). Fixture:
``scripts/r_goldens/make_goldens_seqmeant_rows.R``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from matplotlib.container import ErrorbarContainer  # noqa: E402

import traminex as tx  # noqa: E402

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqmeant_rows.csv",
    dtype=str,
    keep_default_na=False,
)
DATA = {
    "ragged": tx.seqdef(np.array([["A", None, "B", "B"], ["A", "B", None, None]], dtype=object)),
    "complete": tx.seqdef(np.array([["A", "B", "C"], ["B", "B", "A"]], dtype=object)),
}


def _r_rows(case: str, with_missing: bool) -> list[str]:
    rows = GOLDEN[(GOLDEN["case"] == case) & (GOLDEN["with_missing"] == str(with_missing).upper())]
    return rows.sort_values("position", key=lambda c: c.astype(int))["row"].tolist()


@pytest.mark.parametrize("case", sorted(DATA))
@pytest.mark.parametrize("with_missing", [False, True])
def test_bars_are_r_seqmeant_rows(case: str, with_missing: bool) -> None:
    fig, ax = tx.seqmtplot(DATA[case], include_missing=with_missing)
    assert [tick.get_text() for tick in ax.get_xticklabels()] == _r_rows(case, with_missing)
    assert len(ax.patches) == len(_r_rows(case, with_missing))
    plt.close(fig)


def test_no_error_bars_by_default_as_r() -> None:
    assert set(GOLDEN["serr_default"]) == {"FALSE"}
    fig, ax = tx.seqmtplot(DATA["complete"])
    assert not any(isinstance(container, ErrorbarContainer) for container in ax.containers)
    plt.close(fig)


def test_error_bars_when_asked() -> None:
    fig, ax = tx.seqmtplot(DATA["complete"], errbar=True)
    assert sum(isinstance(container, ErrorbarContainer) for container in ax.containers) == 1
    plt.close(fig)
