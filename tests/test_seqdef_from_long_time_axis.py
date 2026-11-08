"""seqdef_from_long keeps its time axis and can reindex onto a grid (TX-AUD-20261010-018).

seqdef_from_long (a traminex constructor; R's seqformat has no long format)
pivots person-period rows to one column per observed time and dropped those
column names, so the object lost its time axis (xtlab None) and a period
no person reports vanished without trace. The time values now become xtlab,
and ``time_grid=`` reindexes onto a given grid (a period no person reports is
not invented otherwise: wave-coded data such as 1, 3, 5 has no period 2).
Date and date-time times are labelled as R names them, the date alone at
midnight; fixture ``scripts/r_goldens/make_goldens_datetime_labels.R``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import traminex as tx

LONG = pd.DataFrame(
    {"id": [1, 1, 1, 2, 2], "t": [2001, 2002, 2005, 2001, 2003], "st": ["X", "X", "X", "Y", "Y"]}
)


LABELS = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "datetime_labels.csv",
    dtype=str,
    keep_default_na=False,
)


def _build(**kwargs) -> tx.STSList:
    return tx.seqdef_from_long(LONG, id_col="id", time_col="t", state_col="st", **kwargs)


def test_time_values_become_xtlab() -> None:
    seqs = _build()
    assert seqs.xtlab == ["2001", "2002", "2003", "2005"]
    assert seqs.to_strings() == ["X-X-*-X", "Y-*-Y-%"]


def test_time_grid_reindexes_onto_the_grid() -> None:
    seqs = _build(time_grid=range(2001, 2006))
    assert seqs.seqlen == 5
    assert seqs.xtlab == ["2001", "2002", "2003", "2004", "2005"]
    assert seqs.to_strings() == ["X-X-*-*-X", "Y-*-Y-%-%"]


def test_time_grid_refuses_a_time_outside_it() -> None:
    with pytest.raises(ValueError, match="2005"):
        _build(time_grid=range(2001, 2004))


def test_an_explicit_xtlab_wins() -> None:
    assert _build(xtlab=["a", "b", "c", "d"]).xtlab == ["a", "b", "c", "d"]


def test_time_grid_refuses_duplicate_entries() -> None:
    # A repeated grid entry repeated the period (batch-3 review M3).
    with pytest.raises(ValueError, match="duplicate.*2001"):
        _build(time_grid=[2001, 2001, 2002, 2003, 2004, 2005])


@pytest.mark.parametrize("case", sorted(set(LABELS["case"])))
def test_date_and_datetime_times_are_labelled_as_r_names_them(case: str) -> None:
    rows = LABELS[LABELS["case"] == case].sort_values("position", key=lambda c: c.astype(int))
    times = pd.to_datetime(rows["time"].tolist())
    if rows["tz"].iloc[0]:
        times = times.tz_localize(rows["tz"].iloc[0])
    long = pd.DataFrame({"id": [1, 1, 1], "t": times, "st": ["A", "B", "A"]})
    seqs = tx.seqdef_from_long(long, id_col="id", time_col="t", state_col="st")
    assert seqs.xtlab == rows["label"].tolist()
