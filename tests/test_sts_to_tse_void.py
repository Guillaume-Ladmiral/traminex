"""STS_to_TSE writes no record for a void cell (TX-AUD-20261019-082).

``STS_to_TSE`` is traminex's state-at-time table (its ledger row: a traminex contract). A void cell is the end of a
shorter sequence, not an observed state, yet it was written as a record whose state is the void marker, and
``TSE_to_STS`` read that back as a real state ``%`` (alphabet ``['%', 'A', 'B']``). Void cells now get no record; a
missing (nr) cell keeps its record, as R's event conversion keeps a transition into missing. Times are
``x.start + position``, which ``TSE_to_STS`` reads back (R's event conversion writes 0-based event times).
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx

RAGGED = np.array(
    [["A", "A", "B", None], ["A", "B", None, None], ["B", np.nan, "B", "A"]], dtype=object
)
ALL_CELLS = [(1, 1, "A"), (1, 2, "A"), (1, 3, "B"), (2, 1, "A"), (2, 2, "B")]
ALL_CELLS += [(3, 1, "B"), (3, 2, "*"), (3, 3, "B"), (3, 4, "A")]
CHANGES = [
    (1, 1, "A"),
    (1, 3, "B"),
    (2, 1, "A"),
    (2, 2, "B"),
    (3, 1, "B"),
    (3, 2, "*"),
    (3, 3, "B"),
    (3, 4, "A"),
]


@pytest.mark.parametrize("changes_only", [False, True])
def test_void_cells_get_no_record_and_missing_cells_keep_theirs(changes_only: bool) -> None:
    x = tx.seqdef(RAGGED)
    table = tx.STS_to_TSE(x, state_changes_only=changes_only)
    assert list(table.itertuples(index=False, name=None)) == (
        CHANGES if changes_only else ALL_CELLS
    )


def test_the_state_table_reads_back_its_void_cells() -> None:
    # One row fills the table: a column of void cells only has no record to read back.
    x = tx.seqdef(
        np.array(
            [["A", "A", "B", "B"], ["A", "B", None, None], ["B", None, None, None]], dtype=object
        )
    )
    back = tx.TSE_to_STS(tx.STS_to_TSE(x, state_changes_only=False))
    assert list(back.alphabet) == list(x.alphabet)
    assert [[back.alphabet[c] for c in row] for row in back.data] == [
        [x.alphabet[c] for c in row] for row in x.data
    ]


def test_times_count_from_the_sequences_start() -> None:
    table = tx.STS_to_TSE(tx.seqdef(RAGGED[:2], start=5), state_changes_only=False)
    assert table.time.tolist() == [5, 6, 7, 5, 6]


def test_an_internal_void_cell_ends_the_run_in_change_mode() -> None:
    # An internal void (only a recode makes one) gets no record either, and the state after it starts a new run.
    x = tx.seqdef(np.array([["A", "A", "A", "B"], ["B", None, None, None]], dtype=object))
    x.data[0, 1] = x.void_code
    table = tx.STS_to_TSE(x)
    assert list(table.itertuples(index=False, name=None)) == [
        (1, 1, "A"),
        (1, 3, "A"),
        (1, 4, "B"),
        (2, 1, "B"),
    ]


@pytest.mark.parametrize("time_col", ["time", "t"])
def test_a_table_without_records_keeps_its_columns(time_col: str) -> None:
    # A sequence of void cells only has no record; a table without records still has its columns, so TSE_to_STS
    # does not fail on missing columns (batch-11 review M-8).
    table = tx.STS_to_TSE(tx.seqdef(np.array([[None, None]], dtype=object)), time_col=time_col)
    assert list(table.columns) == ["id", time_col, "state"]
    assert len(table) == 0


def test_the_docstring_states_when_the_table_reads_back() -> None:
    # The void cells read back only from one record per cell; a change table does not read back (review M-7).
    text = " ".join(tx.STS_to_TSE.__doc__.split())
    assert (
        "With ``state_changes_only=False`` (one record per cell), :func:`TSE_to_STS` reads" in text
    )
    assert "a change table does not read back" in text
