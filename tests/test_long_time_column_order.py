"""Time-column order of the long-to-wide pivot (audit S2-06).

``seqformat(from_="long", to="wide")`` -- and through it ``seqdef_from_long``
and ``TSE_to_STS`` in state mode -- ordered the time columns with Python's
``sorted()`` on the column labels. That sorts string labels as text, so
``t1..t12`` came out ``t1, t10, t11, t12, t2, ...``, and it did so even for an
ordered categorical time column, whose declared category order is the only way
pandas offers to state the order of text labels.

Source of truth: pandas' ordering contract for an ordered ``Categorical`` (its
categories, in declared order). The columns are now sorted with the time
index's own order, which is unchanged for numbers, dates and plain strings
(still text order; see ``seqdef_from_long``'s docstring) and follows the
categories for a categorical time column. There is no R counterpart:
``seqdef_from_long`` is a traminex extension.
"""

from __future__ import annotations

import pandas as pd

from traminex import TSE_to_STS, seqdef_from_long, seqformat

TIMES = [f"t{period}" for period in range(1, 13)]
STATES = [f"S{period:02d}" for period in range(1, 13)]


def _ordered_long(row_order: list[int] | None = None) -> pd.DataFrame:
    frame = pd.DataFrame({"id": [1] * 12, "time": TIMES, "state": STATES})
    if row_order is not None:
        frame = frame.iloc[row_order].reset_index(drop=True)
    return frame.assign(time=pd.Categorical(frame["time"], categories=TIMES, ordered=True))


def _states(x) -> list[str]:
    return [x._state_for_code(int(code)) for code in x.data[0]]


def test_seqdef_from_long_follows_an_ordered_categorical_time_column():
    assert _states(seqdef_from_long(_ordered_long())) == STATES


def test_row_order_does_not_matter_for_an_ordered_categorical():
    shuffled = _ordered_long([5, 0, 11, 3, 1, 2, 4, 6, 7, 8, 9, 10])
    assert _states(seqdef_from_long(shuffled)) == STATES


def test_tse_to_sts_state_mode_follows_an_ordered_categorical_time_column():
    assert _states(TSE_to_STS(_ordered_long())) == STATES


def test_long_to_wide_columns_follow_the_categories():
    wide = seqformat(
        _ordered_long(),
        from_="long",
        to="wide",
        id_col="id",
        time_col="time",
        state_col="state",
    )
    assert list(wide.columns) == ["id", *TIMES]


def test_plain_string_labels_keep_their_documented_text_order():
    plain = pd.DataFrame({"id": [1] * 12, "time": TIMES, "state": STATES})
    expected = [STATES[TIMES.index(label)] for label in sorted(TIMES)]
    assert _states(seqdef_from_long(plain)) == expected
