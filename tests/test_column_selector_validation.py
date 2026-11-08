"""HSPELL_to_STS checks its column selectors as R does (TX-AUD-20261010-036).

An out-of-range numeric selector raised a raw IndexError and a selector that
picked a text column as a begin or end time raised pandas' "Unable to parse
string" without naming the column. R stops with "undefined columns selected"
for a column that does not exist. HSPELL_to_STS's numeric selectors stay
zero-based (docs/API_PARITY.md lists which functions use which base).
"""

from __future__ import annotations

import pandas as pd
import pytest

import traminex as tx

FRAME = pd.DataFrame({"id": [1, 2], "b1": [1, 2], "e1": [3, 4], "s1": ["A", "B"], "extra": [0, 0]})


def test_an_out_of_range_selector_is_undefined() -> None:
    with pytest.raises(ValueError, match="undefined columns selected"):
        tx.HSPELL_to_STS(FRAME, id_col="id", begin=[9], end=[2], status=[3])


def test_a_negative_selector_beyond_the_width_is_undefined() -> None:
    with pytest.raises(ValueError, match="undefined columns selected"):
        tx.HSPELL_to_STS(FRAME, id_col="id", begin=[1], end=[-9], status=[3])


def test_a_text_column_used_as_a_time_is_named() -> None:
    # R-style one-based selectors: 2 and 3 are e1 and s1 here.
    with pytest.raises(ValueError, match="'s1'"):
        tx.HSPELL_to_STS(FRAME, id_col="id", begin=[2], end=[3], status=[4])


def test_zero_based_and_negative_selectors_in_range_still_work() -> None:
    by_name = tx.HSPELL_to_STS(FRAME, id_col="id", begin=["b1"], end=["e1"], status=["s1"])
    by_position = tx.HSPELL_to_STS(FRAME, id_col="id", begin=[1], end=[-3], status=[3])
    assert by_position.to_strings() == by_name.to_strings()


@pytest.mark.parametrize(
    "option",
    [
        {"id_col": "ident"},
        {"tmin_col": "low"},
        {"tmax_col": "high"},
        {"endObs": "eobs"},
        {"pvar": "orign"},
    ],
    ids=lambda option: next(iter(option)),
)
def test_a_misspelled_horizontal_column_name_is_undefined(option: dict) -> None:
    # R stops on seqdata[, pvar], seqdata[, endObs] and seqdata[, id]; a
    # misspelled tmin_col or tmax_col used to be ignored, the others raised a raw
    # KeyError (batch-3 review M6).
    arguments = {"id_col": "id", "begin": ["b1"], "end": ["e1"], "status": ["s1"], **option}
    with pytest.raises(ValueError, match="undefined columns selected"):
        tx.HSPELL_to_STS(FRAME, **arguments)


VERTICAL = pd.DataFrame(
    {"id": ["p1", "p1"], "start": [1, 3], "end": [2, 4], "state": ["A", "B"], "eo": [4, 4]}
)


@pytest.mark.parametrize(
    "option",
    [
        {"id_col": "ident"},
        {"start_col": "begin"},
        {"end_col": "stop"},
        {"state_col": "status"},
        {"tmin_col": "low"},
        {"tmax_col": "high"},
        {"endObs": "eobs"},
        {"pvar": "orign"},
    ],
    ids=lambda option: next(iter(option)),
)
def test_a_misspelled_vertical_column_name_is_undefined(option: dict) -> None:
    arguments = {
        "id_col": "id",
        "start_col": "start",
        "end_col": "end",
        "state_col": "state",
        **option,
    }
    with pytest.raises(ValueError, match="undefined columns selected"):
        tx.HSPELL_to_STS(VERTICAL, **arguments)
