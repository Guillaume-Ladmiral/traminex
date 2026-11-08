"""Vertical spells give the horizontal route's result (TX-AUD-20261011-019).

The vertical layout (one row per spell) is a traminex extension of
TraMineRextras::HSPELL_to_STS, whose horizontal layout (one row per person) the
horizontal route reproduces. The vertical route had its own rules: a "_" filler
that became a real state, its own endObs rule, no time labels, a raw error for
a missing end, and ids sorted as text. Each test describes the same spells
once, builds both layouts, and requires the same object.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import traminex as tx


def _layouts(persons, extras=None):
    """Horizontal and vertical frames for ``persons``: (id, [(begin, end, state), ...])."""
    extras = extras or {}
    width = max(len(spells) for _, spells in persons)
    horizontal = {"id": [pid for pid, _ in persons]}
    for k in range(width):
        horizontal[f"b{k + 1}"] = [s[k][0] if k < len(s) else np.nan for _, s in persons]
        horizontal[f"e{k + 1}"] = [s[k][1] if k < len(s) else np.nan for _, s in persons]
        horizontal[f"s{k + 1}"] = [s[k][2] if k < len(s) else None for _, s in persons]
    for name, values in extras.items():
        horizontal[name] = values
    vertical = {"id": [], "start": [], "end": [], "state": []}
    for name in extras:
        vertical[name] = []
    for index, (pid, spells) in enumerate(persons):
        for begin, end, state in spells:
            vertical["id"].append(pid)
            vertical["start"].append(begin)
            vertical["end"].append(end)
            vertical["state"].append(state)
            for name, values in extras.items():
                vertical[name].append(values[index])
    spell_cols = {
        "begin": [f"b{k + 1}" for k in range(width)],
        "end": [f"e{k + 1}" for k in range(width)],
        "status": [f"s{k + 1}" for k in range(width)],
    }
    return pd.DataFrame(horizontal), spell_cols, pd.DataFrame(vertical)


def _view(x: tx.STSList):
    markers = {x.nr_code: "*", x.void_code: "%"}
    cells = [[markers.get(int(c), x.alphabet[int(c)]) for c in row] for row in x.data]
    return x.ids, None if x.xtlab is None else list(x.xtlab), x.alphabet, cells


def _both(persons, extras=None, **options):
    horizontal, spell_cols, vertical = _layouts(persons, extras)
    h = tx.HSPELL_to_STS(horizontal, id_col="id", **spell_cols, **options)
    v = tx.HSPELL_to_STS(
        vertical, id_col="id", start_col="start", end_col="end", state_col="state", **options
    )
    return _view(h), _view(v)


PERSON = [("p1", [(2000, 2002, "X"), (2003, 2004, "Y")])]


def test_end_of_observation() -> None:
    horizontal, vertical = _both(PERSON, {"eo": [2003]}, endObs="eo")
    assert vertical == horizontal


def test_no_end_of_observation() -> None:
    horizontal, vertical = _both(PERSON)
    assert vertical == horizontal


def test_uncovered_periods_are_missing_not_a_state() -> None:
    persons = [("p1", [(1, 2, "A"), (5, 6, "A")]), ("p2", [(2, 3, "B")])]
    horizontal, vertical = _both(persons)
    assert vertical == horizontal
    assert "_" not in vertical[2]


def test_a_missing_end_runs_the_spell_to_tmax() -> None:
    persons = [("p1", [(2000, np.nan, "X")]), ("p2", [(2000, 2003, "Y")])]
    horizontal, vertical = _both(persons)
    assert vertical == horizontal


def test_numeric_ids_keep_their_input_order() -> None:
    persons = [(2, [(1, 2, "X")]), (10, [(1, 2, "Y")]), (3, [(1, 2, "Z")])]
    horizontal, vertical = _both(persons)
    assert vertical == horizontal
    assert vertical[0] == ["2", "10", "3"]


def test_pvar_shifts_the_times_as_in_the_horizontal_route() -> None:
    persons = [("p1", [(1990, 1993, "A")]), ("p2", [(1991, 1994, "B")])]
    horizontal, vertical = _both(persons, {"origin": [1990, 1992]}, pvar="origin")
    assert vertical == horizontal


@pytest.mark.parametrize("priority", ["later", "earlier", "longer", "shorter"])
def test_overlap_priorities(priority: str) -> None:
    persons = [("p1", [(1, 4, "A"), (2, 3, "B")]), ("p2", [(2, 4, "A"), (1, 3, "B")])]
    horizontal, vertical = _both(persons, overlap_priority=priority)
    assert vertical == horizontal


def test_fillblanks() -> None:
    persons = [("p1", [(1, 2, "A"), (5, 6, "A")]), ("p2", [(2, 3, "B")])]
    horizontal, vertical = _both(persons, fillblanks="Z")
    assert vertical == horizontal


def test_per_person_bounds() -> None:
    persons = [("p1", [(1, 4, "A")]), ("p2", [(1, 4, "B")])]
    horizontal, vertical = _both(
        persons, {"lower": [2, 1], "upper": [3, 4]}, tmin_col="lower", tmax_col="upper"
    )
    assert vertical == horizontal


def test_a_per_person_column_that_varies_within_a_person_is_refused() -> None:
    vertical = pd.DataFrame(
        {"id": ["p1", "p1"], "start": [1, 3], "end": [2, 4], "state": ["A", "B"], "eo": [3, 4]}
    )
    with pytest.raises(ValueError, match="eo"):
        tx.HSPELL_to_STS(
            vertical, id_col="id", start_col="start", end_col="end", state_col="state", endObs="eo"
        )


_COLLIDING = {"endObs": [2004, 2004, 2003], "pvar": [1999, 1999, 2000]}


@pytest.mark.parametrize("argument", sorted(_COLLIDING))
@pytest.mark.parametrize("name", ["id", "b1", "e1", "s1"])
def test_a_person_column_may_share_a_name_with_the_reshaped_columns(
    argument: str, name: str
) -> None:
    # The reshape builds columns id, b<k>, e<k>, s<k>; a person-level column of
    # the same name used to overwrite one of them (batch-3 review M2).
    def build(column: str) -> tx.STSList:
        vertical = pd.DataFrame(
            {
                "pid": ["p1", "p1", "p2"],
                "start": [2000, 2003, 2001],
                "end": [2001, 2004, 2003],
                "state": ["X", "Y", "Z"],
                column: _COLLIDING[argument],
            }
        )
        return tx.HSPELL_to_STS(
            vertical,
            id_col="pid",
            start_col="start",
            end_col="end",
            state_col="state",
            **{argument: column},
        )

    assert _view(build(name)) == _view(build("person_level"))
