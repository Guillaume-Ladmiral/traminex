"""Missing cells in the long-format and event builders (audit S2-05).

``seqdef_from_long`` -- and ``TSE_to_STS`` in state mode, which calls it --
pivoted the long table to wide form and turned every cell into text with
``astype(str)``. A period with no observation therefore became the real state
``"nan"`` (or ``"None"``/``"<NA>"`` for an explicit missing value), and in an
unbalanced panel the NaN the pivot inserts turned integer states into
``"1.0"``. ``TSE_to_STS`` in event mode did the same with the NA state an
``stm`` transition can produce.

Source of truth: TraMineR's ``seqdef`` defaults (``left = NA``, ``gaps = NA``,
``right = "DEL"``; ``SOURCES/TraMineR/R/seqdef.R:7``) make leading and
internal NA cells nonresponse and trailing ones void. R's ``TSE_to_STS``
returns the NA cells and ``seqdef`` decides them the same way. The R values are
in ``tests/r_goldens/audit20260922_s2_cases.txt``, written by
``scripts/r_goldens/make_goldens_audit20260922_s2.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import STSList, TSE_to_STS, alphabet, seqdef, seqdef_from_long

FIXTURE = Path(__file__).resolve().parent / "r_goldens" / "audit20260922_s2_cases.txt"


def _r_case(case: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        name, field, value = line.split("\t")
        if name == case:
            fields[field] = value
    assert fields, case
    return fields


def _r_cells(case: str) -> list[list[str]]:
    return [row.split("|") for row in _r_case(case)["cells"].split(";")]


def _decoded(x: STSList) -> list[list[str]]:
    return [[x._state_for_code(int(code)) for code in row] for row in x.data]


def _long(rows: list[tuple[int, int, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["id", "time", "state"])


# The R table ``wide`` of the generator: leading, internal and trailing NA and
# an unobserved row. Ids 1-3 omit their unobserved periods (an unbalanced
# panel); id 4 lists its three periods with explicit missing values.
NA_DEFAULTS_LONG = _long(
    [
        (1, 1, "A"),
        (1, 3, "B"),
        (2, 2, "A"),
        (2, 3, "B"),
        (3, 1, "A"),
        (3, 2, "A"),
        (4, 1, None),
        (4, 2, np.nan),
        (4, 3, pd.NA),
    ]
)


def test_seqdef_from_long_decides_missing_cells_as_seqdef_does():
    actual = seqdef_from_long(NA_DEFAULTS_LONG)

    assert _decoded(actual) == _r_cells("seqdef_na_defaults")
    assert alphabet(actual) == _r_case("seqdef_na_defaults")["alphabet"].split("|")
    assert actual.alphabet == ["A", "B", "*", "%"]
    assert (actual.nr_code, actual.void_code) == (2, 3)
    assert actual.ids == [1, 2, 3, 4]


def test_seqdef_from_long_matches_seqdef_on_the_equivalent_wide_table():
    wide = pd.DataFrame([["A", None, "B"], [None, "A", "B"], ["A", "A", None]])
    through_long = seqdef_from_long(NA_DEFAULTS_LONG.loc[NA_DEFAULTS_LONG["id"] < 4])
    direct = seqdef(wide)

    assert _decoded(through_long) == _decoded(direct)
    assert through_long.alphabet == direct.alphabet
    assert (through_long.nr_code, through_long.void_code) == (direct.nr_code, direct.void_code)


def test_unbalanced_integer_panel_keeps_integer_state_text():
    long = _long([(1, 1, 1), (1, 2, 2), (1, 3, 2), (2, 1, 1), (2, 2, 1)])

    actual = seqdef_from_long(long)

    assert _decoded(actual) == [["1", "2", "2"], ["1", "1", "%"]]
    assert alphabet(actual) == ["1", "2"]


def test_tse_to_sts_state_mode_decides_missing_cells_as_seqdef_does():
    actual = TSE_to_STS(NA_DEFAULTS_LONG)

    assert _decoded(actual) == _r_cells("seqdef_na_defaults")
    assert alphabet(actual) == ["A", "B"]


def _na_transition_call(**kwargs: object) -> STSList:
    stm = pd.DataFrame(
        {"e1": ["S1", pd.NA], "e2": ["S2", "S2"]},
        index=["None", "S1"],
        dtype=object,
    )
    tse = pd.DataFrame({"id": [1, 1, 2], "time": [1, 2, 1], "event": ["e1", "e1", "e2"]})
    return TSE_to_STS(tse, stm=stm, tmax=4, **kwargs)


def test_tse_to_sts_event_mode_na_transition_becomes_void_as_in_r():
    actual = _na_transition_call()

    assert _decoded(actual) == _r_cells("tse_na_transition")
    assert alphabet(actual) == _r_case("tse_na_transition")["alphabet"].split("|")
    assert actual.void_code is not None and actual.alphabet[actual.void_code] == "%"
    assert actual.nr_code is None


def test_event_mode_refuses_missing_cells_when_void_is_disabled():
    with pytest.raises(ValueError, match="void=None"):
        _na_transition_call(void=None)
