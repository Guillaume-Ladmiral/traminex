"""A real state spelled like a marker stays a state in the long builders (TX-AUD-20261019-005).

Only seqdef's STS route reads a marker spelling as the marker (seqdef ledger row; ruling R3-g
and its 2026-10-19 addendum). seqdef_from_long and TSE_to_STS (state mode) used to look the
nonresponse code up by the marker's text, so a literal ``*`` became the marker and merged with
genuinely missing cells. Expected values are counted by hand from the inputs.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

import traminex as tx


def _quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args, **kwargs)


FULL = pd.DataFrame(
    {"id": [1, 1, 1, 2, 2, 2], "time": [1, 2, 3, 1, 2, 3], "state": ["A", "*", "B", "A", "A", "B"]}
)
# Person 2 has no row at time 2: an internal unobserved period, so nonresponse.
GAP = pd.DataFrame(
    {"id": [1, 1, 1, 2, 2], "time": [1, 2, 3, 1, 3], "state": ["A", "*", "B", "A", "B"]}
)


@pytest.mark.parametrize("build", [tx.seqdef_from_long, tx.TSE_to_STS], ids=["long", "TSE_state"])
def test_literal_star_without_missing_cells_is_a_state(build) -> None:
    x = _quiet(build, FULL)
    # R's order of the real states: '*' sorts before letters in R's collation.
    assert tx.alphabet(x) == ["*", "A", "B"]
    assert x.nr_code is None
    # Every cell is a real state: both sequences have length 3 without missing.
    assert np.asarray(tx.seqlength(x, with_missing=False)).ravel().tolist() == [3, 3]


@pytest.mark.parametrize("build", [tx.seqdef_from_long, tx.TSE_to_STS], ids=["long", "TSE_state"])
def test_literal_star_and_a_missing_cell_keep_different_codes(build) -> None:
    x = _quiet(build, GAP)
    assert tx.alphabet(x) == ["*", "A", "B"]
    assert x.nr_code is not None
    assert x.alphabet[x.nr_code] != "*"
    literal, missing = int(x.data[0, 1]), int(x.data[1, 1])
    assert x.alphabet[literal] == "*"
    assert missing == x.nr_code
    assert literal != missing
    # Person 1: A, *, B (3 real states); person 2: A, missing, B (2 without missing).
    assert np.asarray(tx.seqlength(x, with_missing=False)).ravel().tolist() == [3, 2]
    assert np.asarray(tx.seqlength(x, with_missing=True)).ravel().tolist() == [3, 3]


def test_event_mode_keeps_a_literal_star_state() -> None:
    events = pd.DataFrame({"id": [1, 1], "time": [0, 1], "event": ["*", "B"]})
    x = _quiet(tx.TSE_to_STS, events, tmax=3)
    assert "*" in tx.alphabet(x)
    assert x.nr_code is None


EVENTS = pd.DataFrame({"id": [1, 1, 2], "time": [1, 2, 1], "event": ["A", "B", "A"]})
BUILDERS = {
    "long": lambda markers: tx.seqdef_from_long(GAP, **markers),
    "TSE_state": lambda markers: tx.TSE_to_STS(GAP, **markers),
    "TSE_event": lambda markers: tx.TSE_to_STS(EVENTS, state_col=None, tmax=3, **markers),
}


@pytest.mark.parametrize(
    "markers", [{"nr": "*", "void": "*"}, {"nr": "%"}], ids=["equal", "nr_is_default_void"]
)
@pytest.mark.parametrize("build", sorted(BUILDERS))
def test_equal_markers_are_refused_as_seqdef_refuses_them(build: str, markers: dict) -> None:
    # TX-REP-20261019-012: seqdef refuses nr == void (R reads a missing cell as void then;
    # seqdef ledger row, TX-AUD-20261010-035). The builders respelled one marker silently.
    message = "'nr' and 'void' must be distinct character markers"
    with pytest.raises(ValueError, match=message):
        tx.seqdef([["A", "B"]], **markers)
    with pytest.raises(ValueError, match=message):
        _quiet(BUILDERS[build], markers)


DISABLED = {
    "long": lambda: tx.seqdef_from_long(FULL, nr=None, void=None),
    "TSE_state": lambda: tx.TSE_to_STS(FULL, nr=None, void=None),
    "TSE_event": lambda: tx.TSE_to_STS(
        pd.DataFrame({"id": [1, 2], "time": [1, 1], "event": ["A", "B"]}),
        state_col=None,
        tmax=1,
        nr=None,
        void=None,
    ),
}


@pytest.mark.parametrize("build", sorted(DISABLED))
def test_two_disabled_markers_are_not_refused(build: str) -> None:
    # None disables a marker the table does not need; two disabled markers are not equal ones.
    x = _quiet(DISABLED[build])
    assert x.nr_code is None and x.void_code is None
