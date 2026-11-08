"""SPELL frames survive pandas operations that compare attrs (TX-REP-20261017-004).

STS_to_SPELL and seqformat(to="SPELL") stored the per-spell state codes as a
NumPy array in ``DataFrame.attrs``; pandas compares the attrs of the pieces it
combines (``astype`` with a dict, ``concat``), and comparing two arrays raised
"The truth value of an array with more than one element is ambiguous".
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import traminex as tx

# seqformat(to="SPELL") repeats R's message about the 'missing' codes; it is
# not the subject here.
pytestmark = pytest.mark.filterwarnings("ignore:'missing' set as:UserWarning")

SEQS = tx.seqdef(np.array([list("AAB"), list("ABB")], dtype=object), id=["p1", "p2"])
ROUTES = {
    "STS_to_SPELL": lambda: tx.STS_to_SPELL(SEQS),
    "seqformat": lambda: tx.seqformat(SEQS, from_="STS", to="SPELL"),
}


@pytest.mark.parametrize("route", sorted(ROUTES))
def test_astype_with_a_mapping_works(route: str) -> None:
    frame = ROUTES[route]()
    converted = frame.astype({"begin": float, "end": float})
    assert converted["begin"].tolist() == [1.0, 3.0, 1.0, 2.0]


@pytest.mark.parametrize("route", sorted(ROUTES))
def test_concat_of_two_spell_frames_works(route: str) -> None:
    # Two separately built frames: the same object twice compares equal by identity.
    first, second = ROUTES[route](), ROUTES[route]()
    assert len(pd.concat([first, second], ignore_index=True)) == 2 * len(first)


@pytest.mark.parametrize("route", sorted(ROUTES))
def test_state_codes_are_kept(route: str) -> None:
    assert list(ROUTES[route]().attrs["state_code"]) == [1, 2, 1, 2]
