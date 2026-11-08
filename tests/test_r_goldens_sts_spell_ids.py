"""STS_to_SPELL labels the spells with the object's ids (TX-AUD-20261010-017).

R's STS_to_SPELL takes the ids from the sequence object's row names, which
seqdef(id=) sets (seqformat-STS_to_SPELL.R:31); traminex's STS_to_SPELL
numbered the spells 1..n although seqformat(to="SPELL") already used the ids.
Fixture: ``scripts/r_goldens/make_goldens_sts_spell_ids.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex as tx

# seqformat(to="SPELL") repeats R's message about the 'missing' codes; it is
# not the subject here.
pytestmark = pytest.mark.filterwarnings("ignore:'missing' set as:UserWarning")

GOLDEN = pd.read_csv(Path(__file__).resolve().parent / "r_goldens" / "sts_spell_ids.csv")
CASES = {
    "p1_p2": (np.array([list("AAB"), list("ABB")], dtype=object), ["p1", "p2"]),
    "z_x_y": (np.array([list("AB"), list("BB"), list("BA")], dtype=object), ["z", "x", "y"]),
}


def _r(case: str, route: str) -> pd.DataFrame:
    rows = GOLDEN[(GOLDEN["case"] == case) & (GOLDEN["route"] == route)]
    return rows[["id", "begin", "end", "states"]].reset_index(drop=True)


def _spells(frame: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id": [str(value) for value in frame["id"]],
            "begin": [int(value) for value in frame["begin"]],
            "end": [int(value) for value in frame["end"]],
            "states": [str(value) for value in frame["states"]],
        }
    )


@pytest.mark.parametrize("case", sorted(CASES))
def test_sts_to_spell_uses_the_object_ids(case: str) -> None:
    data, ids = CASES[case]
    spell = tx.STS_to_SPELL(tx.seqdef(data, id=ids))
    assert spell["id"].tolist() == _r(case, "STS_to_SPELL")["id"].tolist()
    pd.testing.assert_frame_equal(
        _spells(spell), _r(case, "STS_to_SPELL"), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )


@pytest.mark.parametrize("case", sorted(CASES))
def test_sts_to_spell_equals_seqformat_spell(case: str) -> None:
    data, ids = CASES[case]
    seqs = tx.seqdef(data, id=ids)
    pd.testing.assert_frame_equal(
        _spells(tx.STS_to_SPELL(seqs)),
        _spells(tx.seqformat(seqs, from_="STS", to="SPELL")),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    pd.testing.assert_frame_equal(
        _spells(tx.STS_to_SPELL(seqs)),
        _r(case, "seqformat"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def _trailing_nr():
    # R's seqdef(right = NA): the trailing missing cells stay nonresponse (TX-AUD-20261019-014).
    data = np.array([["A", "A", None, "B"], ["B", None, None, None]], dtype=object)
    return tx.seqdef(data, id=["q1", "q2"], right=None)


def test_sts_to_spell_keeps_trailing_missing_spells_as_r_s_internal_function() -> None:
    seqs = _trailing_nr()
    pd.testing.assert_frame_equal(
        _spells(tx.STS_to_SPELL(seqs)),
        _r("trailing_nr", "STS_to_SPELL"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
    pd.testing.assert_frame_equal(
        _spells(tx.STS_to_SPELL(seqs, with_missing=False)),
        _r("trailing_nr", "STS_to_SPELL_without_missing"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqformat_spell_still_drops_them_as_r_s_seqformat() -> None:
    seqs = _trailing_nr()
    pd.testing.assert_frame_equal(
        _spells(tx.seqformat(seqs, from_="STS", to="SPELL")),
        _r("trailing_nr", "seqformat"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )
