"""HSPELL_to_STS blanks the periods from the end of observation on (TX-REP-20261017-005).

R's HSPELL_to_STS shifts the begin and end times by ``- tmin + 1``
(HSPELL_to_STS.R:25) but ``endObs`` only by ``- tmin`` (:67) before blanking
``seqindex >= endObs`` (:68), so it also blanks the period ``endObs - 1``,
which lies before the end of observation: a spell X 2000-2002 with endObs 2003
loses 2002. The consistent shift blanks the periods ``>= endObs``; traminex
corrects R here (ledger row HSPELL_to_STS). ``endObs`` is thus the first
unobserved period, by owner decision of 2026-10-18: each test here fails if
the period ``endObs`` is kept (blanking ``> endObs`` instead). R's own rows stay recorded in
``tests/r_goldens/hspell_horizontal_matrix_cases.txt``.
"""

from __future__ import annotations

import pandas as pd

import traminex as tx


def _as_r_frame(x: tx.STSList) -> list[list[str]]:
    markers = {x.nr_code, x.void_code} - {None}
    return [["NA" if code in markers else x.alphabet[code] for code in row] for row in x.data]


SPELLS = {"begin": ["b1", "b2"], "end": ["e1", "e2"], "status": ["s1", "s2"]}
PERSON = pd.DataFrame(
    {"id": ["p1"], "b1": [2000], "e1": [2002], "s1": ["X"], "b2": [2003], "e2": [2004], "s2": ["Y"]}
)


def test_a_period_before_the_end_of_observation_is_kept() -> None:
    result = tx.HSPELL_to_STS(PERSON, id_col="id", endObs=2003, **SPELLS)
    assert _as_r_frame(result) == [["X", "X", "X", "NA", "NA"]]


def test_the_end_of_observation_period_and_later_are_blanked() -> None:
    result = tx.HSPELL_to_STS(PERSON, id_col="id", endObs=2004, **SPELLS)
    assert _as_r_frame(result) == [["X", "X", "X", "Y", "NA"]]


def test_endobs_is_shifted_by_pvar_like_the_spell_times() -> None:
    frame = PERSON.assign(origin=[1999])
    result = tx.HSPELL_to_STS(frame, id_col="id", endObs=2003, pvar="origin", **SPELLS)
    assert _as_r_frame(result) == [["X", "X", "X", "NA", "NA"]]
