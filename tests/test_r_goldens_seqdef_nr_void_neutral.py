"""seqdef refuses nr == void and the NEUTRAL policies, recorded against R (TX-AUD-20261010-035).

R accepts ``nr`` equal to ``void`` and then handles a missing cell as void, so
the missing observation disappears (the row NA, A, B becomes A-B of length 2):
traminex refuses the ambiguous coding, a recorded correction. R's ``"NEUTRAL"``
left/gaps/right policies keep a plain NA cell that is neither missing nor void;
an STSList has no such cell, so traminex refuses them (owner decision of
2026-10-18: the refusal stays).
These refusals already existed: the tests characterize the recorded decisions.
Fixture: ``scripts/r_goldens/make_goldens_seqdef_nr_void_neutral.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqdef_nr_void_neutral.csv",
    keep_default_na=False,
)
DATA = np.array([["A", "B", None], [None, "A", "B"]], dtype=object)


def _r(case: str) -> list[tuple[str, int]]:
    rows = GOLDEN[GOLDEN["case"] == case].sort_values("sequence")
    return list(zip(rows["seqconc"], rows["seqlength"].astype(int)))


def test_nr_equal_to_void_is_refused_where_r_drops_the_missing_cell() -> None:
    with pytest.raises(ValueError, match="distinct"):
        tx.seqdef(DATA, nr="%")
    # R: the leading missing cell of row 2 became void, so three observed
    # periods (missing, A, B) read as a sequence of length 2.
    assert _r("nr_equals_void") == [("A-B", 2), ("A-B", 2)]


@pytest.mark.parametrize("policy", ["left", "right"])
def test_neutral_end_policies_are_refused_where_r_keeps_plain_na_cells(policy: str) -> None:
    with pytest.raises(ValueError, match="NEUTRAL"):
        tx.seqdef(DATA, **{policy: "NEUTRAL"})
    assert any(text.startswith("NA") or text.endswith("NA") for text, _ in _r(f"{policy}_neutral"))


def test_neutral_gaps_are_refused_where_r_keeps_a_plain_na_cell() -> None:
    with pytest.raises(ValueError, match="NEUTRAL"):
        tx.seqdef(np.array([["A", None, "B"]], dtype=object), gaps="NEUTRAL")
    assert _r("gaps_neutral") == [("A-NA-B", 3)]
