"""R-golden tests for ``seqe2stm`` -- the event -> state-machine table.

``TraMineRextras::seqe2stm`` builds the deterministic automaton mapping a set of
already-occurred events (the state) and an incoming event to the next state,
optionally with drop rules. We pin the reachable-state table (row/col names +
cells) against live TraMineRextras 0.6.9.

Generator: scripts/r_goldens/make_goldens_seqe2stm.R.
"""

from __future__ import annotations

import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

from traminex import seqe2stm


def _load_stm(name: str) -> pd.DataFrame:
    path = R_GOLDENS_DIR / f"{name}.tsv"
    lines = [
        ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")
    ]
    cols = lines[0].split("\t")
    idx: list[str] = []
    rows: list[list[str]] = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        idx.append(parts[0])
        rows.append(parts[1:])
    return pd.DataFrame(rows, index=idx, columns=cols)


def test_seqe2stm_abc_matches_r():
    pd.testing.assert_frame_equal(
        seqe2stm(["A", "B", "C"]),
        _load_stm("seqe2stm_abc"),
        check_exact=False,
        rtol=0,
        atol=R_GOLDEN_ATOL,
    )


def test_seqe2stm_droplist_matches_r():
    got = seqe2stm(["A", "B"], drop_list={"B": "A"})
    pd.testing.assert_frame_equal(
        got, _load_stm("seqe2stm_drop"), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqe2stm_first_state_matches_r():
    got = seqe2stm(["X", "Y"], first_state="start")
    pd.testing.assert_frame_equal(
        got, _load_stm("seqe2stm_first"), check_exact=False, rtol=0, atol=R_GOLDEN_ATOL
    )
