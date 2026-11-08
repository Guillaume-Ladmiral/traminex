"""R-golden tests for ``seqeage`` -- the event-age table (TraMineR's unexported
``seqeage`` over the C routine ``tmreventinseq``). Pinned against live TraMineR
TraMineR 2.2-14 on a small TSE with repeated events, simultaneous events, the same event
tied at one time, differing per-sequence counts, and a declared-but-absent level.

Generator: scripts/r_goldens/make_goldens_seqeage.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import assert_matches_golden

from traminex.events import seqeage, seqecreate


def _eseq():
    # Same inline TSE as make_goldens_seqeage.R. The event column is a pandas
    # Categorical with an extra level "D" that never occurs -- mirroring R's
    # factor(levels=c("A","B","C","D")) so the dictionary carries the absent
    # level (seqecreate honours a Categorical's declared level order).
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3", "s3", "s4", "s4"],
            "timestamp": [0, 2, 2, 5, 0, 1, 3, 0, 1, 1, 4, 0, 2],
            "event": pd.Categorical(
                ["A", "B", "C", "A", "A", "A", "C", "B", "A", "A", "C", "A", "B"],
                categories=["A", "B", "C", "D"],
            ),
        }
    )
    return seqecreate(tse)


def test_seqeage_dictionary_carries_absent_level():
    e = _eseq()
    assert e.dictionary == ["A", "B", "C", "D"]


def test_seqeage_A_repeated():
    # A occurs twice in s1/s2/s3, once in s4 -> width 2, -1 fill for s4 col 2.
    assert_matches_golden(seqeage(_eseq(), 1), "seqeage_A")


def test_seqeage_B_single_with_absent():
    # B occurs once in s1/s3/s4, never in s2 -> width 1, -1 for s2.
    assert_matches_golden(seqeage(_eseq(), 2), "seqeage_B")


def test_seqeage_C_single_with_absent():
    assert_matches_golden(seqeage(_eseq(), 3), "seqeage_C")


def test_seqeage_D_never_occurs_one_col_all_minus_one():
    out = seqeage(_eseq(), 4)
    assert out.shape == (4, 1)
    assert np.all(out == -1.0)
    assert_matches_golden(out, "seqeage_D")


def test_seqeage_event_list_uses_first_element():
    # The C routine reads only event.list[1]; seqeage(e, c(1,2)) == seqeage(e, 1).
    e = _eseq()
    np.testing.assert_array_equal(seqeage(e, [1, 2]), seqeage(e, 1))
    assert_matches_golden(seqeage(e, [1, 2]), "seqeage_list12")


def test_seqeage_single_eseq():
    # A single Eseq yields a (1, k) matrix; s1 (index 0) has A at times 0 and 5.
    e = _eseq()
    out = seqeage(e[0], 1)
    assert out.shape == (1, 2)
    np.testing.assert_array_equal(out, np.array([[0.0, 5.0]]))
