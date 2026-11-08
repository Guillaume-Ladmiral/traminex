"""R-golden tests for ``seqecontain`` -- per-sequence event-set containment.

``seqecontain(eseq, event_list)`` is TRUE for a sequence containing **any** of
the listed events; with ``unknown_exclude=True`` it is TRUE only if **all** of
the sequence's events are in the list. Pinned against live TraMineR 2.2-14.

Generator: scripts/r_goldens/make_goldens_eseq_contain.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import load_golden_lines

from traminex import seqecontain, seqecreate


def _eseq():
    # Mirrors scripts/r_goldens/make_goldens_eseq_contain.R.
    return seqecreate(
        pd.DataFrame(
            {
                "id": ["i1", "i1", "i1", "i2", "i2", "i3", "i3", "i4", "i4", "i5"],
                "timestamp": [0, 1, 2, 0, 1, 0, 1, 0, 1, 0],
                "event": ["A", "B", "C", "A", "B", "A", "A", "C", "C", "B"],
            }
        )
    )


def _golden_bool(name):
    return np.array([int(v) for v in load_golden_lines(name)], dtype=bool)


def test_seqecontain_any_single_event_matches_r():
    np.testing.assert_array_equal(seqecontain(_eseq(), ["A"]), _golden_bool("eseq_contain_anyA"))


def test_seqecontain_any_two_events_matches_r():
    np.testing.assert_array_equal(
        seqecontain(_eseq(), ["A", "C"]), _golden_bool("eseq_contain_anyAC")
    )


def test_seqecontain_all_in_set_matches_r():
    np.testing.assert_array_equal(
        seqecontain(_eseq(), ["A", "B"], unknown_exclude=True),
        _golden_bool("eseq_contain_allAB"),
    )


def test_seqecontain_accepts_single_string_and_scalar_eseq():
    e = _eseq()  # sequences in order i1, i2, i3, i4, i5 (indices 0..4)
    # A single label string works like a length-1 list.
    np.testing.assert_array_equal(seqecontain(e, "A"), _golden_bool("eseq_contain_anyA"))
    # A single Eseq returns a scalar bool.
    assert seqecontain(e[1], ["A"]) is True  # i2 contains A
    assert seqecontain(e[3], ["A"]) is False  # i4 is C only


def test_seqecontain_unknown_labels_ignored():
    e = _eseq()
    # "Z" is not in the dictionary -> ignored; matches plain "A".
    np.testing.assert_array_equal(seqecontain(e, ["A", "Z"]), _golden_bool("eseq_contain_anyA"))
