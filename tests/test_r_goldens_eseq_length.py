"""R-golden tests for ``seqelength`` (+ the ``seqelength<-`` setter) -- the
faithful, wrong-target-fixed event-sequence length.

R's ``seqelength`` returns the per-sequence *stored* length (observation time),
which is **-1** for a TSE-created eseq (NOT the number of events, which is what
traminex previously returned). The setter ``eseq.lengths = value`` mirrors
``seqelength<-`` and makes ``as.character`` show a trailing gap to the end of
observation. Pinned against live TraMineR 2.2-14.

Generator: scripts/r_goldens/make_goldens_eseq_length.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, load_golden_lines

from traminex import seqecreate, seqelength


def _eseq():
    # Mirrors scripts/r_goldens/make_goldens_eseq_length.R.
    return seqecreate(
        pd.DataFrame(
            {
                "id": ["i1", "i1", "i2", "i2", "i3"],
                "timestamp": [0, 2, 1, 1, 0],
                "event": ["A", "B", "A", "C", "B"],
            }
        )
    )


def _golden_floats(name):
    return np.array([float(v) for v in load_golden_lines(name)], dtype=float)


def test_seqelength_default_is_minus_one_like_r():
    # Faithful: stored length defaults to -1 for TSE eseq (not the event count).
    np.testing.assert_allclose(
        seqelength(_eseq()), _golden_floats("eseq_len_default"), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqelength_setter_matches_r():
    e = _eseq()
    e.lengths = [5, 4, 2]
    np.testing.assert_allclose(
        seqelength(e), _golden_floats("eseq_len_set"), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_as_character_with_set_length_matches_r():
    # A set length appears as a trailing gap to the end of observation.
    e = _eseq()
    e.lengths = [5, 4, 2]
    assert e.as_character() == load_golden_lines("eseq_len_aschar")


def test_seqelength_setter_length_check():
    e = _eseq()
    try:
        e.lengths = [1.0]  # wrong length (3 sequences)
    except ValueError:
        pass
    else:  # pragma: no cover
        raise AssertionError("expected ValueError for mismatched lengths")
