"""``STSList`` checks its time labels and its alphabet (audit S2-15).

The constructor refused ``labels`` of the wrong length but accepted an
``xtlab`` that does not match the number of positions (time labels silently
shifted against the columns) and an alphabet naming the same state twice (two
codes that decode to the same text, so no consumer can tell them apart).

Source of truth: TraMineR's ``seqdef`` refuses a duplicated alphabet
("duplicated element(s) in alphabet", ``seqdef.R:68``); traminex's ``seqdef``
already refuses ``cnames`` of the wrong length, and the constructor now
applies the same rule (R's ``colnames<-`` pads a short ``cnames`` with NA).

Not changed: a ``cpal`` shorter than the alphabet stays accepted. It is a
supported case -- a palette for the real states only, as R's ``cpal`` has --
that ``visualization.palette.ensure_cmap`` pads with a warning
(``tests/test_audit_fixes_viz_formats.py``), and fixtures such as
``tests/test_r_goldens_seqgen_missing.py`` rely on it.
"""

from __future__ import annotations

import numpy as np
import pytest

from traminex import STSList

DATA = np.array([[0, 1, 1], [1, 0, 0]], dtype=int)


def test_wrong_length_xtlab_is_refused():
    with pytest.raises(ValueError, match="xtlab must have length 3"):
        STSList(DATA, alphabet=["A", "B"], xtlab=["t1", "t2"])


def test_duplicated_alphabet_is_refused():
    with pytest.raises(ValueError, match="duplicated element\\(s\\) in alphabet: A"):
        STSList(DATA, alphabet=["A", "A"])


def test_valid_metadata_is_still_accepted():
    x = STSList(
        np.array([[0, 1, 3], [1, 0, 2]], dtype=int),
        alphabet=["A", "B"],
        labels=["a", "b"],
        cpal=["#111111", "#222222"],
        xtlab=["t1", "t2", "t3"],
        nr_code=2,
        void_code=3,
    )
    assert x.cpal == ["#111111", "#222222"]
    assert x.xtlab == ["t1", "t2", "t3"]
    assert STSList(DATA, alphabet=["A", "B"]).cpal is None
    real_states_only = STSList(
        np.array([[0, 2]], dtype=int), alphabet=["A", "B", "*"], cpal=["#1", "#2"], nr_code=2
    )
    assert real_states_only.cpal == ["#1", "#2"]
