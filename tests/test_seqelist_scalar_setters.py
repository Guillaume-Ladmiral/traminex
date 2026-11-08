"""A single number sets the length or weight of a one-sequence list, as in R.

``seqelength<-`` and ``seqeweight<-`` compare ``length(value)`` with
``length(eseq)``, so R accepts a single number for a one-sequence list and
refuses it, by name, for a longer one. traminex's ``lengths``/``weights``
setters called ``len()`` on the 0-d array and raised
``TypeError: len() of unsized object`` in both cases (2026-09-22 audit S5-10).

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s5.R``
(TraMineR 2.2-14).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s5_event_objects.txt"
R = {
    fields[1]: fields[2:]
    for fields in (line.split("\t") for line in GOLDEN.read_text().splitlines())
    if fields[0] == "scalar_setters"
}


def test_a_single_number_sets_a_one_sequence_list():
    eseq = tx.seqecreate(id=[1, 1], timestamp=[0, 2], event=["a", "b"])
    eseq.lengths = 5
    eseq.weights = 2.5
    length, weight, text = R["one"]
    np.testing.assert_array_equal(tx.seqelength(eseq), [float(length)])
    np.testing.assert_array_equal(tx.seqeweight(eseq), [float(weight)])
    assert eseq.as_character() == [text]


@pytest.mark.parametrize("setter,value", [("length", 5), ("weight", 2.5)])
def test_a_single_number_is_refused_for_a_longer_list(setter, value):
    eseq = tx.seqecreate(id=[1, 2, 3], timestamp=[0, 0, 0], event=["a", "b", "a"])
    (message,) = R[setter]
    with pytest.raises(ValueError, match=re.escape(message.removeprefix("[!] "))):
        setattr(eseq, "lengths" if setter == "length" else "weights", value)
