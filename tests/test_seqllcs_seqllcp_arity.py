"""``seqLLCS``/``seqLLCP`` are two-sequence functions and must refuse a container (F-021).

``SOURCES/TraMineR/R/seqLLCS.R`` and ``seqLLCP.R`` reach
``.C(C_cLCS, as.integer(seq1), as.integer(seq2), ...)``.  A multi-row ``stslist``
is a ``data.frame``, hence a list, which ``as.integer`` cannot coerce, so R
refuses with ``'list' object cannot be coerced to type 'integer'``. R's results
here are TraMineR 2.2-14's, recorded by
``scripts/r_goldens/make_goldens_audit20260922_s10.R``.

traminex previously narrowed to row 0 and returned a plausible scalar.  The
repaired sibling ``seqcomp`` already refuses the same shape through the
``exactly_one`` guard in ``traminex/stats_p2.py::_sequence_codes``; these two use
the identical message so one package does not answer the same question two ways.
"""

from __future__ import annotations

import numpy as np
import pytest
from r_goldens_utils import r_anchor, r_anchor_values

import traminex as tx

_RAW = np.array(
    [
        ["A", "B", "C", "A"],
        ["B", "A", "C", "C"],
    ],
    dtype=object,
)

ARITY_MESSAGE = "sequence objects must contain exactly one sequence"


def _two_row():
    return tx.seqdef(_RAW.copy(), alphabet=["A", "B", "C"])


def _one_row(index: int):
    x = _two_row()
    row = x.copy()
    row.data = np.asarray(x.data)[index : index + 1]
    return row


@pytest.mark.parametrize("function", [tx.seqLLCS, tx.seqLLCP])
def test_refuses_a_multi_sequence_object(function):
    r_case = "llcs_two_rows" if function is tx.seqLLCS else "llcp_two_rows"
    assert r_anchor(r_case, "error") == ["'list' object cannot be coerced to type 'integer'"]
    x = _two_row()
    with pytest.raises(ValueError, match=ARITY_MESSAGE):
        function(x, x)


@pytest.mark.parametrize("function", [tx.seqLLCS, tx.seqLLCP])
def test_refuses_a_multi_sequence_object_in_either_position(function):
    multi = _two_row()
    single = _one_row(0)
    with pytest.raises(ValueError, match=ARITY_MESSAGE):
        function(single, multi)
    with pytest.raises(ValueError, match=ARITY_MESSAGE):
        function(multi, single)


def test_single_row_objects_still_return_the_r_matching_scalar():
    """TraMineR 2.2-14 on this fixture: seqLLCS = 2, seqLLCP = 0."""
    first, second = _one_row(0), _one_row(1)
    assert tx.seqLLCS(first, second) == r_anchor_values("llcs_rows_1_2")[0] == 2
    assert tx.seqLLCP(first, second) == r_anchor_values("llcp_rows_1_2")[0] == 0
    # A sequence against itself: the LCS and the common prefix are the whole row.
    assert tx.seqLLCS(first, first) == r_anchor_values("llcs_row_1_itself")[0] == 4
    assert tx.seqLLCP(first, first) == r_anchor_values("llcp_row_1_itself")[0] == 4


def test_bare_integer_arrays_are_still_accepted_as_single_sequences():
    """The array route carries no container metadata, so 1-D input stays supported."""
    a = np.array([0, 1, 2, 0], dtype=np.int64)
    b = np.array([1, 0, 2, 2], dtype=np.int64)
    assert tx.seqLLCS(a, b) == 2
    assert tx.seqLLCP(a, b) == 0
