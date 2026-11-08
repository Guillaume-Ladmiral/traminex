"""``seqcost`` must not downgrade ``with_missing`` behind ``seqdist``'s back.

TraMineR 2.2-13 carried this block in ``seqcost``::

    if (isTRUE(with.missing) && !any(seqdata == attr(seqdata, "nr"))) {
        with.missing <- FALSE
        msg.warn("seqcost: 'with.missing' set as FALSE ...")
    }

``seqdist.R:112-113`` already decides ``with.missing`` across **both** the source
and an external ``refseq`` (``has.refseq.missings``), so a source without missing
cells and a reference *with* them is a legitimate widening case.  ``seqcost``
then silently undid that decision and returned an undersized ``sm``.  In R the
consequence was an out-of-bounds read: 2.2-13 returned
``5.09888213631624e-310, 5.09888181284790e-310, 2`` on one run and ``0, 0, 2`` on
the next for the same input.  TraMineR **2.2-14 removed the block**
(NEWS: "improved handling of missing and removed turning with.missing to FALSE
when no missing in seqdata, which caused an issue when refseq had missing
states" -- issue reported by Etienne Peyrot).

Every expected value below is live TraMineR 2.2-14's, recorded with ``%.17g``
by ``scripts/r_goldens/make_goldens_audit20260922_s10.R`` in
``tests/r_goldens/inline_r_anchors.txt``.
"""

from __future__ import annotations

import numpy as np
import pytest
from r_goldens_utils import r_anchor_shape, r_anchor_values

import traminex as tx

ALPHABET = ["A", "B", "C"]

#: Distances from each source row to the one-row reference, TraMineR 2.2-14.
R_2214_REFSEQ_DISTANCES = {
    method: r_anchor_values(f"refseq_missing_{method}")
    for method in ("CONSTANT", "TRATE", "INDELS", "INDELSLOG")
}


def _source() -> tx.STSList:
    """Three sequences, no missing cell anywhere."""
    return tx.seqdef(
        np.array([["A", "B", "B", "C"], ["A", "A", "C", "C"], ["B", "B", "A", "A"]], dtype=object),
        alphabet=ALPHABET,
    )


def _reference() -> tx.STSList:
    """One reference sequence that *does* carry a missing cell."""
    return tx.seqdef(np.array([["A", "B", None, "C"]], dtype=object), alphabet=ALPHABET)


def test_the_source_has_no_missing_and_the_reference_does():
    """Guard the fixture: the whole contract rests on this asymmetry."""
    source, reference = _source(), _reference()
    assert source.nr_code is None or not np.any(source.data == source.nr_code)
    assert reference.nr_code is not None and np.any(reference.data == reference.nr_code)


def test_seqcost_keeps_the_requested_missing_state(recwarn):
    """``with_missing=True`` must widen the matrix even on missing-free data.

    TraMineR 2.2-14 returns a 4x4 matrix over ``A,B,C,*`` with indel 1.
    """
    costs = tx.seqcost(_source(), method="CONSTANT", cval=2.0, with_missing=True)
    matrix = np.asarray(costs["sm"] if hasattr(costs, "keys") else costs, dtype=float)

    r_shape = r_anchor_shape("seqcost_constant_with_missing_sm")
    assert r_shape == (4, 4)
    assert matrix.shape == r_shape
    expected = r_anchor_values("seqcost_constant_with_missing_sm").reshape(r_shape, order="F")
    np.testing.assert_allclose(matrix, expected, rtol=0, atol=0)

    downgrades = [w for w in recwarn.list if "with_missing set to False" in str(w.message)]
    assert not downgrades, "seqcost must not downgrade with_missing (removed in 2.2-14)"


@pytest.mark.parametrize("method", sorted(R_2214_REFSEQ_DISTANCES))
def test_seqdist_with_a_missing_bearing_reference_matches_live_r(method):
    """All four string cost methods must reach R 2.2-14's values."""
    distances = tx.seqdist(
        _source(),
        method="OM",
        subm=method,
        indel=1.0,
        refseq=_reference(),
        with_missing=True,
    )
    np.testing.assert_allclose(
        np.asarray(distances, dtype=float).ravel(),
        R_2214_REFSEQ_DISTANCES[method],
        rtol=0,
        atol=1e-10,
    )
