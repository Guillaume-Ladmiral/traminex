"""``disstree`` refuses two of R's weight permutations with an accurate reason.

Audit 2026-09-22, finding S4-13. ``disstree`` refused ``"permute"`` and
``"random-sampling"`` saying that R's ``disstree`` accepts "exactly" the five
other values. In fact TraMineR 2.2-14 accepts both: its split permutation test
(``SOURCES/TraMineR/R/disstree.dissassoc.R:52-64``) sends every value other
than ``"diss"``, ``"group"`` and ``"none"`` to the replicate branch, which
replicates cases by ``as.integer(weights)`` -- real weights truncated, so 0.5
counts as 0 -- and the final adjustment then runs ``dissassoc`` with that value.
An unknown value makes R stop at that adjustment. Verified live
(``reports/audit/full-codebase-audit-20260922/agents/S4/fix/s413_*``). The
refusal is kept, as a recorded divergence, with a message that says so.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from traminex import disstree, seqdef, seqtree

POS = np.array([0.0, 0.1, 0.2, 0.5, 0.6, 0.7, 5.0, 5.1, 5.2, 5.5, 5.6, 5.7])
D = np.abs(POS[:, None] - POS[None, :])
WEIGHTS = np.array([1.5, 0.5, 2.0, 1.0, 0.7, 1.3, 2.2, 0.9, 1.1, 1.6, 0.4, 1.8])
PREDICTORS = pd.DataFrame({"f": ["a"] * 6 + ["b"] * 6})


@pytest.mark.parametrize("method", ["permute", "random-sampling"])
def test_r_only_permutations_are_refused_with_the_reason(method):
    with pytest.raises(ValueError, match="truncated to integers") as refusal:
        disstree(D, PREDICTORS, weights=WEIGHTS, R=5, weight_permutation=method)
    assert "exactly" not in str(refusal.value)


def test_seqtree_gives_the_same_reason():
    # Without weights R and traminex both switch to "none" first, so the
    # refusal concerns weighted sequences only.
    sequences = seqdef(
        np.array([list("AAB"), list("ABB"), list("BBB"), list("BAA")], dtype=object),
        weights=[0.5, 1.5, 1.0, 2.0],
    )
    with pytest.raises(ValueError, match="truncated to integers"):
        seqtree(
            sequences,
            pd.DataFrame({"f": ["a", "a", "b", "b"]}),
            R=5,
            weight_permutation="permute",
        )


def test_unknown_permutation_is_refused():
    with pytest.raises(ValueError, match="unknown weight_permutation='foo'"):
        disstree(D, PREDICTORS, weights=WEIGHTS, R=5, weight_permutation="foo")
