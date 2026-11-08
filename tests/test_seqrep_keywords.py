"""seqrep honours random_state and refuses keywords it cannot use (TX-AUD-20261010-015).

R's seqrep passes ``...`` to seqdist only when it computes the distances
(seqrep.R:24), so with ``diss`` given they are dropped; R has no random_state.
traminex's seqrep dropped every extra keyword with ``diss`` given, including
``random_state``, which dissrep takes: the "random" criterion then followed the
global NumPy stream, and a misspelled option passed silently, also through
seqrplot. seqrep now forwards ``random_state`` and refuses unused keywords.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import traminex as tx  # noqa: E402

ROWS = [list("AABB"), list("ABAB"), list("BBAA"), list("AAAA"), list("ABBB"), list("BABA")] * 3
SEQS = tx.seqdef(np.array(ROWS, dtype=object), alphabet=["A", "B"])
DISS = np.asarray(tx.seqdist(SEQS, method="LCS"), dtype=float)


def test_random_state_makes_random_representatives_reproducible() -> None:
    np.random.seed(1)
    first = tx.seqrep(SEQS, diss=DISS, criterion="random", random_state=7)["rep_indices"]
    np.random.seed(2)
    second = tx.seqrep(SEQS, diss=DISS, criterion="random", random_state=7)["rep_indices"]
    assert list(first) == list(second)


def test_random_state_is_forwarded_by_seqrep_grp() -> None:
    group = ["g", "h"] * 9
    np.random.seed(1)
    first = tx.seqrep_grp(SEQS, group, diss=DISS, criterion="random", ret="rep", random_state=3)
    np.random.seed(2)
    second = tx.seqrep_grp(SEQS, group, diss=DISS, criterion="random", ret="rep", random_state=3)
    assert {k: list(v["rep_indices"]) for k, v in first.items()} == {
        k: list(v["rep_indices"]) for k, v in second.items()
    }


def test_unknown_keyword_with_diss_is_refused() -> None:
    with pytest.raises(TypeError, match="bogus_option"):
        tx.seqrep(SEQS, diss=DISS, bogus_option=1)


def test_distance_keyword_with_diss_is_refused() -> None:
    # With a supplied matrix, a distance option has nothing to act on.
    with pytest.raises(TypeError, match="method"):
        tx.seqrep(SEQS, diss=DISS, method="OM")


def test_distance_keywords_without_diss_still_reach_seqdist() -> None:
    computed = tx.seqrep(SEQS, method="LCS")
    supplied = tx.seqrep(SEQS, diss=DISS)
    assert list(computed["rep_indices"]) == list(supplied["rep_indices"])


def test_seqrplot_refuses_unknown_keywords() -> None:
    with pytest.raises(TypeError, match="nonsense_option"):
        tx.seqrplot(SEQS, diss=DISS, nonsense_option=1)
