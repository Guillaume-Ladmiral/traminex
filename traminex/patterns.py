"""Frequent state-subsequence mining -- a traminex-only extension.

There is **no TraMineR counterpart** to :func:`frequent_subsequences`, so it
cannot be validated by comparison with R. Its authority is the documented
contract below, pinned by ``tests/test_patterns_contract.py`` with expectations
worked out by hand on inputs small enough to check by inspection.

For the R-faithful event-subsequence miner, use
:func:`traminex.seqefsub` instead; that one is golden-pinned against live
TraMineR 2.2-14.
"""

from __future__ import annotations

from typing import List, Tuple

from .sequences import STSList


def frequent_subsequences(
    x: STSList,
    *,
    max_len: int = 3,
    contiguous: bool = True,
    min_support: int = 2,
    support_type: str = "occurrences",  # or "sequences"
) -> List[Tuple[Tuple[str, ...], int, float]]:
    """Mine frequent state subsequences of length 1..``max_len``.

    **traminex-only extension -- no TraMineR counterpart.** This docstring is
    the contract; there is no R oracle to appeal to.

    Contract
    --------
    *What is counted.* A pattern is a tuple of states read left to right.

    - ``contiguous=True`` (default): one occurrence per sliding window. For each
      sequence and each ``k`` in ``1..max_len`` every window ``row[t:t+k]``
      contributes one occurrence.
    - ``contiguous=False``: one occurrence per *combination of positions*. For
      each sequence and each ``k``, every strictly increasing k-tuple of
      positions contributes one occurrence of the states at those positions.
      The same pattern reached through different positions is counted once per
      position set, so a pattern can occur more times than the sequence is
      long. This counts position sets, not distinct embeddings of a pattern.

    *The denominator of ``support_fraction``.* ``count / denominator``, where

    - ``support_type="occurrences"`` and ``contiguous=True``:
      ``nseq * sum(max(seqlen - k + 1, 0) for k in 1..max_len)``, the total
      number of (sequence, k, window) slots;
    - ``support_type="occurrences"`` and ``contiguous=False``:
      ``nseq * sum(comb(seqlen, k) for k in 1..min(max_len, seqlen))``, the
      total number of position sets;
    - ``support_type="sequences"``: ``nseq``. The count is then the number of
      sequences containing the pattern **at least once**, so the fraction is
      the proportion of sequences carrying it and never exceeds 1.

    A zero denominator (an empty object) is replaced by 1 rather than raising.

    *Filtering and order.* A pattern is returned when its count is
    ``>= min_support`` -- an inclusive floor, so ``min_support=2`` keeps a
    pattern seen exactly twice. Results are sorted by count descending, ties
    broken by the state tuple ascending, so ``("A",)`` precedes ``("A", "B")``,
    which precedes ``("B",)``.

    *``max_len`` beyond the sequence length.* Lengths longer than ``seqlen``
    contribute nothing, to the numerator or the denominator.

    *Ragged sequences.* ``STSList`` stores a rectangular matrix, so shorter
    sequences are padded with the void marker (``"%"`` by default). This
    function treats that marker as an ordinary state: it appears in returned
    patterns, and the denominator uses the full matrix width for every
    sequence. Support fractions for a ragged object are therefore computed
    against more slots than the data actually occupies. Drop or recode the
    padding first if that is not what you want.

    Parameters
    ----------
    x : STSList
        State-sequence object.
    max_len : int, default 3
        Longest pattern to mine.
    contiguous : bool, default True
        Count sliding windows rather than position combinations.
    min_support : int, default 2
        Inclusive minimum count for a pattern to be returned.
    support_type : {"occurrences", "sequences"}, default "occurrences"
        Whether the count is total occurrences or number of sequences.

    Returns
    -------
    list of (tuple of str, int, float)
        ``(pattern, count, support_fraction)``, ordered as described above.
    """
    items: dict[Tuple[int, ...], int] = {}
    seq_support: dict[Tuple[int, ...], int] = {}
    for row in x.data:
        L = len(row)
        seen_this_seq: set[Tuple[int, ...]] = set()
        if contiguous:
            for k in range(1, max_len + 1):
                if k > L:
                    break
                for t in range(0, L - k + 1):
                    key = tuple(int(v) for v in row[t : t + k])
                    items[key] = items.get(key, 0) + 1
                    seen_this_seq.add(key)
        else:
            from itertools import combinations

            for k in range(1, max_len + 1):
                if k > L:
                    break
                for idxs in combinations(range(L), k):
                    key = tuple(int(row[i]) for i in idxs)
                    items[key] = items.get(key, 0) + 1
                    seen_this_seq.add(key)
        for key in seen_this_seq:
            seq_support[key] = seq_support.get(key, 0) + 1
    # Compute denominator for support fraction
    denom: float
    if support_type == "occurrences":
        if contiguous:
            denom = x.nseq * sum(max(x.seqlen - k + 1, 0) for k in range(1, max_len + 1))
        else:
            # approximate: number of k-combinations per row: sum C(L,k)
            from math import comb

            denom = x.nseq * sum(comb(x.seqlen, k) for k in range(1, min(max_len, x.seqlen) + 1))
    elif support_type == "sequences":
        denom = float(x.nseq)
    else:
        raise ValueError("support_type must be 'occurrences' or 'sequences'")
    out: List[Tuple[Tuple[str, ...], int, float]] = []
    keys = seq_support if support_type == "sequences" else items
    for key, cnt in keys.items():
        count = seq_support.get(key, 0) if support_type == "sequences" else cnt
        if count >= min_support:
            out.append(
                (
                    tuple(x._state_for_code(i) for i in key),
                    int(count),
                    float(count) / float(denom if denom > 0 else 1),
                )
            )
    out.sort(key=lambda z: (-z[1], z[0]))
    return out
