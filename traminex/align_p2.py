"""Priority 2: Sequence alignment and relative frequency analysis."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .sequences import STSList


class SeqAlign:
    """Result of pairwise sequence alignment.

    Attributes
    ----------
    operation : np.ndarray
        Array of operation codes: 'E' (equal), 'S' (substitution),
        'I' (insertion), 'D' (deletion)
    seq1 : np.ndarray
        Aligned sequence 1 (with '-' for gaps)
    seq2 : np.ndarray
        Aligned sequence 2 (with '-' for gaps)
    cost : np.ndarray
        Cost of each operation
    total_cost : float
        Total alignment cost (OM distance)
    cost_matrix : np.ndarray
        Full DP cost matrix
    operation_matrix : np.ndarray
        Full DP operation tracking matrix
    """

    def __init__(
        self,
        operation: np.ndarray,
        seq1: np.ndarray,
        seq2: np.ndarray,
        cost: np.ndarray,
        cost_matrix: np.ndarray,
        operation_matrix: np.ndarray,
    ):
        self.operation = operation
        self.seq1 = seq1
        self.seq2 = seq2
        self.cost = cost
        self.total_cost = float(cost.sum())
        self.cost_matrix = cost_matrix
        self.operation_matrix = operation_matrix

    def __repr__(self) -> str:
        """Pretty print alignment."""
        lines = []
        lines.append("Sequence Alignment")
        lines.append("=" * 60)

        # Build alignment display
        max_len = max(1, min(len(self.seq1), 80))  # Limit display width
        for i in range(0, len(self.seq1), max_len):
            chunk = slice(i, min(i + max_len, len(self.seq1)))
            lines.append(f"Seq 1:  {' '.join(self.seq1[chunk])}")
            lines.append(f"Op:     {' '.join(self.operation[chunk])}")
            lines.append(f"Cost:   {' '.join(f'{c:.1f}' for c in self.cost[chunk])}")
            lines.append(f"Seq 2:  {' '.join(self.seq2[chunk])}")
            lines.append("")

        lines.append(f"OM distance: {self.total_cost:.3f}")
        return "\n".join(lines)


def seqalign(
    x: STSList,
    indices: Tuple[int, int],
    *,
    indel: float = 1.0,
    sm: Optional[np.ndarray] = None,
    include_missing: bool = False,
) -> SeqAlign:
    """Compute pairwise sequence alignment with operation tracking.

    Performs optimal matching alignment between two sequences, tracking
    all operations (substitutions, insertions, deletions) and their costs.

    Parameters
    ----------
    x : STSList
        Sequence object containing the sequences to align
    indices : tuple of two ints
        0-based indices of the two sequences to align (Python or NumPy integers;
        booleans and fractions are refused)
    indel : float, default 1.0
        Insertion/deletion cost
    sm : np.ndarray, optional
        Substitution cost matrix over the K real states (K×K), or over the
        real states plus the missing state ((K+1)×(K+1)) when
        ``include_missing=True``, even if no cell is missing -- the size
        ``seqsubm(..., with_missing=True)`` builds. If None, uses constant cost
        of 2.0.
    include_missing : bool, default False
        If True, the missing state is an ordinary state (R's ``with.missing``);
        structural void is never a state.

    Returns
    -------
    SeqAlign
        Alignment result object with operations, costs, and matrices

    Examples
    --------
    >>> import numpy as np
    >>> from traminex import seqdef, seqalign, seqsubm
    >>> data = np.array([
    ...     ['A', 'B', 'C'],
    ...     ['A', 'C', 'C'],
    ... ])
    >>> sts = seqdef(data)
    >>> # Align first two sequences
    >>> sm = seqsubm(sts, method='CONSTANT', cval=2.0)
    >>> alignment = seqalign(sts, (0, 1), indel=1.0, sm=sm)
    >>> print(alignment.total_cost)  # OM distance
    2.0
    >>> print(alignment.operation)  # ['E', 'S', 'E']
    ['E' 'S' 'E']
    """
    # Two zero-based positions; R's one-based indices are numbers it coerces
    # (TRUE is 1, 1.7 is 1), which have no counterpart here (TX-AUD-20261011-074).
    # Any non-text sequence of two positions, as the earlier unpacking took
    # (range, a pandas Series; batch-15 review M-3).
    try:
        positions = () if isinstance(indices, (str, bytes)) else tuple(indices)
    except TypeError:
        positions = ()
    if len(positions) != 2 or not all(
        isinstance(i, (int, np.integer)) and not isinstance(i, (bool, np.bool_)) for i in positions
    ):
        raise TypeError(f"indices must be two integer positions, got {indices!r}")
    i1, i2 = int(positions[0]), int(positions[1])
    if i1 < 0 or i1 >= x.nseq or i2 < 0 or i2 >= x.nseq:
        raise ValueError(f"Invalid indices: {indices}")

    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects); fold declared external nr/void codes first.
    x = x._canonicalized()

    # Extract sequences. R aligns over seqlength() positions: trailing structural
    # void is dropped (it is right-padding, never a state), while the missing
    # state nr is kept as an ordinary state.
    seq1_codes = x.data[i1, :]
    seq2_codes = x.data[i2, :]
    if x.void_code is not None:
        seq1_codes = seq1_codes[seq1_codes != x.void_code]
        seq2_codes = seq2_codes[seq2_codes != x.void_code]
    l1 = len(seq1_codes)
    l2 = len(seq2_codes)

    # Build the cost alphabet. R's sm spans the real states, plus nr iff
    # with.missing=TRUE; structural void is never a state. traminex's
    # ``x.alphabet`` may carry the nr/void symbols (when they occur in the
    # data), so strip them and re-append nr last, exactly R's sm layout. An
    # explicit code→row map replaces raw-code indexing: alphabet position and
    # sm row coincide only for the marker-free prefix, so indexing sm by a
    # raw data code mis-addressed costs whenever nr sat anywhere but last.
    kept_codes = [
        code for code, _ in enumerate(x.alphabet) if code != x.nr_code and code != x.void_code
    ]
    if include_missing and x.nr_code is not None:
        kept_codes.append(x.nr_code)
    alph = [x.alphabet[code] for code in kept_codes]
    code_to_row = {code: row for row, code in enumerate(kept_codes)}
    if include_missing and x.nr_code is None:
        # R appends nr to sm's names under with.missing even when no cell is
        # missing (seqalign.R:10-14): sm keeps a row no data cell uses.
        alph.append(x.nr)

    def translate(codes: np.ndarray, row_index: int) -> np.ndarray:
        translated = np.empty(len(codes), dtype=int)
        for position, code in enumerate(codes):
            row = code_to_row.get(int(code))
            if row is None:
                # The only unmapped code after canonicalization and void
                # removal is the missing state under include_missing=False —
                # where R's seqalign stops too (subscript out of bounds).
                raise ValueError(
                    f"sequence {row_index} contains the missing state "
                    f"{x.nr!r}; pass include_missing=True to align over it"
                )
            translated[position] = row
        return translated

    seq1_codes = translate(seq1_codes, i1)
    seq2_codes = translate(seq2_codes, i2)

    K = len(alph)

    # Build substitution cost matrix
    if sm is None:
        # Default: constant cost of 2.0
        sm_array = np.full((K, K), 2.0, dtype=float)
        np.fill_diagonal(sm_array, 0.0)
    else:
        sm_array = np.asarray(sm, dtype=float)
        if sm_array.shape != (K, K):
            raise ValueError(f"sm must be {K}×{K} matrix")

    # Initialize DP matrices
    leven = np.zeros((l1 + 1, l2 + 1), dtype=float)
    operation = np.empty((l1 + 1, l2 + 1), dtype="<U1")

    # Initialize first row and column
    leven[:, 0] = np.arange(l1 + 1) * indel
    leven[0, :] = np.arange(l2 + 1) * indel
    # R initialises the first column and *then* the first row
    # (SOURCES/TraMineR/R/seqalign.R:22-23), so the origin cell ends up "I"
    # rather than empty. traminex overwrote it with "" -- an unrecorded
    # divergence in a public result object, on a cell no traceback visits
    # (cycle-6 section 9). Assignment order now reproduces R exactly.
    operation[:, 0] = "D"
    operation[0, :] = "I"

    # Fill DP matrix
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            c1 = int(seq1_codes[i - 1])
            c2 = int(seq2_codes[j - 1])

            if c1 == c2:
                subst_cost = 0.0
            else:
                subst_cost = sm_array[c1, c2]

            # Three options: substitution, insertion, deletion
            costs = np.array(
                [
                    leven[i - 1, j - 1] + subst_cost,  # Substitution/match
                    leven[i, j - 1] + indel,  # Insertion
                    leven[i - 1, j] + indel,  # Deletion
                ]
            )

            min_idx = np.argmin(costs)
            leven[i, j] = costs[min_idx]

            if min_idx == 0:
                operation[i, j] = "E" if subst_cost == 0 else "S"
            elif min_idx == 1:
                operation[i, j] = "I"
            else:
                operation[i, j] = "D"

    # Backtrack to get alignment
    i = l1
    j = l2
    op_list = []
    seq1_list = []
    seq2_list = []
    cost_list = []

    while i > 0 or j > 0:
        op = operation[i, j]
        old_cost = leven[i, j]

        if op == "S" or op == "E":
            seq1_list.append(alph[seq1_codes[i - 1]])
            seq2_list.append(alph[seq2_codes[j - 1]])
            i -= 1
            j -= 1
        elif op == "I":
            seq1_list.append("-")
            seq2_list.append(alph[seq2_codes[j - 1]])
            j -= 1
        elif op == "D":
            seq1_list.append(alph[seq1_codes[i - 1]])
            seq2_list.append("-")
            i -= 1

        op_list.append(op)
        new_cost = leven[i, j] if i > 0 or j > 0 else 0.0
        cost_list.append(old_cost - new_cost)

    # Reverse lists (backtracking goes backwards)
    op_list.reverse()
    seq1_list.reverse()
    seq2_list.reverse()
    cost_list.reverse()

    return SeqAlign(
        operation=np.array(op_list, dtype="<U1"),
        seq1=np.array(seq1_list, dtype=str),
        seq2=np.array(seq2_list, dtype=str),
        cost=np.array(cost_list, dtype=float),
        cost_matrix=leven,
        operation_matrix=operation,
    )


# The relative-frequency analysis lives in repseq_rf (the R-faithful dissrf /
# seqrf, RP4 T13b). Re-exported here for backward compatibility with the public
# `traminex.align_p2.dissrf` path. The previous simplified dissrf is removed.
from .repseq_rf import dissrf  # noqa: E402,F401
