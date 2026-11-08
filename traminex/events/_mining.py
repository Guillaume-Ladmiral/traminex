"""Faithful (R-parity) event-mining surface -- RP5b.

The R-faithful frequent-subsequence miner and the symbols built on TraMineR's
vendored C++ engine (:mod:`traminex.events._engine`): :func:`seqeconstraint` /
:class:`SeqEConstraint`, :func:`seqefsub` returning a :class:`SubSeqElist`
(the R ``subseqelist`` object), :func:`is_subseqelist`, and the functions built
on them: :func:`seqeapplysub`, :func:`seqerules`, :func:`seqecmpgroup`,
:func:`seqedist` and :func:`seqentrans`.

The numeric results are golden-pinned to live TraMineR 2.2-14, apart from the
corrections of R that the parity ledger records. This is a
BUILD-from-(GPL TraMineR's own C/R) port -- no Sequenzo borrow (its event core is
non-faithful).
"""

from __future__ import annotations

import copy
import math
import re
import warnings
from dataclasses import dataclass, replace
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from .._categorical_order import r_character_rank, stable_sorted_unique
from .._deprecated import MISSING, resolve_deprecated_alias
from .._rformat import r_as_character
from .._rsum import r_sum as _r_sum
from . import _engine
from ._object import SeqElist, is_seqelist

# count.method name -> integer code (seqeconstraint.R).
_COUNT_METHOD_NAMES = {"COBJ": 1, "CDIST_O": 2, "CWIN": 3, "CMINWIN": 4, "CDIST": 5}


@dataclass
class SeqEConstraint:
    """Event-subsequence mining constraints (R ``seqeconstraint``).

    Faithful to ``TraMineR::seqeconstraint`` -- the ``-1`` defaults mean "no
    constraint" (the C engine maps them to +/-infinity). ``count_method`` is the
    integer code 1..6 (``COBJ``=1, ``CDIST_O``=2, ``CWIN``=3, ``CMINWIN``=4,
    ``CDIST``=5; 6 = age, internal).
    """

    max_gap: float = -1
    window_size: float = -1
    age_min: float = -1
    age_max: float = -1
    age_max_end: float = -1
    count_method: int = 1


def seqeconstraint(
    max_gap: float = -1,
    window_size: float = -1,
    age_min: float = -1,
    age_max: float = -1,
    age_max_end: float = -1,
    count_method: Union[int, str] = 1,
    *,
    maxGap=MISSING,
    windowSize=MISSING,
    ageMin=MISSING,
    ageMax=MISSING,
    ageMaxEnd=MISSING,
    countMethod=MISSING,
) -> SeqEConstraint:
    """Define event-subsequence mining constraints (R ``seqeconstraint``).

    Faithful port of ``TraMineR::seqeconstraint`` (``seqeconstraint.R``): builds the
    constraint object consumed by :func:`seqefsub` / ``seqeapplysub``.

    Parameters
    ----------
    max_gap : float, default -1
        Maximum time gap between two consecutive events of a subsequence (``-1`` =
        no constraint).
    window_size : float, default -1
        Maximum total time span of a subsequence occurrence (``-1`` = no
        constraint).
    age_min, age_max, age_max_end : float, default -1
        Age (time from sequence start) constraints on the first / last event of an
        occurrence. If ``age_max_end`` is set and ``age_max`` is not, ``age_max``
        takes ``age_max_end``.
    count_method : int or str, default 1
        Occurrence-counting method: an integer 1..6 or one of the names
        ``"COBJ"`` (1), ``"CDIST_O"`` (2), ``"CWIN"`` (3), ``"CMINWIN"`` (4),
        ``"CDIST"`` (5). As in R, any other value, a fractional number such
        as 2.5 included, raises " [!] unknown count.method input".

    Returns
    -------
    SeqEConstraint
    """
    max_gap = resolve_deprecated_alias(
        max_gap,
        maxGap,
        canonical_name="max_gap",
        alias_name="maxGap",
        function_name="seqeconstraint",
        default=-1,
    )
    window_size = resolve_deprecated_alias(
        window_size,
        windowSize,
        canonical_name="window_size",
        alias_name="windowSize",
        function_name="seqeconstraint",
        default=-1,
    )
    age_min = resolve_deprecated_alias(
        age_min,
        ageMin,
        canonical_name="age_min",
        alias_name="ageMin",
        function_name="seqeconstraint",
        default=-1,
    )
    age_max = resolve_deprecated_alias(
        age_max,
        ageMax,
        canonical_name="age_max",
        alias_name="ageMax",
        function_name="seqeconstraint",
        default=-1,
    )
    age_max_end = resolve_deprecated_alias(
        age_max_end,
        ageMaxEnd,
        canonical_name="age_max_end",
        alias_name="ageMaxEnd",
        function_name="seqeconstraint",
        default=-1,
    )
    count_method = resolve_deprecated_alias(
        count_method,
        countMethod,
        canonical_name="count_method",
        alias_name="countMethod",
        function_name="seqeconstraint",
        default=1,
    )
    if age_max_end != -1 and age_max == -1:
        age_max = age_max_end
    if max_gap != -1 and window_size != -1 and max_gap > window_size:
        raise ValueError(" [!] max.gap is greater than window.size")
    if age_min != -1 and age_max != -1 and age_min > age_max:
        raise ValueError(" [!] age.min is greater than age.max or age.max.end")
    if isinstance(count_method, str):
        if count_method not in _COUNT_METHOD_NAMES:
            raise ValueError(" [!] unknown count.method input")
        cm = _COUNT_METHOD_NAMES[count_method]
    else:
        # R's count.method %in% 1:6 compares values, so 2.5 is refused rather than truncated (seqeconstraint.R:18-22).
        if count_method not in (1, 2, 3, 4, 5, 6):
            raise ValueError(" [!] unknown count.method input")
        cm = int(count_method)
    return SeqEConstraint(
        max_gap=max_gap,
        window_size=window_size,
        age_min=age_min,
        age_max=age_max,
        age_max_end=age_max_end,
        count_method=cm,
    )


class SubSeqElist:
    """A list of frequent event subsequences (R ``subseqelist``), from :func:`seqefsub`.

    Mirrors R's ``subseqelist``: ``subseq`` are the generic subsequence label
    strings (``"(A)"``, ``"(A)-(B)"``, ``"(A,B)"`` -- no gaps shown), ``data`` is a
    :class:`pandas.DataFrame` with ``Support`` (weighted count / total weight) and
    ``Count`` (weighted occurrence count) columns, ``constraint`` the mining
    constraint and ``eseq`` the source event-sequence list. Plot-ready for RP7.
    """

    def __init__(self, eseq, constraint, subseq_objs, data, type="frequent", count_total=None):
        self.eseq = eseq
        self.constraint = constraint
        self._subseq_objs = list(subseq_objs)
        self.subseq: List[str] = [s.sprint() for s in self._subseq_objs]
        self.data: pd.DataFrame = data
        self.type = type
        # The total the counts were taken on: the sum of the weights for weighted
        # mining, the number of sequences otherwise (seqentrans' Avg.occ needs it).
        self.count_total = count_total

    def __len__(self) -> int:
        return len(self.subseq)

    def as_character(self) -> List[str]:
        """The generic subsequence label strings (R ``as.character(x$subseq)``)."""
        return list(self.subseq)

    def __getitem__(self, idx):
        """An integer gives the label; a slice, positions or a mask give a subset.

        A subset is a list of the same class, as R's ``[.subseqelist``
        (subseqelist.R:28-45) returns: ``subseq``, the subsequence objects and
        the ``data`` rows are subset, everything else is carried (TX-AUD-20261010-038).
        An integer keeps Python's scalar-index convention and gives the label
        string, which also makes iteration yield the labels.
        """
        if isinstance(idx, (int, np.integer)) and not isinstance(idx, (bool, np.bool_)):
            return self.subseq[idx]
        size = len(self.subseq)
        if isinstance(idx, slice):
            positions = list(range(size))[idx]
        else:
            if isinstance(idx, (str, bytes)) or not hasattr(idx, "__len__"):
                raise TypeError(
                    "SubSeqElist indices must be an integer, a slice, a sequence of "
                    f"integers or a boolean mask, not {type(idx).__name__}"
                )
            selector = np.asarray(idx)
            if selector.dtype == bool:
                if selector.ndim != 1 or selector.size != size:
                    raise IndexError(f"boolean mask must have length {size}, got {selector.size}")
                positions = np.flatnonzero(selector).tolist()
            elif selector.size == 0 or np.issubdtype(selector.dtype, np.integer):
                if selector.ndim != 1:
                    raise IndexError("SubSeqElist positions must be one-dimensional")
                positions = [int(p) for p in selector]
                for p in positions:
                    if not -size <= p < size:
                        raise IndexError(f"position {p} is out of range for {size} subsequences")
                positions = [p % size for p in positions]
            else:
                raise TypeError(
                    "SubSeqElist positions must be integers or booleans, " f"not {selector.dtype}"
                )
        return self._subset(positions)

    def _subset(self, positions: List[int]) -> "SubSeqElist":
        subset = copy.copy(self)
        subset._subseq_objs = [self._subseq_objs[p] for p in positions]
        subset.subseq = [self.subseq[p] for p in positions]
        subset.data = self.data.iloc[positions]
        return subset

    def __repr__(self) -> str:
        return f"SubSeqElist(n={len(self.subseq)}, type={self.type!r})"


def is_subseqelist(x: object) -> bool:
    """Whether ``x`` is a frequent-subsequence list (R ``is.subseqelist``)."""
    return isinstance(x, SubSeqElist)


def seqefsub(
    eseq=MISSING,
    *,
    str_subseq: Optional[List[str]] = None,
    min_support: Optional[float] = None,
    pmin_support: Optional[float] = None,
    constraint: Optional[SeqEConstraint] = None,
    max_k: int = -1,
    weighted: bool = True,
    seq=MISSING,
    strsubseq=MISSING,
    minSupport=MISSING,
    pMinSupport=MISSING,
    maxK=MISSING,
) -> SubSeqElist:
    """Mine frequent event subsequences (R ``seqefsub``).

    Faithful port of ``TraMineR::seqefsub`` (``seqefsub.R`` + the C engine
    ``C_tmrfindsubsequences``): a prefix-tree miner that enumerates candidate
    subsequences over increasing length, prunes at ``min_support``, then recounts
    each surviving subsequence across all sequences with the constraint's count
    method. Returns a :class:`SubSeqElist` whose ``data`` carries ``Support`` =
    (weighted #sequences with a positive count) / total weight and ``Count`` =
    weighted occurrence count, ordered by ``Count`` descending.

    Parameters
    ----------
    eseq : SeqElist
        Event sequences from :func:`~traminex.events.seqecreate`.
    min_support : float, optional
        Absolute minimum support (weighted sequence count) for a subsequence to be
        kept. Provide either ``min_support`` or ``pmin_support``.
    pmin_support : float, optional
        Proportional minimum support; ``min_support = pmin_support * total weight``.
    constraint : SeqEConstraint, optional
        Mining constraint (default :func:`seqeconstraint`).
    max_k : int, default -1
        Maximum subsequence length (``-1`` = unbounded).
    weighted : bool, default True
        If ``True`` use the event-sequence weights; if ``False`` treat every
        sequence as weight 1 (``total weight = number of sequences``).
    str_subseq : list of str, optional
        Compute support/count for these *specific* subsequences (given as generic
        label strings, e.g. ``"(A)-(B,C)"``) instead of mining frequent ones. No
        mining happens; ``min_support``/``max_k`` are ignored and the subsequences
        are returned in the given order (R ``seqefsub(str.subseq=...)``).

    Returns
    -------
    SubSeqElist
    """
    eseq = resolve_deprecated_alias(
        eseq,
        seq,
        canonical_name="eseq",
        alias_name="seq",
        function_name="seqefsub",
    )
    str_subseq = resolve_deprecated_alias(
        str_subseq,
        strsubseq,
        canonical_name="str_subseq",
        alias_name="strsubseq",
        function_name="seqefsub",
        default=None,
    )
    min_support = resolve_deprecated_alias(
        min_support,
        minSupport,
        canonical_name="min_support",
        alias_name="minSupport",
        function_name="seqefsub",
        default=None,
    )
    pmin_support = resolve_deprecated_alias(
        pmin_support,
        pMinSupport,
        canonical_name="pmin_support",
        alias_name="pMinSupport",
        function_name="seqefsub",
        default=None,
    )
    max_k = resolve_deprecated_alias(
        max_k,
        maxK,
        canonical_name="max_k",
        alias_name="maxK",
        function_name="seqefsub",
        default=-1,
    )
    if not is_seqelist(eseq):
        raise TypeError(" [!] eseq should be a seqelist. See help on seqecreate.")
    if constraint is None:
        constraint = seqeconstraint()
    if constraint.count_method == 3 and constraint.window_size == -1:
        raise ValueError(" [!] count method 3 or 'CWIN' requires a window.size.")
    if constraint.count_method == 6:
        # seqefsub.R:26-28 (2026-09-22 audit S5-09); mining still proceeds.
        warnings.warn(" [!] count method 6 is only for internal use.", UserWarning, stacklevel=2)

    orig_w = [float(s.weight) for s in eseq.sequences]
    n = len(eseq.sequences)
    if weighted:
        use_w = orig_w
        totseq = float(sum(orig_w))
    else:
        use_w = [1.0] * n
        totseq = float(n)

    if str_subseq is not None:
        # Specific (not frequent) subsequences: no mining. Build them from their
        # label strings (seqecreatesub), then Count = sum_i w_i*count(cMethod) and
        # Support = (sum_i w_i*COBJ-count)/totseq, in the GIVEN order (no reorder).
        return _seqefsub_str(eseq, str_subseq, constraint, use_w, totseq)

    if min_support is None:
        if pmin_support is None:
            raise ValueError(
                " [!] You should specify a minimum support through min_support or "
                "pmin_support argument"
            )
        min_support = float(pmin_support) * float(sum(use_w))

    eng = [_engine.eseq_to_engine(s) for s in eseq.sequences]
    for s_eng, w in zip(eng, use_w):
        s_eng.weight = w

    supports, counts, subseqs = _engine.tmrfindsubsequences(
        eng,
        constraint.max_gap,
        constraint.window_size,
        constraint.age_min,
        constraint.age_max,
        constraint.age_max_end,
        constraint.count_method,
        float(min_support),
        int(max_k),
    )

    # Order by Count descending, ties keeping enumeration order (R order(decreasing=TRUE)).
    counts_arr = np.asarray(counts, dtype=float)
    order = np.argsort(-counts_arr, kind="stable")
    sub_ord = [subseqs[i] for i in order]
    support_ord = np.asarray(supports, dtype=float)[order]
    count_ord = counts_arr[order]
    data = pd.DataFrame({"Support": support_ord / totseq, "Count": count_ord})

    return SubSeqElist(eseq, constraint, sub_ord, data, type="frequent", count_total=totseq)


def _seqecreatesub(str_subseq: List[str], eseq) -> List["_engine._Sequence"]:
    """Parse generic subsequence label strings into engine sequences (R ``seqecreatesub``).

    A label like ``"(A)-(B,C)"`` is split into groups on ``)-(`` (each group is one
    transition); within a group, comma-separated tokens are simultaneous events. Each
    token maps to its 1-based code via ``which(token == levels(eseq))`` (an unknown or
    ambiguous token is an error), the group index is the timestamp, and the events are
    sorted by ``(time, code)`` -- exactly ``seqecreate.R::seqecreatesub``.
    """
    codebase = list(eseq.dictionary)
    objs: List[_engine._Sequence] = []
    for label in str_subseq:
        inner = re.sub(r"^\(|\)$", "", str(label).strip())
        groups = re.split(r"\)\s*-\s*\(", inner)
        timestamps: List[float] = []
        codes: List[int] = []
        tindex = 1
        for grp in groups:
            for tok in grp.split(","):
                tok = tok.strip()
                matches = [k for k, c in enumerate(codebase) if c == tok]
                if len(matches) != 1:
                    raise ValueError(
                        f" [!] '{tok}' not found in the alphabet of events: " + ",".join(codebase)
                    )
                timestamps.append(float(tindex))
                codes.append(matches[0] + 1)  # 1-based factor code
            tindex += 1
        order = sorted(range(len(codes)), key=lambda i: (timestamps[i], codes[i]))
        s = _engine._Sequence(-1, list(codebase), generic=True)
        prev_t = 0.0
        head: Optional[_engine._Node] = None
        tail: Optional[_engine._Node] = None
        for i in order:
            t = timestamps[i]
            node = _engine._Node(codes[i], t if tail is None else t - prev_t)
            if head is None:
                head = node
            else:
                assert tail is not None
                tail.next = node
            tail = node
            prev_t = t
        s.event = head
        objs.append(s)
    return objs


def _seqefsub_str(eseq, str_subseq, constraint, use_w, totseq) -> SubSeqElist:
    """seqefsub for user-specified subsequences (R ``seqefsub`` str.subseq branch)."""
    subseq_objs = _seqecreatesub(str_subseq, eseq)
    ret = SubSeqElist(eseq, constraint, subseq_objs, pd.DataFrame(), type="user")
    w = np.asarray(use_w, dtype=float)
    # Count under the user constraint's count method; Support numerator under COBJ (=1).
    m_count = seqeapplysub(ret, constraint=constraint).to_numpy()
    m_support = seqeapplysub(ret, constraint=replace(constraint, count_method=1)).to_numpy()
    count_col = (w[:, None] * m_count).sum(axis=0)
    support_col = (w[:, None] * m_support).sum(axis=0)
    data = pd.DataFrame({"Support": support_col / totseq, "Count": count_col})
    return SubSeqElist(eseq, constraint, subseq_objs, data, type="user", count_total=totseq)


# method name -> count.method integer (seqeapplysub.R).
_APPLY_METHOD_MAP = {
    "count": 2,
    "presence": 1,
    "age": 6,
    "COBJ": 1,
    "CDIST_O": 2,
    "CWIN": 3,
    "CMINWIN": 4,
    "CDIST": 5,
}


def seqeapplysub(
    subseq: SubSeqElist,
    method: Optional[Union[int, str]] = None,
    constraint: Optional[SeqEConstraint] = None,
    rules: bool = False,
) -> pd.DataFrame:
    """Count each subsequence in each event sequence (R ``seqeapplysub``).

    Faithful port of ``TraMineR::seqeapplysub`` (the C engine
    ``C_tmrmatrixsubseqinseq``): returns the ``(n_seq x n_subseq)`` matrix whose
    ``[i, j]`` cell is the number of occurrences of subsequence ``j`` in sequence
    ``i`` under the chosen count method (rows = the event sequences' label strings,
    columns = the generic subsequence labels).

    Parameters
    ----------
    subseq : SubSeqElist
        A :func:`seqefsub` result.
    method : int or str, optional
        Count method override: ``"count"`` (2), ``"presence"`` (1), ``"age"`` (6),
        or a method name/number ``"COBJ"``/1 .. ``"CDIST"``/5. Defaults to the
        subsequence list's own constraint. As in R, any other value, such as 2.5
        or 6, leaves the constraint's count method; a number written as text,
        such as ``"2"``, counts as that number.
    constraint : SeqEConstraint, optional
        Constraint to count under (defaults to ``subseq.constraint``).
    rules : bool, default False
        If ``True`` count subsequences within the subsequences themselves (the
        ``(n_subseq x n_subseq)`` matrix used by :func:`seqerules`).

    Returns
    -------
    pandas.DataFrame
        The count matrix, indexed by sequence labels with subsequence-label columns.
    """
    if not is_subseqelist(subseq):
        raise TypeError("subseq should be a subseqelist. See help on seqefsub.")
    if constraint is None:
        constraint = subseq.constraint
    if method is not None:
        if isinstance(method, str) and method in _APPLY_METHOD_MAP:
            constraint = replace(constraint, count_method=_APPLY_METHOD_MAP[method])
        # R's method %in% 1:5 matches by value, text with text (seqeapplysub.R:28): 2.5 and "2.0" fall through, "2" is 2.
        elif method in (("1", "2", "3", "4", "5") if isinstance(method, str) else (1, 2, 3, 4, 5)):
            constraint = replace(constraint, count_method=int(method))
        # else: leave the constraint's count.method unchanged (matches seqeapplysub.R)
    if constraint.count_method == 3 and constraint.window_size == -1:
        # seqefsub's refusal (seqefsub.R:23-25). seqeapplysub.R has no check and
        # R aborts the session; the engine hit a NaN conversion (audit S5-11).
        raise ValueError(" [!] count method 3 or 'CWIN' requires a window.size.")

    sub_objs = subseq._subseq_objs
    if rules:
        targets = sub_objs
        row_labels = list(subseq.subseq)
    else:
        targets = [_engine.eseq_to_engine(s) for s in subseq.eseq.sequences]
        row_labels = subseq.eseq.as_character()

    matrix = _engine.tmrmatrixsubseqinseq(
        sub_objs,
        targets,
        constraint.max_gap,
        constraint.window_size,
        constraint.age_min,
        constraint.age_max,
        constraint.age_max_end,
        constraint.count_method,
    )
    return pd.DataFrame(matrix, index=row_labels, columns=list(subseq.subseq))


def _pnorm(z: float) -> float:
    """Standard normal CDF (R ``pnorm``); ``nan``/``+-inf`` propagate as in R."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _implicativestat(x, y, kind: str = "intensity", resid: str = "standard") -> np.ndarray:
    """Implicative statistic of a rule (R internal ``implicativestat``).

    Faithful port of ``TraMineR``'s internal ``implicativestat`` (``implicativestat.R``)
    used by :func:`seqerules`. NOTE: the R source carries a level-coercion quirk --
    ``if (length(levels(x)==1))`` is always truthy -- so the factors are *always*
    forced to levels ``c("0", "1")``; this is reproduced here (the inputs are 0/1
    presence indicators). Returns the ``2 x 2`` matrix over (x-level, y-level);
    ``seqerules`` reads cell ``[1, 1]`` (R ``[2, 2]``).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    levels = (0.0, 1.0)
    result = np.zeros((2, 2), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        for i, xv in enumerate(levels):
            condi = x == xv
            for j, yv in enumerate(levels):
                condj = y == yv
                nnbj = float(np.sum(condi & ~condj))
                nexpnbj = float(np.sum(condi)) * float(np.sum(~condj)) / n
                if resid == "standard":
                    indice = (nnbj - nexpnbj) / np.sqrt(nexpnbj)
                elif resid == "deviance":
                    indice = np.sign(nnbj - nexpnbj) * np.sqrt(
                        abs(2 * nnbj * np.log(nnbj / nexpnbj))
                    )
                elif resid == "Freeman-Tukey":
                    indice = np.sqrt(nnbj) + np.sqrt(1 + nnbj) - np.sqrt(4 * nexpnbj + 1)
                elif resid == "adjusted":
                    indice = (nnbj - nexpnbj) / np.sqrt(
                        nexpnbj * (np.sum(condj) / n) * (1 - np.sum(condi) / n)
                    )
                else:
                    raise ValueError(
                        "resid should be one of standard, deviance, Freeman-Tukey or adjusted"
                    )
                result[i, j] = indice if kind == "indice" else _pnorm(-float(indice))
    return result


def seqerules(fsubseq: SubSeqElist, sortv: Optional[str] = None, decreasing: bool = False):
    """Association rules from frequent event subsequences (R ``seqerules``).

    Faithful port of ``TraMineR::seqerules`` (``seqerules.R``): for every mined
    subsequence with at least two transitions, every prefix split ``A => B`` yields
    a rule with the full TraMineR statistics computed from the presence matrix
    (:func:`seqeapplysub` ``method="presence"``):

    - ``Support`` = number of sequences containing the whole subsequence,
    - ``Conf`` = P(AB)/P(A), ``Lift`` = P(AB)/(P(A)P(B)),
    - ``Standardlift`` = (Lift - lmb)/(ups - lmb) with ``lmb = max(P(A)+P(B)-1, 1/n)``
      and ``ups = 1/max(P(A), P(B))``,
    - ``JMeasure``, ``ImplicStat`` (``implicativestat(A, B)[2, 2]``),
    - ``p.value`` = ``pnorm(ImplicStat)``, ``p.valueB1`` / ``p.valueB2`` =
      Bonferroni corrections over the rule / discriminating-subsequence counts.

    Degenerate rules produce R's exact ``NaN``/``+-Inf`` IEEE results.

    Parameters
    ----------
    fsubseq : SubSeqElist
        A :func:`seqefsub` result.
    sortv : str, optional
        Sort the rules by one of the (R-internal) statistic names ``"support"``,
        ``"confidence"``, ``"lift"``, ``"slift"``, ``"jmeasure"``, ``"implicative"``,
        ``"implicativep"``, ``"implicativepb1"``, ``"implicativepb2"`` or ``"rule"``.
    decreasing : bool, default False
        Sort direction when ``sortv`` is given.

    Returns
    -------
    pandas.DataFrame
        Columns ``Rules, Support, Conf, Lift, Standardlift, JMeasure, ImplicStat,
        p.value, p.valueB1, p.valueB2``.
    """
    if not is_subseqelist(fsubseq):
        raise TypeError("fsubseq should be a subseqelist. See help on seqefsub.")

    matcount = seqeapplysub(fsubseq, method="presence")
    # Column-name -> summed presence vector (and the raw 0/1 vector, for implicativestat).
    col_sum = {lab: matcount[lab].to_numpy() for lab in matcount.columns}
    ntotal = len(fsubseq.eseq)

    cols = ["Rules", "Support", "Conf", "Lift", "Standardlift", "JMeasure", "ImplicStat"]
    rows: List[dict] = []
    pvals: List[float] = []
    nsubseq = 0
    with np.errstate(divide="ignore", invalid="ignore"):
        for label in fsubseq.subseq:
            splitted = label.split("-")
            if len(splitted) <= 1:
                continue
            nsubseq += 1
            for z in range(1, len(splitted)):
                a = "-".join(splitted[:z])
                b = "-".join(splitted[z:])
                counta = float(np.sum(col_sum[a]))
                countb = float(np.sum(col_sum[b]))
                countab = float(np.sum(col_sum[label]))
                pa = np.float64(counta) / ntotal
                pb = np.float64(countb) / ntotal
                pab = np.float64(countab) / ntotal
                confidence = pab / pa
                lift = pab / (pa * pb)
                lmb = max(pa + pb - 1.0, np.float64(1.0) / ntotal)
                ups = np.float64(1.0) / max(pa, pb)
                slift = (lift - lmb) / (ups - lmb)
                invconf = 1.0 - confidence
                invpb = 1.0 - pb
                jmeasure = confidence * np.log2(confidence / pb) + invconf * np.log2(
                    invconf / invpb
                )
                icstat = float(_implicativestat(col_sum[a], col_sum[b], kind="indice")[1, 1])
                icstatp = _pnorm(icstat)
                rows.append(
                    {
                        "Rules": f"{a} => {b}",
                        "Support": float(countab),
                        "Conf": float(confidence),
                        "Lift": float(lift),
                        "Standardlift": float(slift),
                        "JMeasure": float(jmeasure),
                        "ImplicStat": icstat,
                    }
                )
                pvals.append(icstatp)

    # R quirk: nrule ends at (#rules + 1) after the loop; p.valueB2 uses nsubseq.
    nrule = len(rows) + 1
    for r, p in zip(rows, pvals):
        r["p.value"] = float(p)
        r["p.valueB1"] = float(1.0 - (1.0 - p) ** nrule)
        r["p.valueB2"] = float(1.0 - (1.0 - p) ** nsubseq)

    all_cols = cols + ["p.value", "p.valueB1", "p.valueB2"]
    out = pd.DataFrame(rows, columns=all_cols)

    if sortv is not None:
        sort_map = {
            "rule": "Rules",
            "support": "Support",
            "confidence": "Conf",
            "lift": "Lift",
            "slift": "Standardlift",
            "jmeasure": "JMeasure",
            "implicative": "ImplicStat",
            "implicativep": "p.value",
            "implicativepb1": "p.valueB1",
            "implicativepb2": "p.valueB2",
        }
        key = sort_map.get(sortv, sortv)
        # R's order() (seqerules.R:77-78): text in R's collation, not code points, so the
        # rule "(a) => (B)" sorts before "(B) => (a)" (TX-AUD-20261019-083); a pandas string
        # dtype (pandas 3's default) is text too (batch-11 review M-1).
        is_text = pd.api.types.is_object_dtype(out[key]) or pd.api.types.is_string_dtype(out[key])
        text_rank = (
            (lambda column: column.map(r_character_rank(column.tolist()))) if is_text else None
        )
        out = out.sort_values(
            key, ascending=not decreasing, kind="stable", key=text_rank
        ).reset_index(drop=True)
    return out


class SubSeqElistChisq(SubSeqElist):
    """Discriminating subsequences with chi-square statistics (R ``subseqelistchisq``).

    The result of :func:`seqecmpgroup`: a :class:`SubSeqElist` whose ``data`` carries
    the per-group ``Freq``/``Resid`` plus ``p.value`` / ``statistic`` / ``index``,
    ordered by statistic. ``labels`` are the group levels and ``bonferroni`` records
    whether the Bonferroni correction was applied and over how many tests.
    """

    def __init__(
        self, eseq, constraint, subseq_objs, data, type="chisq", labels=None, bonferroni=None
    ):
        super().__init__(eseq, constraint, subseq_objs, data, type=type)
        self.labels = list(labels) if labels is not None else []
        self.bonferroni = dict(bonferroni) if bonferroni is not None else {}

    def __repr__(self) -> str:
        return f"SubSeqElistChisq(n={len(self.subseq)}, groups={self.labels})"


def _group_presence_table(weights, group_codes, presence, nlevels: int) -> np.ndarray:
    """The ``nlevels x 2`` weighted group-by-presence table (R ``xtabs(ww ~ group + p)``).

    ``xtabs`` fills each cell with ``tapply(w, by, sum)``: the cell's weights
    in input order, summed as R's ``sum`` does (:func:`_r_sum`).
    """
    table = np.zeros((nlevels, 2), dtype=float)
    for gi in range(nlevels):
        mask = group_codes == gi
        table[gi, 0] = _r_sum(weights[mask & (presence == 0)])
        table[gi, 1] = _r_sum(weights[mask & (presence == 1)])
    return table


def _chisq_test(table: np.ndarray):
    """R ``chisq.test`` on a (possibly weighted) contingency table, bit for bit.

    Returns ``(statistic, p_value, residuals, expected, row_sums)``. The
    arithmetic is that of ``stats::chisq.test`` (R 4.6.1, no simulation):
    ``sr <- rowSums(x)``, ``sc <- colSums(x)``, ``n <- sum(x)``,
    ``E <- outer(sr, sc)/n``; for a 2 x 2 table one scalar Yates correction
    ``YATES <- min(0.5, abs(x - E))`` (0 otherwise); then
    ``STATISTIC <- sum((abs(x - E) - YATES)^2/E)`` and residuals
    ``(x - E)/sqrt(E)``. The sums follow R's order and long-double precision
    (:func:`_r_sum`), so statistics that tie in exact arithmetic keep R's own
    last-bit differences -- seqecmpgroup orders subsequences on them
    (2026-09-22 audit S5-16). The p-value comes from SciPy's ``chi2.sf``,
    which agrees with R's ``pchisq`` to rounding.
    """
    from scipy.stats import chi2 as _chi2

    table = np.asarray(table, dtype=float)
    nr, nc = table.shape
    n = _r_sum(table)
    sr = np.array([_r_sum(table[row, :]) for row in range(nr)], dtype=float)
    sc = np.array([_r_sum(table[:, column]) for column in range(nc)], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        E = np.outer(sr, sc) / n
        deviation = np.abs(table - E)
        # One correction for the whole table; np.min, unlike Python's min,
        # propagates a NaN as R's min does.
        yates = float(np.min(np.append(deviation, 0.5))) if nr == 2 and nc == 2 else 0.0
        stat = _r_sum((deviation - yates) ** 2 / E)
        resid = (table - E) / np.sqrt(E)
    df = (nr - 1) * (nc - 1)
    pval = float(_chi2.sf(stat, df))
    return stat, pval, resid, E, sr


def _is_missing(value: object) -> bool:
    missing = pd.isna(value)
    return isinstance(missing, (bool, np.bool_)) and bool(missing)


def _group_factor(group) -> tuple[List[str], np.ndarray]:
    """R's ``factor(group)``: level names and a 0-based code per value (-1 = missing).

    A pandas categorical keeps its category order with unused categories
    dropped; other values take R's factor order
    (:func:`~traminex._categorical_order.stable_sorted_unique`: numbers by
    value, text in R's collation, which needs ICU). Levels are named as
    R's ``as.character`` prints them (``2``, ``1e+05``, ``TRUE``), and missing
    values (``None``, ``NaN``, ``pd.NA``) get no level, as R's ``NA`` does.
    """
    if isinstance(group, pd.Categorical) or isinstance(
        getattr(group, "dtype", None), pd.CategoricalDtype
    ):
        factor = pd.Categorical(group).remove_unused_categories()
        levels = list(factor.categories)
        codes = np.asarray(factor.codes, dtype=int)
    else:
        values = list(group)
        levels = stable_sorted_unique(values)
        position = {value: code for code, value in enumerate(levels)}
        codes = np.asarray(
            [-1 if _is_missing(value) else position[value] for value in values], dtype=int
        )
    return [str(r_as_character(level)) for level in levels], codes


def seqecmpgroup(
    subseq: SubSeqElist,
    group,
    method: str = "chisq",
    pvalue_limit: Optional[float] = None,
    weighted: bool = True,
) -> SubSeqElistChisq:
    """Most discriminating subsequences between groups (R ``seqecmpgroup``).

    Faithful port of ``TraMineR::seqecmpgroup`` (``seqecmpgroup.R``): for each mined
    subsequence, cross-tabulate (weighted) ``group x presence`` and run a chi-square
    test; report the per-group frequency of the subsequence (``Freq.<level>``), the
    Pearson residuals (``Resid.<level>``), the chi-square ``statistic`` and
    ``p.value`` (Bonferroni-corrected when ``method="bonferroni"``), reordered by
    statistic descending.

    Parameters
    ----------
    subseq : SubSeqElist
        A :func:`seqefsub` result.
    group : sequence
        Group label per event sequence, converted as R's ``factor(group)``: a
        pandas categorical keeps its category order (unused categories are
        dropped), numbers sort as numbers and text in the package's factor
        order, and the levels are named as R prints the values (``2``,
        ``1e+05``, ``TRUE``). A missing value (``None``, ``NaN``, ``pd.NA``)
        leaves its sequence out of every group-by-presence table, as R leaves
        out ``NA``, while it still counts in the presence total.
    method : {"chisq", "bonferroni"}, default "chisq"
        ``"bonferroni"`` applies ``1 - (1 - p)^ntest`` to each p-value.
    pvalue_limit : float, optional
        Keep the (by-statistic top) subsequences; the number kept is the count of
        subsequences whose p-value is ``<= pvalue_limit`` (default 2, i.e. all).
        When none passes, the result is empty. R keeps ``1:0`` rows there
        (``seqecmpgroup.R:88``), which still returns the top subsequence
        although its p-value fails the limit; traminex deliberately does not.
    weighted : bool, default True
        Use the event-sequence weights in the contingency table.

    Returns
    -------
    SubSeqElistChisq
    """
    if not is_subseqelist(subseq):
        raise TypeError(" [!] subseq should be a subseqelist")
    bonferroni = method == "bonferroni"
    method_eff = "chisq" if bonferroni else method
    if method_eff != "chisq":
        raise ValueError(" [!] This method is not (yet) implemented")
    if pvalue_limit is None:
        pvalue_limit = 2

    eseqs = subseq.eseq.sequences
    if weighted:
        ww = np.array([float(s.weight) for s in eseqs], dtype=float)
    else:
        ww = np.ones(len(eseqs), dtype=float)

    glevels, group_codes = _group_factor(group)
    if len(group_codes) != len(eseqs):
        # R stops in xtabs with "variable lengths differ (found for 'group')";
        # without this the mismatch surfaced as a numpy broadcast error raised
        # from inside the private counting loop.
        raise ValueError(
            f" [!] group has {len(group_codes)} values but subseq covers "
            f"{len(eseqs)} event sequences"
        )
    if len(glevels) < 2:
        # A one-level grouping gives a 1 x 2 table, so the chi-square has zero
        # degrees of freedom and every p-value is NaN. `NaN <= pvalue_limit` is
        # False, which silently dropped every subsequence from the result. R
        # aborts with "incorrect number of dimensions"; name the argument
        # instead of returning an empty frame that looks like a finding.
        raise ValueError(
            " [!] seqecmpgroup compares groups: group must have at least two "
            f"distinct levels, found {len(glevels)}"
        )
    seqmatrix = seqeapplysub(subseq, method="presence").to_numpy()
    # R's length(group): sequences with a missing group count here and in the
    # presence total, although xtabs leaves them out of the table.
    n = len(group_codes)
    ntest = len(subseq.subseq)

    rows: List[dict] = []
    for i in range(seqmatrix.shape[1]):
        col = seqmatrix[:, i]
        sp = float(np.sum(col))
        row: dict = {}
        if sp > 0 and sp < n:
            table = _group_presence_table(ww, group_codes, col, len(glevels))
            stat, pval, resid, _E, sr = _chisq_test(table)
            if bonferroni:
                pval = 1.0 - (1.0 - pval) ** ntest
            with np.errstate(divide="ignore", invalid="ignore"):
                freq = table[:, 1] / sr
            row["p.value"] = float(pval)
            row["statistic"] = float(stat)
            row["index"] = i + 1
            for gi, gl in enumerate(glevels):
                row[f"Freq.{gl}"] = float(freq[gi])
            for gi, gl in enumerate(glevels):
                row[f"Resid.{gl}"] = float(resid[gi, 1])
        else:
            row["p.value"] = 1.0
            row["statistic"] = 0.0
            row["index"] = i + 1
            for gl in glevels:
                row[f"Freq.{gl}"] = sp / n
            for gl in glevels:
                row[f"Resid.{gl}"] = 0.0
        rows.append(row)

    res = pd.DataFrame(rows)
    pvalues = res["p.value"].to_numpy(dtype=float)
    undefined = np.flatnonzero(np.isnan(pvalues))
    if undefined.size:
        # A group level -- or a presence class -- whose sequences all carry zero
        # weight makes an expected count zero, so the statistic is NaN. Comparing
        # NaN against pvalue_limit is False, which used to drop those rows
        # without a word. R aborts with "NA/NaN argument"; say what happened.
        raise ValueError(
            " [!] the group x presence table has zero total weight in a row or "
            "column, so the chi-square statistic is undefined for subsequence(s) "
            + ", ".join(repr(subseq.subseq[int(i)]) for i in undefined[:5])
            + (" ..." if undefined.size > 5 else "")
        )
    keep = int(np.sum(pvalues <= pvalue_limit))
    order = np.argsort(-res["statistic"].to_numpy(), kind="stable")[:keep]

    support_ord = subseq.data["Support"].to_numpy()[order]
    res_ord = res.iloc[order].reset_index(drop=True)
    data = pd.concat([pd.DataFrame({"Support": support_ord}), res_ord], axis=1)
    sub_objs_ord = [subseq._subseq_objs[j] for j in order]
    return SubSeqElistChisq(
        subseq.eseq,
        subseq.constraint,
        sub_objs_ord,
        data,
        type=method,
        labels=glevels,
        bonferroni={"used": bonferroni, "ntest": ntest},
    )


def _seqedist_option_code(
    value: str | int | bool,
    *,
    name: str,
    public_choices: tuple[str, ...],
    public_codes: tuple[int, ...],
    numeric_codes: tuple[int, ...],
) -> int:
    """Translate a public string or retained Python scalar code."""
    if isinstance(value, str):
        normalized = value.lower()
        if not normalized:
            raise ValueError(f"{name} must not be empty")
        exact = [i for i, choice in enumerate(public_choices) if choice == normalized]
        matches = exact or [
            i for i, choice in enumerate(public_choices) if choice.startswith(normalized)
        ]
        if len(matches) != 1:
            detail = "ambiguous" if matches else "unknown"
            choices = ", ".join(repr(choice) for choice in public_choices)
            raise ValueError(f"{name} value {value!r} is {detail}; expected one of {choices}")
        return public_codes[matches[0]]

    if isinstance(value, (bool, np.bool_)):
        return int(value)
    if isinstance(value, (int, np.integer)):
        code = int(value)
        if code not in numeric_codes:
            raise ValueError(f"{name} numeric code must be one of {numeric_codes}")
        return code
    raise TypeError(f"{name} must be a scalar string, integer code, or boolean")


def seqedist(
    eseq: SeqElist,
    idcost,
    vparam: float,
    interval: str | int | bool = "No",
    norm: str | int | bool = "YujianBo",
) -> np.ndarray:
    """Pairwise event-sequence edit distance (``TraMineRextras::seqedist``).

    The public contract follows the TraMineRextras wrapper, including its
    case-insensitive unique partial matching and string defaults. The translated
    codes are passed to TraMineR's ``C_tmrseqedist`` engine, a dynamic-programming
    edit distance where each event has an indel cost and ``vparam`` penalises age
    differences between aligned transitions.

    Parameters
    ----------
    eseq : SeqElist
        Event sequences from :func:`~traminex.events.seqecreate`.
    idcost : array-like
        Per-event indel cost, length = size of the event dictionary (event code
        ``k`` uses ``idcost[k - 1]``).
    vparam : float
        Time-difference weight ``V``.
    interval : {"No", "previous", "next"} or int or bool, default "No"
        Public strings map respectively to engine codes ``0``, ``2``, and ``3``.
        Case-insensitive unique prefixes are accepted. Integer codes ``0`` through
        ``3`` and booleans (``False -> 0``, ``True -> 1``) are retained Python
        extensions.
    norm : {"none", "max", "YujianBo"} or int or bool, default "YujianBo"
        Public strings map respectively to engine codes ``0``, ``1``, and ``2``.
        Case-insensitive unique prefixes are accepted. Integer codes ``0`` through
        ``2`` and booleans (``False -> 0``, ``True -> 1``) are retained Python
        extensions.

    Returns
    -------
    numpy.ndarray
        The symmetric ``(n_seq x n_seq)`` distance matrix (zero diagonal).

    Notes
    -----
    Two departures correct defects of TraMineR's engine (semantic ledger,
    ``seqedist`` row): an event sequence with no events is at the cost of
    inserting every event of the other sequence, whatever the row order (R gives
    0 when the empty sequence comes first); and a raw distance of 0 normalises
    to 0 under ``"YujianBo"`` (R returns NaN when the bound is 0 too, e.g. two
    empty sequences or every ``idcost`` 0).
    """
    if not is_seqelist(eseq):
        raise TypeError(" [!] eseq should be a seqelist. See help on seqecreate.")
    idcost_array = np.asarray(idcost, dtype=float)
    if idcost_array.ndim != 1 or idcost_array.size != len(eseq.dictionary):
        raise ValueError("You should specify one indel cost per event in the event alphabet")

    interval_code = _seqedist_option_code(
        interval,
        name="interval",
        public_choices=("no", "previous", "next"),
        public_codes=(0, 2, 3),
        numeric_codes=(0, 1, 2, 3),
    )
    norm_code = _seqedist_option_code(
        norm,
        name="norm",
        public_choices=("none", "max", "yujianbo"),
        public_codes=(0, 1, 2),
        numeric_codes=(0, 1, 2),
    )
    eng = [_engine.eseq_to_engine(s) for s in eseq.sequences]
    return _engine.tmrseqedist(eng, idcost_array, float(vparam), norm_code, interval_code)


def _subseq_to_event_string(subseq) -> str:
    """Render a subsequence to a TraMineR event-label string.

    A generic subsequence label string (``"(A)-(B)"``, possibly with comma-grouped
    simultaneous events) is returned unchanged; a tuple of single event labels is
    rendered as ``"(A)-(B)"``.
    """
    if isinstance(subseq, str):
        return subseq
    return "-".join(f"({lab})" for lab in subseq)


def _ntrans_nevent(subseq):
    """The ``seqentrans`` ntrans/nevent formula (seqentrans.R lines 7-10)."""
    s = _subseq_to_event_string(subseq)
    ntrans = len(s.split("-"))
    nevent = ntrans - 1 + len(s.split(","))
    return ntrans, nevent


def seqentrans(fsubseq, *, avg_occ: bool = False):
    """Augment a frequent-subsequence list with ntrans / nevent / Avg.occ (R ``seqentrans``).

    Faithful port of ``TraMineRextras::seqentrans`` (``seqentrans.R``): adds, per
    subsequence, ``ntrans`` (number of ``-``-separated transitions of the label),
    ``nevent`` (``ntrans - 1 + `` number of ``,``-separated tokens) and -- when
    ``avg_occ`` -- ``Avg.occ``. Performs no mining.

    ``Avg.occ`` is the average number of occurrences by sequence: ``Count`` over the total
    it was taken on (the total weight after weighted mining, the number of sequences after
    ``weighted=False``), where R always divides by the total weight
    (TX-AUD-20261019-016).

    Accepts the faithful :class:`SubSeqElist` from :func:`seqefsub` (a copy with the
    columns added to its ``data`` is returned and the argument is left unchanged, as R
    returns a modified copy of ``fsubseq``),
    or -- for convenience -- a plain DataFrame carrying a ``subseq`` column (and, for
    ``avg_occ``, a ``count`` column plus a ``seqe_weight_total`` attr, which the caller
    sets to the total the counts were taken on).

    Returns
    -------
    SubSeqElist or pandas.DataFrame
        An augmented copy of the input (a :class:`SubSeqElist` for a SubSeqElist input).
    """
    if is_subseqelist(fsubseq):
        # R returns a modified copy and leaves its argument unchanged.
        result = copy.copy(fsubseq)
        result.data = fsubseq.data.copy()
        nt_ne = [_ntrans_nevent(s) for s in result.subseq]
        result.data["ntrans"] = [nt for nt, _ in nt_ne]
        result.data["nevent"] = [ne for _, ne in nt_ne]
        if avg_occ:
            # The average by sequence divides Count by the total it was taken on; R
            # divides by the total weight even after weighted=False mining
            # (seqentrans.R:11-14; TX-AUD-20261019-016).
            total = getattr(fsubseq, "count_total", None)
            if total is None:
                total = float(sum(s.weight for s in fsubseq.eseq))
            if float(total) == 0.0:
                raise ValueError("seqentrans(avg_occ=True) needs a positive total weight")
            result.data["Avg.occ"] = result.data["Count"].astype(float) / float(total)
        return result

    if not isinstance(fsubseq, pd.DataFrame):
        raise TypeError(
            "fsubseq must be a SubSeqElist (seqefsub) or a DataFrame with a subseq column"
        )
    if "subseq" not in fsubseq.columns:
        raise ValueError("seqentrans expects a SubSeqElist or a DataFrame with a 'subseq' column")
    out = fsubseq.copy()
    out.attrs = dict(fsubseq.attrs)
    nt_ne = [_ntrans_nevent(s) for s in out["subseq"]]
    out["ntrans"] = [nt for nt, _ in nt_ne]
    out["nevent"] = [ne for _, ne in nt_ne]
    if avg_occ:
        if "count" not in out.columns:
            raise ValueError("seqentrans(avg_occ=True) requires a 'count' column")
        wtot = fsubseq.attrs.get("seqe_weight_total")
        if wtot is None or float(wtot) == 0.0:
            raise ValueError(
                "seqentrans(avg_occ=True) needs a positive 'seqe_weight_total' attr, set by "
                "the caller to the total the counts were taken on (avg_occ = count / total)"
            )
        out["avg_occ"] = out["count"].astype(float) / float(wtot)
    return out
