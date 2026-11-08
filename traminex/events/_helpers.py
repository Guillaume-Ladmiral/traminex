"""Non-R event convenience helpers (no TraMineR equivalent).

These traminex-only helpers (``events_transition_matrix`` /
``events_transition_degrees``, ``frequent_events`` / ``frequent_ngrams`` /
``frequent_event_pairs`` / ``frequent_event_pairs_timegap``,
``events_to_state_sequences`` and the flat :class:`EventSequence` container)
provide lightweight exploratory tools that have no R counterpart. The faithful,
R-parity event-mining surface (``seqeconstraint`` / ``seqefsub`` / ``seqeapplysub``
/ ``seqerules`` / ``seqecmpgroup`` / ``seqedist`` / ``seqentrans``) lives in
:mod:`traminex.events._mining` on the ported C++ engine
(:mod:`traminex.events._engine`).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from numbers import Real
from typing import Any, Dict, List, NoReturn, Optional, SupportsIndex, Tuple

import numpy as np
import pandas as pd

from .._categorical_order import r_seqstatl_order
from ..sequences import STSList, _require_stslist


class _ReadOnlyList(list):
    """A stored event vector: reads like a ``list``, refuses in-place writes.

    ``EventSequence`` validates its (id, time, label) triple on assignment and
    ``SeqElist``'s flat views must describe the same events as its sequences;
    an in-place write (``append``, item assignment, ...) would slip past both,
    as ``e.time[1] = nan`` did before its array was frozen (cycle-10 N-19;
    2026-09-22 audit S5-05). ``+=`` builds a new list, so on an
    ``EventSequence`` it goes through the validated assignment.
    """

    __slots__ = ()

    @staticmethod
    def _refuse() -> NoReturn:
        raise TypeError("this stored event vector is read-only; assign a new vector to change it")

    def __setitem__(self, *args: Any) -> None:
        self._refuse()

    def __delitem__(self, *args: Any) -> None:
        self._refuse()

    def append(self, *args: Any) -> None:
        self._refuse()

    def extend(self, *args: Any) -> None:
        self._refuse()

    def insert(self, *args: Any) -> None:
        self._refuse()

    def pop(self, *args: Any) -> Any:
        self._refuse()

    def remove(self, *args: Any) -> None:
        self._refuse()

    def clear(self) -> None:
        self._refuse()

    def sort(self, *args: Any, **kwargs: Any) -> None:
        self._refuse()

    def reverse(self) -> None:
        self._refuse()

    def __iadd__(self, other: Iterable[Any]) -> List[Any]:  # type: ignore[override,misc]
        return list(self) + list(other)

    def __imul__(self, count: SupportsIndex) -> List[Any]:  # type: ignore[override]
        return list(self) * count

    def __reduce__(self) -> Tuple[type, Tuple[List[Any]]]:
        return (type(self), (list(self),))


def _events_by_id_preserve_ties(
    events: "EventSequence",
) -> Dict[str, List[Tuple[float, str]]]:
    """Group events per id, ordering by time and preserving original order on ties."""
    by_id: Dict[str, List[Tuple[float, int, str]]] = {}
    for pos, (pid, t, lab) in enumerate(zip(events.id, events.time, events.label)):
        by_id.setdefault(pid, []).append((float(t), pos, lab))
    grouped: Dict[str, List[Tuple[float, str]]] = {}
    for pid, items in by_id.items():
        items.sort(key=lambda item: (item[0], item[1]))
        grouped[pid] = [(t, lab) for t, _, lab in items]
    return grouped


def frequent_events(events: "EventSequence", min_support: int = 1) -> List[Tuple[str, int, float]]:
    """Count event labels across all ids.

    Returns a list of (label, count, support_fraction) sorted by count desc.
    ``support_fraction``'s denominator is the total number of event rows,
    which differs from the pair-based denominators of
    :func:`frequent_event_pairs` and :func:`frequent_event_pairs_timegap`.
    """
    _require_event_sequence(events)
    counts: Dict[str, int] = {}
    for lab in events.label:
        counts[lab] = counts.get(lab, 0) + 1
    total = len(events.label)
    out: List[Tuple[str, int, float]] = []
    for k, v in counts.items():
        if v >= min_support:
            out.append((k, int(v), (float(v) / float(total)) if total > 0 else 0.0))
    out.sort(key=lambda z: (-z[1], z[0]))
    return out


def events_transition_matrix(
    events: "EventSequence",
    normalize: bool = False,
    id_weights: Mapping[str, float] | None = None,
) -> Tuple[List[str], np.ndarray]:
    """Build adjacency matrix of event label transitions (ordered pairs within id).

    Returns (labels, M) where M[i,j] counts transitions from labels[i] -> labels[j].
    If normalize=True, row-normalizes to probabilities.
    If id_weights provided, transitions from an id are weighted accordingly.
    An id without an entry defaults to weight 1 and unmatched entries are
    inert (reusable mappings) -- mind typos, which silently default.
    """
    _require_event_sequence(events)
    weights = _validate_id_weights(id_weights)
    labs = sorted(set(events.label))
    idx = {lab: i for i, lab in enumerate(labs)}
    K = len(labs)
    M = np.zeros((K, K), dtype=float)
    by_id = _events_by_id_preserve_ties(events)
    for pid, items in by_id.items():
        w = 1.0 if weights is None else weights.get(pid, 1.0)
        labs_seq = [lab for _, lab in items]
        for a, b in zip(labs_seq[:-1], labs_seq[1:]):
            M[idx[a], idx[b]] += w
    if normalize:
        row = M.sum(axis=1, keepdims=True)
        M = np.divide(M, row, out=np.zeros_like(M), where=row > 0)
    return labs, M


def _require_event_sequence(events: object) -> None:
    """Refuse input without flat event views before an attribute is read (TX-REP-20261018-001).

    The helpers read the flat ``id``, ``time`` and ``label`` views that an
    :class:`EventSequence` and a ``SeqElist`` (from ``seqecreate``) both expose.
    """
    if not all(hasattr(events, view) for view in ("id", "time", "label")):
        raise TypeError(
            "events must be an EventSequence or a SeqElist (flat id, time and label views)"
        )


@dataclass
class EventSequence:
    id: List[str]
    time: np.ndarray  # per-event times
    label: List[str]  # per-event labels

    def __post_init__(self) -> None:
        """Validate and defensively normalize the three parallel event vectors."""
        normalized = _validate_event_vectors(self.id, self.time, self.label)
        for name, value in zip(("id", "time", "label"), normalized):
            object.__setattr__(self, name, value)
        object.__setattr__(self, "_validated", True)

    def __setattr__(self, name: str, value: object) -> None:
        """Re-validate the triple on every post-construction field assignment.

        Without this, mutating one vector after construction silently restores
        the truncated-``zip`` behavior the constructor validation prevents.
        """
        if name in ("id", "time", "label") and self.__dict__.get("_validated"):
            candidate = {"id": self.id, "time": self.time, "label": self.label}
            candidate[name] = value
            normalized = _validate_event_vectors(
                candidate["id"], candidate["time"], candidate["label"]
            )
            for field_name, field_value in zip(("id", "time", "label"), normalized):
                object.__setattr__(self, field_name, field_value)
            return
        object.__setattr__(self, name, value)


def _validate_event_vectors(
    id_values: object,
    time_values: object,
    label_values: object,
) -> tuple[List[str], np.ndarray, List[str]]:
    """Validate and normalize the (id, time, label) triple as one unit."""
    normalized_id = _validate_string_vector(id_values, name="id")
    normalized_label = _validate_string_vector(label_values, name="label")

    raw_time = np.asarray(time_values)
    if raw_time.ndim != 1:
        raise ValueError("time must be one-dimensional")
    if any(
        isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) for value in raw_time
    ):
        raise TypeError("time must contain numeric values")
    normalized_time = np.asarray([float(value) for value in raw_time], dtype=float)
    if not np.all(np.isfinite(normalized_time)):
        raise ValueError("time must contain only finite values")
    # ``EventSequence.__setattr__`` re-validates on assignment, but an in-place
    # write went straight past it: ``e.time[1] = nan`` left the object holding a
    # value this very check rejects (cycle-10 N-19). This array is already a
    # fresh copy built from the comprehension above, so freezing it constrains
    # only the stored vector, never the caller's argument. Reassignment stays
    # the supported way to change it.
    normalized_time.flags.writeable = False

    if not (len(normalized_id) == len(normalized_time) == len(normalized_label)):
        raise ValueError("id, time, and label must have the same length")
    return normalized_id, normalized_time, normalized_label


def _validate_string_vector(value: object, *, name: str) -> List[str]:
    """Return a defensive, read-only list copy of a non-string iterable of strings."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise TypeError(f"{name} must be a collection of strings")
    try:
        normalized = list(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be a collection of strings") from exc
    if any(not isinstance(item, str) for item in normalized):
        raise TypeError(f"{name} must be a collection of strings")
    # Frozen for the same reason as ``time`` (cycle-10 N-19): an in-place
    # ``append`` or item write would bypass the triple validation (S5-05).
    return _ReadOnlyList(normalized)


def _validate_id_weights(
    id_weights: Mapping[str, float] | None,
) -> Optional[Dict[str, float]]:
    """Validate every supplied id weight before transition calculation.

    Keys need not be observed ids: an id without an entry defaults to weight
    1 and an unmatched entry is inert, so one weights mapping can serve
    several datasets (pinned by the event-sequence contract tests). Mind
    typos: a misspelled id silently leaves its transitions at weight 1.
    """
    if id_weights is None:
        return None
    if not isinstance(id_weights, Mapping):
        raise TypeError("id_weights must be a mapping")

    normalized: Dict[str, float] = {}
    for pid, value in id_weights.items():
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
            raise TypeError(f"id_weights[{pid!r}] must be a real number")
        numeric = float(value)
        if not np.isfinite(numeric):
            raise ValueError(f"id_weights[{pid!r}] must be finite")
        if numeric < 0:
            raise ValueError(f"id_weights[{pid!r}] must be non-negative")
        normalized[pid] = numeric
    return normalized


def events_to_state_sequences(
    events: EventSequence,
    *,
    states: Optional[List[str]] = None,
    time_grid: Optional[np.ndarray] = None,
    fill: str = "_",
    align_missing: bool = False,
) -> STSList:
    """Convert event sequences to regular state sequences (simple spells-to-grid).

    - If time_grid is None, builds an integer-step grid from ``floor(min time)``
      to ``ceil(max time)``, so a state entered at a fractional time is still
      represented at the following grid point.
    - Each subject is in the last seen event label until the next event; starts with fill.
    - Rows are ordered by sorted id (the pinned contract); note the counting
      helpers group by first appearance instead.
    - Without ``states``, the alphabet is R's ``seqdef`` order of the grid's states
      (numbers by value, R's collation otherwise; owner decision of 2026-10-19).
    """
    _require_event_sequence(events)
    ids = sorted(set(events.id))
    # Build per-id event lists
    per = _events_by_id_preserve_ties(events)
    # Determine global time grid
    if time_grid is None:
        if len(events.time) == 0:
            grid = np.empty(0, dtype=float)
        else:
            tmin = math.floor(min(events.time))
            tmax = math.ceil(max(events.time))
            grid = np.arange(tmin, tmax + 1)
    else:
        grid = np.asarray(time_grid)
    L = len(grid)
    rows: list[list[str]] = []
    for pid in ids:
        ev = per.get(pid, [])
        cur = fill
        row = []
        k = 0
        for t in grid:
            while k < len(ev) and ev[k][0] <= t:
                cur = ev[k][1]
                k += 1
            row.append(cur)
        rows.append(row)
    matrix = np.asarray(rows, dtype=object).reshape(len(rows), L)
    # Create STSList via seqdef-like process
    observed_states = [str(s) for s in pd.unique(matrix.ravel(order="C").astype(str)).tolist()]
    if states is None:
        # R's seqdef order of the grid's states (owner decision D2 of 2026-10-19).
        all_states = r_seqstatl_order(observed_states)
    else:
        all_states = list(states)
        if fill in observed_states and fill not in all_states:
            all_states.append(fill)
        unknown = [state for state in observed_states if state not in all_states]
        if unknown:
            raise ValueError(f"states does not include observed states: {unknown}")
    # If aligning missing, ensure the fill symbol is part of the alphabet
    if align_missing and fill not in all_states:
        all_states = list(all_states) + [fill]
    code = {s: i for i, s in enumerate(all_states)}
    if matrix.size == 0:
        mat = np.empty(matrix.shape, dtype=int)
    else:
        mat = np.vectorize(code.__getitem__, otypes=[int])(matrix.astype(str))
    if align_missing:
        nr_code = code.get(fill, None)
        return STSList(
            data=mat,
            alphabet=list(all_states),
            nr=fill,
            nr_code=nr_code,
            ids=list(ids),
        )
    return STSList(data=mat, alphabet=list(all_states), ids=list(ids))


@_require_stslist
def frequent_ngrams(
    x: STSList, n: int = 2, min_support: int = 1
) -> List[Tuple[Tuple[str, ...], int, float]]:
    """Count contiguous state n-grams across all sequences.

    **traminex-only extension -- no TraMineR counterpart.** This docstring is
    the contract; it is pinned by ``tests/test_patterns_contract.py`` with
    hand-computed expectations. It consumes an :class:`STSList` (state
    sequences), not an :class:`EventSequence`; it stays exported from
    ``traminex.events`` for compatibility with the frozen 226-name API.

    Contract
    --------
    Every window ``row[t:t+n]`` of every sequence contributes one occurrence.
    ``support_fraction`` is ``count / (nseq * (seqlen - n + 1))``, the total
    number of window positions. An n-gram is returned when its count is
    ``>= min_support`` (an inclusive floor); results are sorted by count
    descending, ties broken by the n-gram ascending.

    ``n <= 0`` or ``n > seqlen`` returns an empty list rather than raising.

    As in :func:`traminex.patterns.frequent_subsequences`, the void padding of
    a ragged object counts as an ordinary state and the denominator uses the
    full matrix width.
    """
    if n <= 0 or n > x.seqlen:
        return []
    counts: dict[Tuple[int, ...], int] = {}
    for row in x.data:
        for t in range(0, x.seqlen - n + 1):
            gram = tuple(int(v) for v in row[t : t + n])
            counts[gram] = counts.get(gram, 0) + 1
    total_positions = x.nseq * (x.seqlen - n + 1)
    items: List[Tuple[Tuple[str, ...], int, float]] = []
    for k, v in counts.items():
        if v >= min_support:
            items.append(
                (tuple(x._state_for_code(i) for i in k), int(v), float(v) / float(total_positions))
            )
    items.sort(key=lambda z: (-z[1], z[0]))
    return items


def frequent_event_pairs(
    events: EventSequence, min_support: int = 1
) -> List[Tuple[Tuple[str, str], int, float]]:
    """Count ordered adjacent event label pairs per id.

    For each id, sort events by time and count consecutive label pairs (lab_t, lab_{t+1}).
    Returns a list of (pair, count, support_fraction) sorted by count desc.
    ``support_fraction``'s denominator is the total number of consecutive
    pairs, which differs from :func:`frequent_events` (event rows) and
    :func:`frequent_event_pairs_timegap` (gap-bounded ordered pairs).
    """
    _require_event_sequence(events)
    # build per-id sorted labels
    by_id = _events_by_id_preserve_ties(events)
    total_pairs = 0
    counts: Dict[Tuple[str, str], int] = {}
    for pid, items in by_id.items():
        labs = [lab for _, lab in items]
        for a, b in zip(labs[:-1], labs[1:]):
            counts[(a, b)] = counts.get((a, b), 0) + 1
            total_pairs += 1
    out: List[Tuple[Tuple[str, str], int, float]] = []
    for k, v in counts.items():
        if v >= min_support:
            out.append((k, int(v), (float(v) / float(total_pairs)) if total_pairs > 0 else 0.0))
    out.sort(key=lambda z: (-z[1], z[0]))
    return out


def events_transition_degrees(
    events: EventSequence,
    normalize: bool = False,
    id_weights: Mapping[str, float] | None = None,
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """Compute in/out degrees from event transition matrix.

    Returns (labels, indegree, outdegree). If normalize=True, divides by total transitions.
    If id_weights provided, weights apply to transitions from an id.
    """
    _require_event_sequence(events)
    labs, M = events_transition_matrix(events, normalize=False, id_weights=id_weights)
    indeg = M.sum(axis=0)
    outdeg = M.sum(axis=1)
    if normalize and M.sum() > 0:
        indeg = indeg / M.sum()
        outdeg = outdeg / M.sum()
    return labs, indeg, outdeg


def frequent_event_pairs_timegap(
    events: EventSequence, max_gap: int, min_support: int = 1
) -> List[Tuple[Tuple[str, str], int, float]]:
    """Count ordered event pairs per id with time gap <= max_gap.

    Returns (pair, count, support_fraction) sorted by count desc.
    ``support_fraction``'s denominator is the total number of ordered pairs
    within ``max_gap`` -- it grows with ``max_gap`` and is not comparable to
    the denominators of :func:`frequent_events` (event rows) or
    :func:`frequent_event_pairs` (consecutive pairs only).
    """
    _require_event_sequence(events)
    by_id = _events_by_id_preserve_ties(events)
    total = 0
    counts: Dict[Tuple[str, str], int] = {}
    for pid, items in by_id.items():
        for i in range(len(items) - 1):
            t1, a = items[i]
            for j in range(i + 1, len(items)):
                t2, b = items[j]
                if t2 - t1 <= max_gap:
                    counts[(a, b)] = counts.get((a, b), 0) + 1
                    total += 1
                else:
                    break
    out: List[Tuple[Tuple[str, str], int, float]] = []
    for k, v in counts.items():
        if v >= min_support:
            out.append((k, int(v), (float(v) / float(total)) if total > 0 else 0.0))
    out.sort(key=lambda z: (-z[1], z[0]))
    return out
