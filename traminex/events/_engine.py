"""Faithful pure-Python port of TraMineR's vendored C++ event-mining engine (RP5b).

This module is **internal** -- it ports the algorithms in
``SOURCES/TraMineR/src/`` that the R event functions ``.Call`` into, with the
corrections of R that the parity ledger records (its ``seqedist`` and
``seqefsub`` rows):

- ``eventseq.cpp`` -- ``Sequence`` / ``SequenceEventNode``: the event-node chain,
  ``addEvent`` insertion, ``sprint`` label rendering, ``contain``, and the **six
  count methods** (``count`` = COBJ/CDIST_O inner matcher, ``count3`` = CWIN,
  ``count4`` = CMINWIN, ``count5`` = CDIST, ``first_occurence`` = age);
- ``constraint.cpp`` -- the ``-1`` "no bound" defaults, which R sets to
  ``+/-DBL_MAX`` and this port to ``+/-inf`` (as CMINWIN's running minimum
  starts at ``inf`` where R's starts at ``DBL_MAX``). The bounds are only
  compared with, or added to, finite gaps and ages, so both give the same
  results for event times below ``2**970`` (about 1e292), above which
  ``DBL_MAX`` plus a time overflows; CWIN and CMINWIN with default bounds were
  checked against R (TX-AUD-20261010-078);
- ``prefixtree.cpp`` / ``treeeventnode.cpp`` / ``treeeventmap.cpp`` -- the
  prefix-tree frequent-subsequence miner;
- ``tmrsequence.cpp`` -- ``tmrfindsubsequences`` (the miner driver) and
  ``tmrmatrixsubseqinseq`` (the subsequence-count matrix);
- ``seqedist.cpp`` -- ``tmrseqedist`` (the event-sequence edit distance).

The results are golden-pinned via the public wrappers in
:mod:`traminex.events._mining` against live TraMineR 2.2-14. A pybind11/C++
acceleration, if ever wanted, is an RP9 ``atol=0`` byte-parity concern -- never a
faithfulness change. There is no Sequenzo borrow here (its event core is
non-faithful); this is a BUILD-from-(GPL TraMineR's own C/R) port.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from .._rformat import r_number_text

INF = float("inf")


# --------------------------------------------------------------------------- #
# Constraint (constraint.cpp): map the R -1 sentinels onto +/- infinity.
# --------------------------------------------------------------------------- #
class _Constraint:
    """Port of ``Constraint`` (constraint.cpp): ``-1`` -> ``+/-inf``, where R sets ``+/-DBL_MAX`` (module docstring)."""

    __slots__ = (
        "max_gap",
        "window_size",
        "age_min_begin",
        "age_max_begin",
        "age_max_end",
        "count_method",
    )

    def __init__(self, mg, ws, aminb, amaxb, amaxe, cmethod):
        self.max_gap = INF if mg == -1 else float(mg)
        self.window_size = INF if ws == -1 else float(ws)
        self.age_min_begin = -INF if aminb == -1 else float(aminb)
        self.age_max_begin = INF if amaxb == -1 else float(amaxb)
        self.age_max_end = INF if amaxe == -1 else float(amaxe)
        self.count_method = 1 if cmethod == -1 else int(cmethod)


# --------------------------------------------------------------------------- #
# Event-node chain (eventseq.cpp): Sequence / SequenceEventNode.
# --------------------------------------------------------------------------- #
class _Node:
    """One event in a sequence (``SequenceEventNode``): type, gap to previous, next."""

    __slots__ = ("type", "gap", "next")

    def __init__(self, t: int, g: float):
        self.type = int(t)
        self.gap = float(g)
        self.next = None  # type: Optional[_Node]

    def greater_than(self, e: int, g: float) -> bool:
        if self.gap > g:
            return True
        if self.gap == g and self.type > e:
            return True
        return False

    def add_event(self, e: int, t: float) -> None:
        current = self
        remaining = t
        while current.next is not None:
            delta = remaining - current.gap
            if current.next.greater_than(e, delta):
                new_node = _Node(e, delta)
                current.next.gap -= delta
                new_node.next = current.next
                current.next = new_node
                return
            current = current.next
            remaining = delta
        current.next = _Node(e, remaining - current.gap)

    def next_without_gap(self) -> Optional["_Node"]:
        if self.next is not None and self.next.gap == 0:
            return self.next
        return None

    def next_with_gap(self) -> Optional["_Node"]:
        # ``SequenceEventNode::getNextWithGap`` (eventseq.h:100) recurses once per
        # zero-gap node, so its depth is the length of a run of SIMULTANEOUS
        # events. C recurses on the machine stack; CPython would raise
        # RecursionError from inside seqedist(interval="next") -- its only caller
        # -- at roughly 1,000 simultaneous events, and sooner when the caller was
        # already deep. Walking the run keeps call-stack use constant, exactly as
        # the tree traversals below already do. Same node returned, same order.
        node = self
        while node.next is not None:
            if node.next.gap == 0:
                node = node.next
                continue
            return node.next
        return None


class _Sequence:
    """An event sequence (``Sequence``): an event-node chain + weight / obs_time."""

    __slots__ = ("idpers", "event", "dictionary", "obs_time", "weight", "generic")

    def __init__(self, idpers: int, dictionary: List[str], *, generic: bool = False):
        self.idpers = int(idpers)
        # A generic pattern (a subsequence), not a person's sequence. TraMineR's C++
        # marks it with idpers == -1, a valid id, so a sequence whose id is -1 lost its
        # gaps; the flag keeps the two apart (TX-AUD-20261019-006).
        self.generic = bool(generic)
        self.event = None  # type: Optional[_Node]
        self.dictionary = dictionary
        self.obs_time = -1.0
        self.weight = 1.0

    def is_generic(self) -> bool:
        return self.generic

    def add_event(self, e: int, t: float) -> None:
        if self.event is not None:
            if self.event.greater_than(e, t):
                self.event.gap -= t
                s = _Node(e, t)
                s.next = self.event
                self.event = s
            else:
                self.event.add_event(e, t)
        else:
            self.event = _Node(e, t)

    def copy(self) -> "_Sequence":
        s = _Sequence(self.idpers, self.dictionary, generic=self.generic)
        s.obs_time = self.obs_time
        s.weight = self.weight
        node = self.event
        if node is not None:
            head = _Node(node.type, node.gap)
            s.event = head
            node = node.next
            while node is not None:
                head.next = _Node(node.type, node.gap)
                head = head.next
                node = node.next
        return s

    def sprint(self) -> str:
        """Faithful ``Sequence::sprint`` -- gaps shown iff this is not generic."""
        if self.event is None:
            return ""
        return _sprint_node(
            self.event,
            start=True,
            print_gap=not self.is_generic(),
            dictionary=self.dictionary,
            remaining_time=self.obs_time,
        )

    # ----- the six count methods (Sequence::count switch) ------------------ #
    def count(
        self,
        s: "_Sequence",
        max_gap: float,
        window_size: float,
        age_min: float,
        age_max: float,
        age_max_end: float,
        c_method: int,
    ) -> float:
        return _seq_count(self, s, max_gap, window_size, age_min, age_max, age_max_end, c_method)

    def contain(self, event_set: set, exclude: bool) -> bool:
        if self.event is None:
            return False
        sen: Optional[_Node] = self.event
        while sen is not None:
            if sen.type in event_set:
                if not exclude:
                    return True
            elif exclude:
                return False
            sen = sen.next
        return exclude


def _fmt_gap(g: float) -> str:
    """Render a gap as R's ``format()`` does (``TMRNumberFormat``, eventseq.cpp:90).

    Seven significant digits, then fixed or scientific notation, whichever is
    narrower: ``1e+05`` and ``1e-04``, but ``123456`` and ``2.9999``.
    """
    return r_number_text(float(g), digits=7)


def _sprint_node(node: _Node, start, print_gap, dictionary, remaining_time) -> str:
    """Faithful ``SequenceEventNode::sprint`` (eventseq.cpp).

    ``remaining_time`` is the observation time, decremented by each node's gap as
    the C++ recursion descends (``next.sprint(..., remainingTime - this.gap)``); a
    positive value at the leaf appends a trailing gap to the end of observation.
    Generic subsequences pass ``print_gap=False`` (no gaps shown) and
    ``remaining_time=-1`` (no trailing gap).
    """
    parts: List[str] = []
    cur: Optional[_Node] = node
    is_start = start
    rem = remaining_time
    while cur is not None:
        lab = dictionary[cur.type - 1]
        if is_start:
            if cur.gap > 0 and print_gap:
                parts.append(f"{_fmt_gap(cur.gap)}-({lab}")
            else:
                parts.append(f"({lab}")
        elif cur.gap > 0:
            if print_gap:
                parts.append(f")-{_fmt_gap(cur.gap)}-({lab}")
            else:
                parts.append(f")-({lab}")
        else:
            parts.append(f",{lab}")
        is_start = False
        if cur.next is None:
            if rem is not None and rem > 0:
                parts.append(")-" + _fmt_gap(rem - cur.gap))
            else:
                parts.append(")")
        elif rem is not None:
            rem = rem - cur.gap
        cur = cur.next
    return "".join(parts)


# --------------------------------------------------------------------------- #
# Inner count matchers (SequenceEventNode::count / count3 / count4 / count5).
# --------------------------------------------------------------------------- #
def _iter_node_matches(
    this,
    s,
    max_gap,
    window_size,
    age_max_end,
    gap_consumed,
    current_age,
    consumed_limit=None,
):
    """Yield completed matches in the recursive C++ traversal order."""
    stack = [(this, s, gap_consumed, current_age)]
    while stack:
        pattern, sequence, consumed, age = stack.pop()
        if consumed_limit is not None and consumed > consumed_limit[0]:
            continue
        next_pattern = pattern.next
        if next_pattern is None:
            yield consumed, age
            continue

        sen = sequence.next
        candidates = []
        if next_pattern.gap == 0:
            while sen is not None and sen.gap == 0:
                if next_pattern.type == sen.type:
                    candidates.append((next_pattern, sen, consumed, age))
                sen = sen.next
        else:
            while sen is not None and sen.gap == 0:
                sen = sen.next
            g = 0.0
            while sen is not None:
                g += sen.gap
                if (
                    g > max_gap
                    or g + consumed > window_size
                    or age + g > age_max_end
                    or (consumed_limit is not None and g + consumed > consumed_limit[0])
                ):
                    break
                if next_pattern.type == sen.type:
                    candidates.append((next_pattern, sen, consumed + g, age + g))
                sen = sen.next

        stack.extend(reversed(candidates))


def _node_count(this, s, max_gap, window_size, age_max_end, gap_consumed, current_age) -> int:
    """``SequenceEventNode::count`` -- COBJ/CDIST_O inner matcher (eventseq.cpp:391)."""
    return sum(
        1
        for _ in _iter_node_matches(
            this, s, max_gap, window_size, age_max_end, gap_consumed, current_age
        )
    )


def _node_count3(
    this, s, max_gap, window_size, age_max_end, gap_consumed, current_age, Win, tWin, lWin
) -> int:
    """``SequenceEventNode::count3`` -- CWIN windowed marking (eventseq.cpp:437)."""
    c = 0
    for consumed, age in _iter_node_matches(
        this, s, max_gap, window_size, age_max_end, gap_consumed, current_age
    ):
        for i in range(lWin):
            if (tWin[i] + window_size) <= age_max_end:
                if age - consumed >= tWin[i] and age <= tWin[i] + window_size and Win[i] == 0:
                    Win[i] = 1
            else:
                if age - consumed >= tWin[i] and age <= age_max_end and Win[i] == 0:
                    Win[i] = 1
        c += 1
    return c


def _node_count4(this, s, max_gap, window_size, age_max_end, gap_consumed, current_age, min_win):
    """``SequenceEventNode::count4`` -- CMINWIN (eventseq.cpp:527). Returns (c, min_win)."""
    c = 0
    consumed_limit = [min_win]
    for consumed, _ in _iter_node_matches(
        this,
        s,
        max_gap,
        window_size,
        age_max_end,
        gap_consumed,
        current_age,
        consumed_limit,
    ):
        c += 1
        min_win = min(min_win, consumed)
        consumed_limit[0] = min_win
    return c, min_win


def _count5_event_index(event_type, age, type_sen, t_sen, length):
    """Mirror count5's index-zero sentinel lookup exactly."""
    index = 0
    cursor = 0
    while index == 0 and cursor < length:
        if t_sen[cursor] == age and type_sen[cursor] == event_type:
            index = cursor
        cursor += 1
    return index


def _count5_candidates(
    this, s, max_gap, window_size, age_max_end, gap_consumed, current_age, type_sen, t_sen, length
):
    """Return count5 child calls in the order used by eventseq.cpp."""
    next_pattern = this.next
    if next_pattern is None:
        return []

    sen = s.next
    candidates = []
    if next_pattern.gap == 0:
        while sen is not None and sen.gap == 0:
            if next_pattern.type == sen.type:
                index = _count5_event_index(next_pattern.type, current_age, type_sen, t_sen, length)
                candidates.append((sen, gap_consumed, current_age, index))
            sen = sen.next
    else:
        while sen is not None and sen.gap == 0:
            sen = sen.next
        g = 0.0
        while sen is not None:
            g += sen.gap
            if g > max_gap or g + gap_consumed > window_size or current_age + g > age_max_end:
                break
            if next_pattern.type == sen.type:
                age = current_age + g
                index = _count5_event_index(next_pattern.type, age, type_sen, t_sen, length)
                candidates.append((sen, gap_consumed + g, age, index))
            sen = sen.next
    return candidates


def _node_count5(
    this,
    s,
    max_gap,
    window_size,
    age_max_end,
    gap_consumed,
    current_age,
    typeSen,
    tSen,
    lSen,
    flagSen,
) -> int:
    """``SequenceEventNode::count5`` -- CDIST non-overlapping (eventseq.cpp:613)."""
    index = _count5_event_index(this.type, current_age, typeSen, tSen, lSen)
    if flagSen[index] != 0:
        return 0
    flagSen[index] = 1
    if this.next is None:
        for i in range(lSen):
            if flagSen[i] == 1:
                flagSen[i] = 2
        return 1

    def make_frame(pattern, sequence, consumed, age):
        return {
            "pattern": pattern,
            "candidates": _count5_candidates(
                pattern,
                sequence,
                max_gap,
                window_size,
                age_max_end,
                consumed,
                age,
                typeSen,
                tSen,
                lSen,
            ),
            "cursor": 0,
            "count": 0,
            "active_index": 0,
        }

    stack = [make_frame(this, s, gap_consumed, current_age)]
    returned = None
    while stack:
        if returned is not None:
            parent = stack[-1]
            parent["count"] += returned
            if flagSen[parent["active_index"]] == 2:
                returned = parent["count"]
                stack.pop()
            else:
                returned = None
            continue

        frame = stack[-1]
        if frame["cursor"] >= len(frame["candidates"]):
            for i in range(lSen):
                if flagSen[i] == 1:
                    flagSen[i] = 0
            returned = frame["count"]
            stack.pop()
            continue

        sen, consumed, age, child_index = frame["candidates"][frame["cursor"]]
        frame["cursor"] += 1
        if flagSen[child_index] != 0:
            continue

        next_pattern = frame["pattern"].next
        frame["active_index"] = child_index
        flagSen[child_index] = 1
        if next_pattern.next is None:
            for i in range(lSen):
                if flagSen[i] == 1:
                    flagSen[i] = 2
            returned = 1
        else:
            stack.append(make_frame(next_pattern, sen, consumed, age))

    return 0 if returned is None else returned


def _first_occurence(this_seq, s, max_gap, window_size, age_min, age_max, age_max_end) -> float:
    """``Sequence::first_occurence`` -- age of first occurrence, else -1 (eventseq.cpp:152)."""
    if this_seq.event is None or s.event is None:
        return -1.0
    age = 0.0
    sen = s.event
    while sen is not None:
        age += sen.gap
        if age > age_max:
            return -1.0
        if age >= age_min and this_seq.event.type == sen.type:
            if _node_count(this_seq.event, sen, max_gap, window_size, age_max_end, 0, age) > 0:
                return age
        sen = sen.next
    return -1.0


def _seq_count(this_seq, s, max_gap, window_size, age_min, age_max, age_max_end, c_method) -> float:
    """``Sequence::count`` -- the six-method dispatch (eventseq.cpp:169)."""
    if this_seq.event is None or s.event is None:
        return 0
    c = 0
    if c_method == 1 or c_method not in (2, 3, 4, 5, 6):
        # COBJ (and the default branch): 1 if the subsequence occurs at all.
        sen = s.event
        age = 0.0
        while sen is not None:
            age += sen.gap
            if age > age_max:
                break
            if age >= age_min and this_seq.event.type == sen.type:
                cs = _node_count(this_seq.event, sen, max_gap, window_size, age_max_end, 0, age)
                if cs > 0:
                    c = 1
            sen = sen.next
    elif c_method == 2:
        # CDIST_O: total number of distinct occurrences.
        sen = s.event
        age = 0.0
        while sen is not None:
            age += sen.gap
            if age > age_max:
                break
            if age >= age_min and this_seq.event.type == sen.type:
                c += _node_count(this_seq.event, sen, max_gap, window_size, age_max_end, 0, age)
            sen = sen.next
    elif c_method == 3:
        # CWIN: number of windows of width window_size containing the subsequence.
        if age_min == -INF:
            start_win = s.event.gap - window_size
        else:
            start_win = age_min
        if age_max_end == INF:
            sen = s.event
            age = 0.0
            while sen is not None:
                age += sen.gap
                sen = sen.next
            stop_win = age + window_size
        else:
            stop_win = age_max_end
        l_win = int(round(stop_win - start_win - window_size + 1))
        if l_win < 0:
            l_win = 0
        tWin = [start_win + i for i in range(l_win)]
        Win = [0] * l_win
        sen = s.event
        age = 0.0
        while sen is not None:
            age += sen.gap
            if age > age_max:
                break
            if age >= age_min and this_seq.event.type == sen.type:
                _node_count3(
                    this_seq.event, sen, max_gap, window_size, age_max_end, 0, age, Win, tWin, l_win
                )
            sen = sen.next
        c = sum(Win)
    elif c_method == 4:
        # CMINWIN: occurrences within the minimum-width window.
        min_win = INF
        min_win_old = min_win
        sen = s.event
        age = 0.0
        while sen is not None:
            age += sen.gap
            if age > age_max:
                break
            if age >= age_min and this_seq.event.type == sen.type:
                cs, min_win = _node_count4(
                    this_seq.event, sen, max_gap, window_size, age_max_end, 0, age, min_win
                )
                if min_win_old == min_win and cs > 0:
                    c += 1
                if min_win_old > min_win:
                    c = 1
                    min_win_old = min_win
            sen = sen.next
    elif c_method == 5:
        # CDIST: distinct occurrences with no event-timestamp overlap.
        l_sen = 0
        sen = s.event
        while sen is not None:
            l_sen += 1
            sen = sen.next
        tSen = [0.0] * l_sen
        typeSen = [0] * l_sen
        sen = s.event
        age = 0.0
        for i in range(l_sen):
            age += sen.gap
            tSen[i] = age
            typeSen[i] = sen.type
            sen = sen.next
        flagSen = [0] * l_sen
        sen = s.event
        age = 0.0
        while sen is not None:
            age += sen.gap
            if age > age_max:
                break
            if age >= age_min and this_seq.event.type == sen.type:
                c += _node_count5(
                    this_seq.event,
                    sen,
                    max_gap,
                    window_size,
                    age_max_end,
                    0,
                    age,
                    typeSen,
                    tSen,
                    l_sen,
                    flagSen,
                )
            sen = sen.next
    elif c_method == 6:
        # age of first appearance. C++ Sequence::count accumulates into an `int c`,
        # so the (double) first-occurrence age is truncated toward zero (and the
        # -1.0 "not found" sentinel stays -1). int() truncates toward zero to match.
        c += int(_first_occurence(this_seq, s, max_gap, window_size, age_min, age_max, age_max_end))
    return c


# --------------------------------------------------------------------------- #
# Prefix tree (treeeventnode.cpp / treeeventmap.cpp / prefixtree.cpp).
# --------------------------------------------------------------------------- #
class _NodeCounter:
    """Live-node count of one prefix tree (R's ``TreeEventNode::nodeCount``).

    R keeps one static counter, which only a single-threaded caller can share;
    each tree owns its count here so that concurrent mining calls cannot move
    one another's stop rule (TX-AUD-20261011-021).
    """

    __slots__ = ("value",)

    def __init__(self) -> None:
        self.value = 0


class _TreeEventNode:
    """Port of ``TreeEventNode``: type, support, lastID, brother & child maps."""

    __slots__ = ("type", "support", "lastID", "brother", "child", "counter")

    def __init__(self, t: int, counter: _NodeCounter):
        self.type = int(t)
        self.support = 0.0
        # The id of the last sequence counted; None, not -1 as in the C++ (an id may
        # be -1: TX-AUD-20261019-006).
        self.lastID: int | None = None
        self.brother = {}  # type: Dict[int, _TreeEventNode]
        self.child = {}  # type: Dict[int, _TreeEventNode]
        self.counter = counter
        counter.value += 1

    def add_sequence_internal(self, s, en, cst, gap_consumed, current_age, k, current_k):
        # Each frame mirrors one C++ addSequenceInternal invocation. Re-adding
        # the parent before its descendant preserves the recursive depth-first
        # traversal while keeping call-stack use constant.
        stack = [(self, en, gap_consumed, current_age, current_k, False, en, 0.0)]
        while stack:
            node, start, consumed, age, depth, entered, cursor, current_gap = stack.pop()
            if not entered:
                if cst.count_method == 2 or node.lastID != s.idpers:
                    node.support += s.weight
                    node.lastID = s.idpers
                if start.next is None or depth > k:
                    continue

            next_node = cursor.next
            if next_node is None:
                continue
            next_gap = current_gap + next_node.gap
            if (
                consumed + next_gap > cst.window_size
                or next_gap > cst.max_gap
                or next_gap + age > cst.age_max_end
            ):
                continue

            stack.append((node, start, consumed, age, depth, True, next_node, next_gap))
            branches = node.child if next_gap > 0 else node.brother
            descendant = branches.get(next_node.type)
            if descendant is None and k == depth:
                descendant = _TreeEventNode(next_node.type, node.counter)
                branches[next_node.type] = descendant
            if descendant is not None:
                stack.append(
                    (
                        descendant,
                        next_node,
                        consumed + next_gap,
                        age + next_gap,
                        depth + 1,
                        False,
                        next_node,
                        0.0,
                    )
                )

    def simplify_tree(self, min_sup):
        _simplify_map(self.brother, min_sup)
        _simplify_map(self.child, min_sup)

    def clear_support(self):
        stack = [self]
        while stack:
            node = stack.pop()
            node.support = 0.0
            node.lastID = None
            stack.extend(node.child.values())
            stack.extend(node.brother.values())

    def get_subsequences(self, results, supports, s, step, ed):
        _get_subseq_map(self.brother, results, supports, s, step, ed)
        _get_subseq_map(self.child, results, supports, s, step + 1, ed)


def _free_subtree(node: _TreeEventNode) -> None:
    """Mirror the C++ destructor's node-count decrement without recursion."""
    stack = [node]
    while stack:
        current = stack.pop()
        stack.extend(current.brother.values())
        stack.extend(current.child.values())
        current.brother.clear()
        current.child.clear()
        current.counter.value -= 1


def _simplify_map(m: Dict[int, _TreeEventNode], min_sup: float) -> None:
    pending = [m]
    while pending:
        current_map = pending.pop()
        retained = []
        for key in list(current_map.keys()):
            node = current_map[key]
            if node.support < min_sup:
                _free_subtree(node)
                del current_map[key]
            else:
                retained.append(node)
        for node in reversed(retained):
            pending.append(node.child)
            pending.append(node.brother)


def _get_subseq_map(m: Dict[int, _TreeEventNode], results, supports, s2, step, ed) -> None:
    # Stack insertion is reversed so pop order matches recursive std::map
    # traversal: each node, its brother map, its child map, then the next node.
    stack = [(m[key], s2, step) for key in reversed(sorted(m.keys()))]
    while stack:
        node, parent, node_step = stack.pop()
        sequence = parent.copy() if parent is not None else _Sequence(-1, ed, generic=True)
        sequence.add_event(node.type, node_step)
        results.append(sequence)
        supports.append(node.support)
        for key in reversed(sorted(node.child.keys())):
            stack.append((node.child[key], sequence, node_step + 1))
        for key in reversed(sorted(node.brother.keys())):
            stack.append((node.brother[key], sequence, node_step))


class _PrefixTree:
    """Port of ``PrefixTree``: the top-level event map + the mining loop."""

    __slots__ = ("child", "counter")

    def __init__(self):
        self.child = {}  # type: Dict[int, _TreeEventNode]
        self.counter = _NodeCounter()

    def add_sequence(self, s: _Sequence, cst: _Constraint, k: int) -> None:
        if s.event is None:
            return
        e: Optional[_Node] = s.event
        age = 0.0
        while e is not None:
            age += e.gap
            if age > cst.age_max_begin:
                break
            if age >= cst.age_min_begin:
                ten = self.child.get(e.type)
                if ten is not None:
                    ten.add_sequence_internal(s, e, cst, 0, age, k, 2)
                elif k == 1:
                    ten = _TreeEventNode(e.type, self.counter)
                    self.child[e.type] = ten
                    ten.add_sequence_internal(s, e, cst, 0, age, k, 2)
            e = e.next

    def simplify_tree(self, min_sup: float) -> None:
        _simplify_map(self.child, min_sup)

    def clear_support(self) -> None:
        for n in self.child.values():
            n.clear_support()

    def get_subsequences(self, ed) -> Tuple[List[_Sequence], List[float]]:
        results = []  # type: List[_Sequence]
        supports = []  # type: List[float]
        _get_subseq_map(self.child, results, supports, None, 0, ed)
        return results, supports


# --------------------------------------------------------------------------- #
# Driver: tmrfindsubsequences (tmrsequence.cpp).
# --------------------------------------------------------------------------- #
def tmrfindsubsequences(
    seqs,
    max_gap,
    window_size,
    age_min_begin,
    age_max_begin,
    age_max_end,
    count_method,
    min_support,
    max_k,
):
    """Port of ``tmrfindsubsequences`` -- mine + recount frequent subsequences.

    Returns ``(supports, counts, subsequences)`` where ``subsequences`` is a list
    of generic ``_Sequence`` objects, ``counts[j]`` = sum over seqs of
    ``weight*count(subseq_j, cMethod)``, ``supports[j]`` = weighted number of
    sequences with a positive count.
    """
    cst = _Constraint(max_gap, window_size, age_min_begin, age_max_begin, age_max_end, count_method)
    m_gap = cst.max_gap
    w_size = cst.window_size
    a_min = cst.age_min_begin
    a_max = cst.age_max_begin
    a_max_end = cst.age_max_end
    c_method = int(count_method)
    max_k_eff = 2**31 - 1 if max_k == -1 else int(max_k)
    m_support = float(min_support)

    root = _PrefixTree()
    ed = seqs[0].dictionary if seqs else []
    try:
        k = 1
        while k <= max_k_eff:
            root.clear_support()
            last_node_count = root.counter.value
            for s in seqs:
                root.add_sequence(s, cst, k)
            root.simplify_tree(m_support)
            if root.counter.value - last_node_count == 0:
                break
            k += 1

        subsequences, supports = root.get_subsequences(ed)
        counts = [0.0] * len(subsequences)
        for j, sub in enumerate(subsequences):
            c = 0.0
            sp = 0.0
            for s in seqs:
                counting = s.weight * sub.count(s, m_gap, w_size, a_min, a_max, a_max_end, c_method)
                c += counting
                if counting > 0:
                    sp += s.weight
            counts[j] = c
            supports[j] = sp
        return supports, counts, subsequences
    finally:
        for node in list(root.child.values()):
            _free_subtree(node)
        root.child.clear()


# --------------------------------------------------------------------------- #
# Count matrix: tmrmatrixsubseqinseq (tmrsequence.cpp) -> seqeapplysub.
# --------------------------------------------------------------------------- #
def tmrmatrixsubseqinseq(
    subseqs, seqs, max_gap, window_size, age_min_begin, age_max_begin, age_max_end, count_method
):
    """Port of ``tmrmatrixsubseqinseq`` -- the ``(n_seq x n_subseq)`` count matrix.

    Faithful to the C++: this maps the ``-1`` sentinels inline (and does *not*
    remap ``count_method``); counts are unweighted.
    """
    w_size = INF if window_size == -1 else float(window_size)
    m_gap = INF if max_gap == -1 else float(max_gap)
    a_min = -INF if age_min_begin == -1 else float(age_min_begin)
    a_max = INF if age_max_begin == -1 else float(age_max_begin)
    a_max_end = INF if age_max_end == -1 else float(age_max_end)
    c_method = int(count_method)

    ns = len(seqs)
    nsub = len(subseqs)
    matrix = np.zeros((ns, nsub), dtype=float)
    for i, sub in enumerate(subseqs):
        for j, s in enumerate(seqs):
            matrix[j, i] = sub.count(s, m_gap, w_size, a_min, a_max, a_max_end, c_method)
    return matrix


# --------------------------------------------------------------------------- #
# Event-sequence distance: tmrseqedist (seqedist.cpp).
# --------------------------------------------------------------------------- #
def _indel_subst_cost(sen1: _Node, sen2: _Node, idcost) -> Tuple[float, float, float]:
    """Port of ``getIndelSubstitutionCost`` (seqedist.cpp): indel1, indel2, subst."""
    indel1 = indel2 = subst = 0.0
    s1: Optional[_Node] = sen1
    while s1 is not None:
        indel1 += idcost[s1.type - 1]
        s1 = s1.next_without_gap()
    s2: Optional[_Node] = sen2
    while s2 is not None:
        indel2 += idcost[s2.type - 1]
        s2 = s2.next_without_gap()

    s1 = sen1
    s2 = sen2
    while True:
        assert s1 is not None and s2 is not None
        if s1.type == s2.type:
            s1 = s1.next_without_gap()
            s2 = s2.next_without_gap()
        elif s1.type < s2.type:
            subst += idcost[s1.type - 1]
            s1 = s1.next_without_gap()
        else:
            subst += idcost[s2.type - 1]
            s2 = s2.next_without_gap()
        if s1 is None:
            while s2 is not None:
                subst += idcost[s2.type - 1]
                s2 = s2.next_without_gap()
            return indel1, indel2, subst
        if s2 is None:
            while s1 is not None:
                subst += idcost[s1.type - 1]
                s1 = s1.next_without_gap()
            return indel1, indel2, subst


def tmrseqedist(seqs, idcost, vparam, norm, interval):
    """Port of ``tmrseqedist`` (seqedist.cpp) -- the event-sequence edit distance."""
    V = float(vparam)
    interval = int(interval)
    norm = int(norm)
    ns = len(seqs)
    matrix = np.zeros((ns, ns), dtype=float)

    def gap_groups(seq: _Sequence):
        """Heads of each transition (gap>0 boundary), like ``getNextWithGap`` walks."""
        heads = []
        node = seq.event
        first = True
        while node is not None:
            if first or node.gap > 0:
                heads.append(node)
                first = False
            node = node.next
        return heads

    for i in range(ns):
        s = seqs[i]
        matrix[i, i] = 0.0
        s_heads = gap_groups(s)
        ni = len(s_heads)
        for j in range(i + 1, ns):
            s2 = seqs[j]
            s2_heads = gap_groups(s2)
            nj = len(s2_heads)
            maxevent = max(ni, nj) + 2
            Fmat = np.zeros((maxevent, maxevent), dtype=float)
            maxcost = 0.0
            # First column: indel of each transition of s.
            ei = 1
            for head in s_heads:
                cost = 0.0
                node = head
                while node is not None:
                    cost += idcost[node.type - 1]
                    node = node.next_without_gap()
                maxcost += cost
                Fmat[ei, 0] = Fmat[ei - 1, 0] + cost
                ei += 1
            # First row: indel of each transition of s2.
            ej = 1
            for head in s2_heads:
                cost = 0.0
                node = head
                while node is not None:
                    cost += idcost[node.type - 1]
                    node = node.next_without_gap()
                maxcost += cost
                Fmat[0, ej] = Fmat[0, ej - 1] + cost
                ej += 1
            # DP fill.
            ei = 1
            age1 = 0.0
            for ihead, sen in enumerate(s_heads):
                if interval == 2:
                    age1 = sen.gap
                elif interval == 3:
                    nxt = sen.next_with_gap()
                    if nxt is not None:
                        age1 += nxt.gap
                    # else: the C++ writes obsTime-age1 to a temp that is then
                    # overwritten by `tmpage1 = age1`, so age1 is left unchanged.
                else:
                    age1 += sen.gap
                tmpage1 = age1
                age2 = 0.0
                ej = 1
                for jhead, sen2 in enumerate(s2_heads):
                    if interval == 2:
                        age2 = sen2.gap
                    elif interval == 3:
                        nxt2 = sen2.next_with_gap()
                        if nxt2 is not None:
                            age2 += nxt2.gap
                        # else: dead obsTime temp in the C++ (tmpage2 unused) ->
                        # age2 unchanged.
                    else:
                        age2 += sen2.gap
                    indel1, indel2, scost = _indel_subst_cost(sen, sen2, idcost)
                    scosttot = indel1 + indel2
                    up = Fmat[ei - 1, ej] + indel1
                    left = Fmat[ei, ej - 1] + indel2
                    kcost = V * abs(tmpage1 - age2) + scost
                    if kcost > scosttot:
                        kcost = scosttot
                    kcost += Fmat[ei - 1, ej - 1]
                    if interval == 1:
                        age2 = tmpage1
                    if kcost > up:
                        kcost = up
                    if kcost > left:
                        kcost = left
                    Fmat[ei, ej] = kcost
                    ej += 1
                ei += 1
            dist = Fmat[ei - 1, ej - 1]
            if norm > 0:
                if norm == 1:
                    if maxcost > 0:
                        dist /= maxcost
                    elif dist > 0:
                        dist = 1.0
                elif norm == 2 and dist != 0:
                    # A raw distance of 0 stays 0, as TraMineR's own normaliser
                    # returns (distancecalculator.h: if (rawdist==0) return 0).
                    # seqedist.cpp:231 divides 0 by 0 instead and returns NaN for
                    # two empty sequences or an all-zero idcost: an R defect
                    # corrected here (TX-AUD-20261010-037; ledger seqedist row).
                    dist = (2 * dist) / (maxcost + dist)
            matrix[i, j] = dist
            matrix[j, i] = dist
    return matrix


# Integers up to 2**52 in magnitude subtract exactly in binary64, so for such
# times every order of insertion produces the same chain.
_EXACT_INTEGER_LIMIT = float(2**52)


def eseq_to_engine(eseq) -> _Sequence:
    """Build an engine ``_Sequence`` from an RP5a :class:`~traminex.events._object.Eseq`.

    TraMineR builds the node chain with ``Sequence::addEvent``, inserting the
    events one at a time in their input order (``tmrsequence.cpp:113-127``,
    ``eventseq.cpp:21-48``); each insertion subtracts the gaps it walks past, so
    with fractional times the stored gaps -- and the ages re-added from them --
    carry that order's rounding (``1.9999999999999998`` rather than ``2``).
    The events are replayed in the order :func:`seqecreate` recorded; an
    ``Eseq`` built by hand is replayed in its array order, as R's
    ``C_tmrsequence`` inserts the vectors it is given.

    When every time is a whole number no subtraction rounds, so every insertion
    order yields the chain of consecutive differences of the (time, code)-sorted
    events; that chain is built directly. Weight / obs_time are carried over.
    """
    s = _Sequence(int(eseq.id_code), list(eseq.dictionary))
    s.weight = float(eseq.weight)
    s.obs_time = float(eseq.obs_time)
    times = np.asarray(eseq.times, dtype=float)
    codes = np.asarray(eseq.codes, dtype=int)
    count = len(times)
    exact = bool(
        np.all(np.isfinite(times))
        and np.all(np.abs(times) <= _EXACT_INTEGER_LIMIT)
        and np.all(times == np.round(times))
    )
    if not exact:
        order = getattr(eseq, "_insertion_order", None)
        if order is None or not np.array_equal(np.sort(order), np.arange(count)):
            order = np.arange(count)
        for k in order:
            s.add_event(int(codes[k]), float(times[k]))
        return s

    head = None  # type: Optional[_Node]
    tail = None  # type: Optional[_Node]
    previous = 0.0
    for k in np.lexsort((codes, times)):
        time = float(times[k])
        node = _Node(int(codes[k]), time if head is None else time - previous)
        if head is None:
            head = node
        else:
            assert tail is not None
            tail.next = node
        tail = node
        previous = time
    s.event = head
    return s
