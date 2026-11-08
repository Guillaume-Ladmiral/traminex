"""Private TraMineR-compatible state-order canonicalization.

This module centralizes the quirks shared by ``seqprecstart``, ``seqibad``,
``seqidegrad``, and the ranked ``seqindic`` routes.  It intentionally remains
private: the public API is the analytical functions that consume the contract.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .sequences import STSList
from .sequences import alphabet as sequence_alphabet


@dataclass(frozen=True)
class StateOrderContract:
    """Immutable, alphabet-aligned result of TraMineR state canonicalization.

    ``groups`` contains comparable groups in rank order.  A group retains every
    member of its R equivalence class so consumers can collapse transitions to
    its first member.  ``code_to_group`` maps comparable alphabet positions to
    those groups and is ``-1`` for incomparable states, including an
    incomparable member retained in a group solely for equivalence collapse.
    ``precedence`` is aligned with ``labels``.
    """

    labels: tuple[str, ...]
    groups: tuple[tuple[int, ...], ...]
    precedence: np.ndarray
    code_to_group: np.ndarray


def _normalize_equivalence(
    state_equiv: list[list[str]] | None,
    labels: tuple[str, ...],
) -> list[list[str]]:
    if state_equiv is None:
        return []
    if not isinstance(state_equiv, list) or any(
        isinstance(group, str) or not isinstance(group, (list, tuple, np.ndarray))
        for group in state_equiv
    ):
        raise TypeError("state_equiv must be a list of equivalence classes")
    groups = [list(group) for group in state_equiv]
    flat = [state for group in groups for state in group]
    if len(flat) != len(set(flat)):
        raise ValueError("Multiple occurrence of same state in state_equiv")
    unknown = [state for state in flat if state not in labels]
    if unknown:
        raise ValueError(f"Bad state_equiv, states not in the alphabet: {unknown}")
    return groups


def _check_and_promote_order(
    x: STSList,
    order: list[str | None],
    equiv: list[list[str]],
    labels: tuple[str, ...],
) -> list[str | None]:
    if len(order) != len(set(order)):
        raise ValueError("Multiple occurrences of same state in state_order")
    unknown = [state for state in order if state is not None and state not in labels]
    if unknown:
        raise ValueError(f"Bad state_order, states not in the alphabet: {unknown}")

    # states.check uses alphabet(seqdata), without missing, for this pass even
    # when with.missing=TRUE.  Its inner test examines only the first listed
    # member of a class, and every insertion occurs at that same position.
    ordinary = list(sequence_alphabet(x, with_missing=False))
    initially_noncomp = [state for state in ordinary if state not in order]
    for state in initially_noncomp:
        for group in equiv:
            if state not in group:
                continue
            if group and group[0] in order:
                insert_at = order.index(group[0]) + 1
                order[insert_at:insert_at] = [state]
            break
    return order


def _explicit_order(
    raw: np.ndarray,
    labels: tuple[str, ...],
    state_equiv: list[list[str]] | None,
) -> tuple[list[str | None], list[list[str]], np.ndarray]:
    ordered = np.argsort(raw, kind="stable")
    negative_count = int(np.count_nonzero(raw < 0))
    if negative_count == raw.size:
        # R's `(1 + n):n` descending-index quirk selects an out-of-range NA and
        # the final ordered state when every value is negative.
        order: list[str | None] = [None, labels[int(ordered[-1])]]
    else:
        order = [labels[int(position)] for position in ordered[negative_count:]]

    values, counts = np.unique(raw, return_counts=True)
    derived_equiv = [
        [labels[index] for index in range(len(labels)) if raw[index] == value]
        for value, count in zip(values, counts)
        if count > 1
    ]
    # Source quirk: with no explicit-cost ties, caller state.equiv survives the
    # warning that says state.order and state.equiv were overridden.
    equiv = derived_equiv if derived_equiv else _normalize_equivalence(state_equiv, labels)

    normalized = raw.copy()
    negative = normalized < 0
    if np.any(negative):
        if negative_count == raw.size:
            replacement = np.nan
        else:
            replacement = float(np.mean(normalized[ordered[negative_count:]]))
        normalized[negative] = replacement
    with np.errstate(invalid="ignore", divide="ignore"):
        normalized = (normalized - np.min(normalized)) / (np.max(normalized) - np.min(normalized))
    return order, equiv, normalized


def _canonical_groups(
    labels: tuple[str, ...],
    order: list[str | None],
    equiv: list[list[str]],
) -> tuple[tuple[tuple[int, ...], ...], np.ndarray]:
    groups: list[tuple[int, ...]] = []
    assigned: set[str] = set()
    for state in order:
        if state is None or state in assigned:
            continue
        containing = next((group for group in equiv if state in group), None)
        members = list(containing) if containing is not None else [state]
        if containing is not None and members and members[0] in order and state != members[0]:
            continue
        members = [member for member in members if member in labels]
        groups.append(tuple(labels.index(member) for member in members))
        assigned.update(members)

    code_to_group = np.full(len(labels), -1, dtype=int)
    comparable = {state for state in order if state is not None}
    for group_index, group in enumerate(groups):
        for code in group:
            if labels[code] in comparable:
                code_to_group[code] = group_index
    return tuple(groups), code_to_group


def _resolve_state_order(
    x: STSList,
    *,
    state_order: Sequence[str] | None = None,
    state_equiv: list[list[str]] | None = None,
    stprec: Sequence[float] | np.ndarray | None = None,
    with_missing: bool = False,
) -> StateOrderContract:
    """Resolve TraMineR ``states.check``/``seqprecstart`` into one contract."""
    if not isinstance(x, STSList):
        raise TypeError("x must be a state sequence object")
    labels = tuple(sequence_alphabet(x, with_missing=with_missing))

    if stprec is not None:
        supplied = np.asarray(stprec, dtype=float).reshape(-1).copy()
        if supplied.size != len(labels):
            raise ValueError("length(stprec) should equal length of alphabet")
        if np.any(np.isnan(supplied)):
            raise ValueError("NA/NaN argument")
        order, equiv, precedence = _explicit_order(supplied, labels, state_equiv)
    else:
        if state_order is None:
            order = list(labels)
        elif isinstance(state_order, str):
            order = [state_order]
        else:
            order = list(state_order)
        equiv = _normalize_equivalence(state_equiv, labels)
        precedence = np.empty(len(labels), dtype=float)

    order = _check_and_promote_order(x, order, equiv, labels)

    if stprec is None:
        if not order:
            raise ValueError("wrong sign in 'by' argument")
        if len(order) == 1:
            ranked_costs = np.asarray([np.nan])
        else:
            ranked_costs = np.linspace(0.0, 1.0, len(order))
        noncomp = [state for state in labels if state not in order]
        noncomp_cost = float(np.mean(ranked_costs))
        by_state = {
            state: float(cost) for state, cost in zip(order, ranked_costs) if state is not None
        }
        by_state.update({state: noncomp_cost for state in noncomp})
        precedence = np.asarray([by_state[state] for state in labels], dtype=float)

    for group in equiv:
        positions = [labels.index(state) for state in group]
        if positions:
            precedence[positions] = float(np.mean(precedence[positions]))

    groups, code_to_group = _canonical_groups(labels, order, equiv)
    immutable_precedence = np.asarray(precedence, dtype=float).copy()
    immutable_mapping = np.asarray(code_to_group, dtype=int).copy()
    immutable_precedence.setflags(write=False)
    immutable_mapping.setflags(write=False)
    return StateOrderContract(
        labels=labels,
        groups=groups,
        precedence=immutable_precedence,
        code_to_group=immutable_mapping,
    )
