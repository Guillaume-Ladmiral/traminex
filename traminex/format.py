from __future__ import annotations

import math
from numbers import Integral
from typing import Any, Mapping, Optional, Sequence, Union, cast

import numpy as np
import pandas as pd

from ._categorical_order import r_character_rank, r_seqstatl_order
from ._rformat import r_as_character, r_number_text
from ._sequence_builders import (
    _default_seqdef_palette,
    build_sts_from_state_matrix,
    validate_seq_kwargs,
)
from .seqformat import (
    _cell_is_na,
    _collision_free_marker,
    _dataframe_to_r_matrix,
    _r_state_text,
    _r_strsplit,
    seqdef_from_long,
)
from .seqformat import (
    seqformat as _seqformat,
)
from .sequences import STSList, _require_stslist, _validate_sequence_weights


def _as_rectangular_object(data: object) -> np.ndarray:
    def padded_rows() -> np.ndarray:
        if isinstance(data, (str, bytes)):
            raise ValueError("seqdef expects a 2D table of states")
        try:
            source_rows = list(cast(Sequence[Sequence[object]], data))
            if any(isinstance(row, (str, bytes)) for row in source_rows):
                raise TypeError
            rows = [list(row) for row in source_rows]
        except (TypeError, ValueError):
            raise ValueError("seqdef expects a 2D table of states") from None
        # A cell that is itself a sequence (rows given as matrices, say) is not a state;
        # without this, its printed form became a state (2026-11-08, gap part 4).
        if any(isinstance(cell, (list, tuple, np.ndarray)) for row in rows for cell in row):
            raise ValueError("seqdef expects a 2D table of states")
        width = max((len(row) for row in rows), default=0)
        padded: np.ndarray = np.full((len(rows), width), None, dtype=object)
        for row_index, row in enumerate(rows):
            padded[row_index, : len(row)] = row
        return padded

    try:
        matrix = np.asarray(data, dtype=object)
    except ValueError:
        matrix = padded_rows()
    if matrix.ndim == 1:
        matrix = padded_rows()
    if matrix.ndim != 2:
        raise ValueError("seqdef expects a 2D table of states")
    return matrix


def _selectors(var: Sequence[int | str] | int | str | None) -> list[int | str] | None:
    if var is None:
        return None
    if isinstance(var, str):
        return [var]
    if isinstance(var, Integral) and not isinstance(var, bool):
        return [int(var)]
    return list(cast(Sequence[int | str], var))


def _extract_seqdef_table(
    data: object,
    var: Sequence[int | str] | int | str | None,
) -> tuple[
    np.ndarray,
    list[object] | None,
    list[str] | None,
    tuple[np.ndarray, np.ndarray] | None,
]:
    selected = _selectors(var)
    source_marker_masks: tuple[np.ndarray, np.ndarray] | None = None
    if isinstance(data, STSList):
        matrix = np.empty(data.data.shape, dtype=object)
        for index, code in np.ndenumerate(data.data):
            matrix[index] = data._state_for_code(int(code))
        nr_mask = (
            np.zeros(data.data.shape, dtype=bool)
            if data.nr_code is None
            else data.data == data.nr_code
        )
        void_mask = (
            np.zeros(data.data.shape, dtype=bool)
            if data.void_code is None
            else data.data == data.void_code
        )
        columns = None if data.xtlab is None else list(data.xtlab)
        ids = None if data.ids is None else list(data.ids)
        if selected is not None:
            sts_positions: list[int] = []
            sts_ncols = cast(tuple[int, int], matrix.shape)[1]
            for selector in selected:
                if isinstance(selector, str):
                    if columns is None or selector not in columns:
                        raise ValueError("undefined columns selected")
                    sts_positions.append(columns.index(selector))
                elif -sts_ncols <= selector < sts_ncols:
                    sts_positions.append(selector % sts_ncols)
                else:
                    raise ValueError("undefined columns selected")
            matrix = matrix[:, sts_positions]
            nr_mask = nr_mask[:, sts_positions]
            void_mask = void_mask[:, sts_positions]
            if columns is not None:
                columns = [columns[position] for position in sts_positions]
        source_marker_masks = (nr_mask, void_mask)
    elif isinstance(data, pd.DataFrame):
        frame = data
        if selected is not None:
            columns_to_use: list[object] = []
            for selector in selected:
                if isinstance(selector, str):
                    if selector not in frame.columns:
                        raise ValueError("undefined columns selected")
                    columns_to_use.append(selector)
                elif -frame.shape[1] <= selector < frame.shape[1]:
                    columns_to_use.append(frame.columns[selector])
                else:
                    raise ValueError("undefined columns selected")
            frame = frame.loc[:, columns_to_use]
        matrix = _dataframe_to_r_matrix(frame)
        columns = [str(column) for column in frame.columns]
        default_index = isinstance(frame.index, pd.RangeIndex) and frame.index.equals(
            pd.RangeIndex(len(frame))
        )
        ids = None if default_index else list(frame.index)
    else:
        matrix = _as_rectangular_object(data)
        if selected is not None:
            array_positions: list[int] = []
            ncols = cast(tuple[int, int], matrix.shape)[1]
            for selector in selected:
                if not isinstance(selector, int) or isinstance(selector, bool):
                    raise ValueError("undefined columns selected")
                if selector < -ncols or selector >= ncols:
                    raise ValueError("undefined columns selected")
                array_positions.append(selector % ncols)
            matrix = matrix[:, array_positions]
        columns = None
        ids = None
    return matrix, ids, columns, source_marker_masks


def _render_seqdef_matrix(matrix: np.ndarray) -> np.ndarray:
    rendered = np.empty(matrix.shape, dtype=object)
    for index, value in np.ndenumerate(matrix):
        rendered[index] = None if _cell_is_na(value) else _r_state_text(value)
    return rendered


def _decompress_seqdef_sts(
    matrix: np.ndarray,
    *,
    separator: str,
) -> tuple[np.ndarray, list[str]]:
    observed = [value for value in matrix.ravel(order="C") if not _cell_is_na(value)]
    if any(not isinstance(value, str) for value in observed):
        raise ValueError("non-character argument")
    if matrix.shape[1] != 1:
        raise ValueError("compressed STS input must contain exactly one column")
    rows: list[list[object]] = []
    for value in matrix[:, 0]:
        if _cell_is_na(value):
            rows.append([None])
            continue
        parts = _r_strsplit(str(value), separator)
        rows.append([None if part == "NA" else part for part in parts])
    width = max((len(row) for row in rows), default=0)
    if width == 0:
        raise ValueError("length of 'dimnames' [2] not equal to array extent")
    out: np.ndarray = np.full((len(rows), width), None, dtype=object)
    for row_index, row in enumerate(rows):
        out[row_index, : len(row)] = row
    return out, [f"[{position}]" for position in range(1, width + 1)]


def _stable_states(matrix: np.ndarray, excluded: np.ndarray) -> list[str]:
    states: list[str] = []
    seen: set[str] = set()
    for index, value in np.ndenumerate(matrix):
        if excluded[index] or value is None:
            continue
        state = str(value)
        if state not in seen:
            seen.add(state)
            states.append(state)
    return states


def _state_option(value: object) -> str | None:
    if value is None or _cell_is_na(value):
        return None
    return _r_state_text(value)


def _apply_missing_policies(
    matrix: np.ndarray,
    missing_mask: np.ndarray,
    nr_mask: np.ndarray,
    void_mask: np.ndarray,
    *,
    left: str | None,
    gaps: str | None,
    right: str | None,
    void: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not bool(missing_mask.any()):
        return matrix, nr_mask.copy(), void_mask.copy()
    out = matrix.copy()
    out_nr = nr_mask.copy()
    out_void = void_mask.copy()
    width = matrix.shape[1]

    def replace(row_index: int, positions: slice | np.ndarray, policy: str | None) -> None:
        if policy is None:
            return
        if policy == "NEUTRAL":
            raise ValueError("'NEUTRAL' missing policies cannot be encoded in an STSList")
        if policy == "DEL":
            out[row_index, positions] = void
            out_void[row_index, positions] = True
        else:
            out[row_index, positions] = policy
            out_void[row_index, positions] = False
        out_nr[row_index, positions] = False

    for row_index in range(matrix.shape[0]):
        row_missing = missing_mask[row_index]
        if not bool(row_missing.any()):
            continue
        observed = np.flatnonzero(~row_missing)
        if observed.size == 0:
            if right is not None:
                replace(row_index, slice(None), right)
        else:
            first = int(observed[0])
            last = int(observed[-1])
            if first > 0 and left is not None:
                replace(row_index, slice(None, first), left)
            if last + 1 < width and right is not None:
                replace(row_index, slice(last + 1, None), right)
            if gaps is not None and last - first > 1:
                interior = np.arange(first + 1, last)
                gap_positions = interior[row_missing[interior]]
                replace(row_index, gap_positions, gaps)

        # seqprep moves structural voids to the tail.  Follow provenance rather
        # than marker spelling so a code-aware literal ``%`` remains a real state.
        structural_void = out_void[row_index]
        void_count = int(structural_void.sum())
        if void_count:
            kept = ~structural_void
            kept_values = list(out[row_index, kept])
            kept_nr = out_nr[row_index, kept].copy()
            kept_count = width - void_count
            out[row_index, :] = kept_values + [void] * void_count
            out_nr[row_index, :] = False
            out_nr[row_index, :kept_count] = kept_nr
            out_void[row_index, :] = False
            out_void[row_index, kept_count:] = True
    return out, out_nr, out_void


def _resolved_ids(
    requested: str | Sequence[object] | None,
    inherited: list[object] | None,
    nseq: int,
) -> list[object] | None:
    if requested is None:
        values = inherited
    elif isinstance(requested, str) and requested == "auto":
        values = [f"[{index}]" for index in range(1, nseq + 1)]
    elif isinstance(requested, str):
        values = [requested]
    else:
        values = list(requested)
    if values is None:
        return None
    if len(values) != nseq:
        raise ValueError("invalid 'row.names' length")
    if any(_cell_is_na(value) for value in values):
        raise ValueError("missing values in 'row.names' are not allowed")
    if pd.Index(values).has_duplicates:
        raise ValueError("duplicate 'row.names' are not allowed")
    return values


def _rename_states(
    raw: np.ndarray, excluded: np.ndarray, old_states: Sequence[str], new_states: Sequence[str]
) -> None:
    """Rename every state of ``raw`` once, in place, outside ``excluded`` cells.

    All cells are selected before any is renamed. Renaming one state at a time, as R
    does (seqdef.R:78-81), renames again a cell whose new name is a later old one and
    merges distinct states, against seqdef.Rd:42, which gives one label per alphabet
    state (2026-10-25 study, finding F-01; a DIVERGENCE row of
    docs/API_PARITY_DECISIONS.csv).
    """
    selected = [(raw == old_state) & ~excluded for old_state in old_states]
    for cells, new_state in zip(selected, new_states):
        raw[cells] = new_state


def seqdef(
    data: Union[
        STSList,
        pd.DataFrame,
        np.ndarray,
        Sequence[Sequence[Union[str, int, float, None]]],
    ],
    var: Sequence[int | str] | int | str | None = None,
    *,
    informat: str = "STS",
    stsep: str | None = None,
    alphabet: str | Sequence[str] | None = None,
    states: str | Sequence[str] | None = None,
    id: str | Sequence[object] | None = None,
    weights: Sequence[float] | None = None,
    start: int | None = 1,
    left: object = None,
    right: object = "DEL",
    gaps: object = None,
    missing: object = None,
    void: str = "%",
    nr: str = "*",
    cnames: str | Sequence[str] | None = None,
    xtstep: int = 1,
    tick_last: bool = False,
    cpal: str | Sequence[str] | None = None,
    missing_color: str = "darkgrey",
    labels: str | Sequence[str] | None = None,
    **seqformat_kwargs: Any,
) -> STSList:
    """Define categorical state trajectories using TraMineR-style preprocessing.

    When ``alphabet`` is omitted, the real states take the order of TraMineR's
    ``seqstatl``, as in R: by numeric value when every value in the data (a
    ``missing`` code included) is an R number, R's collation otherwise (ICU
    root collation, see docs/DEV_BUILD.md). ``states``, ``labels``, ``cpal``
    and the default palette are paired with the states in that order, so every
    state gets the name, label and colour R gives it. An explicit ``alphabet``
    keeps the caller's order. The collation is the ICU library that is loaded,
    as R's is its own, and ICU versions change over time: a few characters
    order differently between them, so for an order that matches R's whatever
    the versions, pass ``alphabet=`` explicitly. A state spelled as a non-canonical number
    (``"01"``, ``"0x10"``) is kept and ordered by its value, where R turns such
    cells into ``NA``.
    Internally observed nonresponse/void markers are appended to that real-state
    alphabet because :class:`STSList` uses one integer codebook. Integer ``var``
    selectors retain the established Python zero-based (and negative-index)
    convention; column names are portable between Python and R.

    On the STS route (text or array input), a cell or an ``alphabet`` entry
    spelled like ``nr`` or ``void`` is read as that marker: ``alphabet=["A",
    "*"]`` with ``nr="*"`` declares the nonresponse marker, it does not add a
    state. R lists such a spelling as a state, but its object then carries the
    text twice (a state level and the marker level) and ``seqlength(...,
    with.missing=FALSE)`` counts those cells as missing. The SPS and SPELL
    routes, and an :class:`STSList` input, tell markers apart by their codes
    and keep a literal ``"*"`` or ``"%"`` as a real state. To keep such a
    state on the STS route, give ``nr``/``void`` other spellings.

    Time labels (``xtlab``) come only from ``cnames``. Unlike R, which keeps
    the source column names (``y2001``, ...) or names matrix columns ``T1``,
    ``T2``, ..., the source column names are not copied: without ``cnames``,
    ``xtlab`` is ``None`` and plots label the axis from ``start``. Pass
    ``cnames=list(frame.columns)`` to keep them.
    """
    if informat not in {"STS", "SPS", "SPELL"}:
        raise ValueError("'informat' must be one of: STS, SPS, SPELL")
    left_state = _state_option(left)
    right_state = _state_option(right)
    gap_state = _state_option(gaps)
    if not isinstance(void, str) or void in {left_state, gap_state, right_state}:
        raise ValueError("'void' must be a character different from left, gaps, and right!")
    if not isinstance(nr, str) or nr == void:
        raise ValueError("'nr' and 'void' must be distinct character markers")
    if any(policy == nr for policy in (left_state, gap_state, right_state)):
        raise ValueError("missing replacement states must differ from the 'nr' marker")

    source, inherited_ids, columns, source_marker_masks = _extract_seqdef_table(data, var)
    source_nr_marker = data.nr if isinstance(data, STSList) else nr
    source_void_marker = data.void if isinstance(data, STSList) else void
    code_aware_real_states: set[str] | None = None
    if isinstance(data, STSList):
        marker_codes = {code for code in (data.nr_code, data.void_code) if code is not None}
        code_aware_real_states = {
            _r_state_text(state)
            for code, state in enumerate(data.alphabet)
            if code not in marker_codes
        }
    normalized_issues = (
        None
        if not isinstance(data, STSList) or data.issues is None
        else {category: list(indexes) for category, indexes in data.issues.items()}
    )
    if informat == "STS":
        if seqformat_kwargs:
            unknown_arguments = ", ".join(sorted(seqformat_kwargs))
            raise TypeError(f"unsupported STS seqdef arguments: {unknown_arguments}")
        separator = stsep
        if separator is not None and not isinstance(separator, str):
            raise ValueError("'stsep' must be a character")
        if separator is None and source.shape[1] == 1:
            observed = [value for value in source[:, 0] if not _cell_is_na(value)]
            if observed and all(isinstance(value, str) for value in observed):
                if any("-" in value for value in observed):
                    separator = "-"
                elif any(":" in value for value in observed):
                    separator = ":"
        if separator is not None:
            source, columns = _decompress_seqdef_sts(source, separator=separator)
            source_marker_masks = None
            code_aware_real_states = None
        raw = _render_seqdef_matrix(source)
    else:
        if isinstance(data, STSList):
            raise ValueError("SPS/SPELL input must be a data frame or matrix")
        missing_for_format = "_!NA!_" if _cell_is_na(missing) else _r_state_text(missing)
        selected_for_format = _selectors(var)
        source_width = (
            data.shape[1]
            if isinstance(data, (pd.DataFrame, np.ndarray)) and data.ndim == 2
            else cast(tuple[int, int], source.shape)[1]
        )
        format_var = (
            None
            if selected_for_format is None
            else [
                (
                    selector + 1
                    if isinstance(selector, int) and selector >= 0
                    else source_width + selector + 1 if isinstance(selector, int) else selector
                )
                for selector in selected_for_format
            ]
        )
        normalized = _seqformat(
            data,
            var=format_var,
            from_=informat,
            to="STS",
            stsep=stsep,
            missing=missing_for_format,
            **seqformat_kwargs,
        )
        assert isinstance(normalized, STSList)
        raw = np.empty(normalized.data.shape, dtype=object)
        special_codes = {
            code for code in (normalized.nr_code, normalized.void_code) if code is not None
        }
        for index, code in np.ndenumerate(normalized.data):
            raw[index] = None if int(code) in special_codes else normalized.alphabet[int(code)]
        # A zero-valued mask pair marks this route as code-aware: cells decoded
        # from special codes are already None, while literal marker spellings
        # that came through as real states must not be reclassified by text.
        source_marker_masks = (
            np.zeros(raw.shape, dtype=bool),
            np.zeros(raw.shape, dtype=bool),
        )
        code_aware_real_states = {
            _r_state_text(state)
            for code, state in enumerate(normalized.alphabet)
            if code not in special_codes
        }
        columns = None if normalized.xtlab is None else list(normalized.xtlab)
        inherited_ids = None if normalized.ids is None else list(normalized.ids)
        normalized_issues = (
            None
            if normalized.issues is None
            else {category: list(indexes) for category, indexes in normalized.issues.items()}
        )
        if informat == "SPS" and inherited_ids is None:
            inherited_ids = [f"[{index}]" for index in range(1, len(raw) + 1)]

    # R's seqstatl decides numeric or character order over every value in the
    # data -- a missing code and marker spellings included -- before seqdef
    # drops the missing code (seqdef.R:59-62).
    seqstatl_seen = [value for value in raw.ravel() if value is not None]
    missing_state = (
        None if informat in {"SPS", "SPELL"} or _cell_is_na(missing) else _r_state_text(missing)
    )
    missing_mask = np.zeros(raw.shape, dtype=bool)
    preencoded_nr = np.zeros(raw.shape, dtype=bool)
    preencoded_void = np.zeros(raw.shape, dtype=bool)
    code_aware_markers = source_marker_masks is not None
    for index, value in np.ndenumerate(raw):
        if source_marker_masks is not None and source_marker_masks[0][index]:
            if missing_state is not None and missing_state in {value, source_nr_marker}:
                missing_mask[index] = True
                raw[index] = None
                continue
            raw[index] = nr
            preencoded_nr[index] = True
            continue
        if source_marker_masks is not None and source_marker_masks[1][index]:
            if missing_state is not None and missing_state in {value, source_void_marker}:
                missing_mask[index] = True
                raw[index] = None
                continue
            raw[index] = void
            preencoded_void[index] = True
            continue
        is_missing = value is None or (missing_state is not None and value == missing_state)
        missing_mask[index] = is_missing
        if is_missing:
            raw[index] = None
        elif source_marker_masks is None:
            preencoded_nr[index] = value == nr
            preencoded_void[index] = value == void

    excluded = missing_mask | preencoded_nr | preencoded_void
    observed_states = _stable_states(raw, excluded)
    declared_alphabet: list[str] | None = None
    declared_real_positions: list[int] | None = None
    declared_nr = False
    declared_void = False
    if alphabet is not None:
        alphabet_values = [alphabet] if isinstance(alphabet, str) else alphabet
        declared_alphabet = [_r_state_text(value) for value in alphabet_values]
        duplicated = pd.Index(declared_alphabet)[pd.Index(declared_alphabet).duplicated()].tolist()
        if duplicated:
            raise ValueError(f"\n [!] duplicated element(s) in alphabet: {duplicated[0]}")
        # A marker spelling in a code-aware source is an ordinary declared state
        # when that spelling was observed at a non-marker code.  Actual marker
        # cells retain their separate masks and will receive a collision-free
        # spelling during materialization.
        declared_nr = nr in declared_alphabet and not (
            code_aware_markers
            and code_aware_real_states is not None
            and nr in code_aware_real_states
        )
        declared_void = void in declared_alphabet and not (
            code_aware_markers
            and code_aware_real_states is not None
            and void in code_aware_real_states
        )
        declared_real_positions = [
            position
            for position, state in enumerate(declared_alphabet)
            if not ((state == nr and declared_nr) or (state == void and declared_void))
        ]
        source_alphabet = [declared_alphabet[position] for position in declared_real_positions]
        unknown_states = [state for state in observed_states if state not in source_alphabet]
        if unknown_states:
            preview = ", ".join(unknown_states[:5])
            suffix = "" if len(unknown_states) <= 5 else f", ... (+{len(unknown_states) - 5} more)"
            raise ValueError(
                f"Data contains states not present in provided alphabet: {preview}{suffix}"
            )
    else:
        # Without an alphabet, R takes the seqstatl order of every value it
        # sees (seqdef.R:59-66; seqstatl.R:21-30) and pairs states=, labels=,
        # cpal= and its default palette with it (seqdef.R:78-82, 167-176,
        # 215-226): owner decision D2 of 2026-10-19.
        source_alphabet = r_seqstatl_order(observed_states, seen=seqstatl_seen)

    if states is None:
        real_alphabet = list(source_alphabet)
    else:
        state_values = [states] if isinstance(states, str) else states
        mapped_states = [_r_state_text(value) for value in state_values]
        if (
            declared_alphabet is not None
            and declared_real_positions is not None
            and len(mapped_states) == len(declared_alphabet)
        ):
            mapped_states = [mapped_states[position] for position in declared_real_positions]
        if len(mapped_states) != len(source_alphabet):
            raise ValueError(
                "'states' must have the same length as the alphabet "
                f"({len(source_alphabet)} states), got {len(mapped_states)}"
            )
        real_alphabet = mapped_states
        if len(set(real_alphabet)) != len(real_alphabet):
            raise ValueError("'states' must contain unique state codes")
        _rename_states(raw, excluded, source_alphabet, real_alphabet)

    if missing_state is not None:
        mapped_missing = (raw == missing_state) & ~preencoded_nr & ~preencoded_void
        missing_mask |= mapped_missing
        raw[mapped_missing] = None

    for policy in (left_state, gap_state, right_state):
        if (
            isinstance(policy, str)
            and policy not in ("DEL", "NEUTRAL")
            and policy not in real_alphabet
        ):
            real_alphabet.append(policy)
    raw, preencoded_nr, preencoded_void = _apply_missing_policies(
        raw,
        missing_mask,
        preencoded_nr,
        preencoded_void,
        left=left_state,
        gaps=gap_state,
        right=right_state,
        void=void,
    )
    final_nr = preencoded_nr.copy()
    final_void = preencoded_void.copy()
    for index, value in np.ndenumerate(raw):
        unresolved = value is None or (
            missing_state is not None
            and value == missing_state
            and not final_nr[index]
            and not final_void[index]
        )
        if unresolved:
            final_nr[index] = True
            final_void[index] = False

    used_nr = declared_nr or bool(final_nr.any())
    used_void = declared_void or bool(final_void.any())
    effective_nr = _collision_free_marker(nr, real_alphabet, "NR") if used_nr else nr
    occupied = list(real_alphabet)
    if used_nr:
        occupied.append(effective_nr)
    effective_void = _collision_free_marker(void, occupied, "VOID") if used_void else void
    raw[final_nr] = effective_nr
    raw[final_void] = effective_void
    internal_alphabet = list(real_alphabet)
    if used_nr:
        internal_alphabet.append(effective_nr)
    if used_void:
        internal_alphabet.append(effective_void)
    codebook = {state: index for index, state in enumerate(internal_alphabet)}
    encoded = np.empty(raw.shape, dtype=int)
    for index, value in np.ndenumerate(raw):
        state = str(value)
        if state not in codebook:
            raise ValueError(f"state {state!r} is not represented in the final alphabet")
        encoded[index] = codebook[state]

    if weights is None:
        w = None
    else:
        w = _validate_sequence_weights(weights, raw.shape[0])

    nr_label = effective_nr
    void_label = effective_void
    added_policy_count = len(real_alphabet) - len(source_alphabet)
    if labels is None:
        real_labels = list(real_alphabet)
    else:
        label_values = [labels] if isinstance(labels, str) else labels
        supplied_labels = [str(label) for label in label_values]
        if (
            declared_alphabet is not None
            and declared_real_positions is not None
            and len(supplied_labels) == len(declared_alphabet) + added_policy_count
        ):
            real_labels = [supplied_labels[position] for position in declared_real_positions]
            real_labels.extend(supplied_labels[len(declared_alphabet) :])
            if declared_nr:
                nr_label = supplied_labels[declared_alphabet.index(nr)]
            if declared_void:
                void_label = supplied_labels[declared_alphabet.index(void)]
        elif len(supplied_labels) == len(real_alphabet) + int(used_nr) + int(used_void):
            real_labels = supplied_labels[: len(real_alphabet)]
            marker_labels = iter(supplied_labels[len(real_alphabet) :])
            if used_nr:
                nr_label = next(marker_labels)
            if used_void:
                void_label = next(marker_labels)
        else:
            real_labels = supplied_labels
        if len(real_labels) != len(real_alphabet):
            raise ValueError("\n [!] number of labels must equal number of states in the alphabet")
    nr_color = missing_color
    void_color = "white"
    if cpal is None:
        real_palette = _default_seqdef_palette(len(real_alphabet))
    else:
        palette_values = [cpal] if isinstance(cpal, str) else cpal
        supplied_palette = [str(color) for color in palette_values]
        if (
            declared_alphabet is not None
            and declared_real_positions is not None
            and len(supplied_palette) == len(declared_alphabet) + added_policy_count
        ):
            real_palette = [supplied_palette[position] for position in declared_real_positions]
            real_palette.extend(supplied_palette[len(declared_alphabet) :])
            if declared_nr:
                nr_color = supplied_palette[declared_alphabet.index(nr)]
            if declared_void:
                void_color = supplied_palette[declared_alphabet.index(void)]
        elif len(supplied_palette) == len(real_alphabet) + int(used_nr) + int(used_void):
            real_palette = supplied_palette[: len(real_alphabet)]
            marker_colors = iter(supplied_palette[len(real_alphabet) :])
            if used_nr:
                nr_color = next(marker_colors)
            if used_void:
                void_color = next(marker_colors)
        else:
            real_palette = supplied_palette
        if len(real_palette) != len(real_alphabet):
            raise ValueError("\n [!] number of colors in 'cpal' must equal length of alphabet")
    internal_labels = list(real_labels)
    internal_palette = list(real_palette)
    if used_nr:
        internal_labels.append(nr_label)
        internal_palette.append(nr_color)
    if used_void:
        internal_labels.append(void_label)
        internal_palette.append(void_color)

    start_value = 1 if start is None else int(start)
    if cnames is not None:
        column_names = [cnames] if isinstance(cnames, str) else cnames
        xtlab = [str(name) for name in column_names]
        if len(xtlab) != raw.shape[1]:
            raise ValueError(f"cnames must have length {raw.shape[1]}")
    else:
        xtlab = None
    final_ids = _resolved_ids(id, inherited_ids, raw.shape[0])

    return STSList(
        data=np.asarray(encoded, dtype=int),
        alphabet=internal_alphabet,
        labels=internal_labels,
        cpal=internal_palette,
        xtlab=xtlab,
        xtstep=xtstep,
        tick_last=bool(tick_last),
        start=start_value,
        weights=w,
        nr=effective_nr,
        void=effective_void,
        nr_code=codebook.get(effective_nr) if used_nr else None,
        void_code=codebook.get(effective_void) if used_void else None,
        ids=final_ids,
        missing_color=missing_color,
        issues=normalized_issues,
    )


def _fce_argument_values(value: object) -> list[object]:
    if isinstance(value, (str, Integral)) and not isinstance(value, bool):
        return [value]
    try:
        return list(cast(Sequence[object], value))
    except TypeError:
        return [value]


def _fce_column(seqdata: pd.DataFrame, selector: object) -> object:
    if isinstance(selector, str):
        if selector not in seqdata.columns:
            raise ValueError("undefined columns selected")
        return selector
    if isinstance(selector, Integral) and not isinstance(selector, bool):
        position = int(selector) - 1
        if position < 0 or position >= seqdata.shape[1]:
            raise ValueError("undefined columns selected")
        return seqdata.columns[position]
    raise ValueError("undefined columns selected")


def FCE_to_TSE(
    seqdata: pd.DataFrame,
    *,
    id: int | str | Sequence[object] | None = None,
    cols: Sequence[int | str] | int | str = (),
    eventlist: Optional[Sequence[object]] = None,
    firstEvent: object | None = None,
    time_col: str = "time",
) -> pd.DataFrame:
    """Convert fixed-column events to an ``id``/time/event table.

    Numeric selectors are one-based, as in TraMineRextras 0.6.9. Events are
    sorted by time and then label *within each source row*; person blocks retain
    source order. Labels (and character times) sort in R's collation order,
    ``_x, a, b, B, C``, not by code points. ``firstEvent`` is always placed at
    time zero. Set
    ``time_col="timestamp"`` for compatibility with the former Python schema.

    Two unsafe upstream indexing quirks are deliberately corrected: eventless
    people contribute no rows (and an all-eventless input returns a true empty
    three-column table), and a generated ID for a one-row input remains ``1``.
    Event labels must also match ``cols`` one-for-one instead of silently
    producing missing labels or ignoring extras.
    """
    if not isinstance(seqdata, pd.DataFrame):
        raise ValueError("'seqdata' must be a data frame")
    selectors = _fce_argument_values(cols)
    if not selectors or len(seqdata) == 0:
        raise ValueError("argument is of length zero")
    event_columns = [_fce_column(seqdata, selector) for selector in selectors]

    if eventlist is None:
        event_names: list[object] = [_r_state_text(column) for column in event_columns]
    else:
        event_names = _fce_argument_values(eventlist)
        if len(event_names) != len(event_columns):
            raise ValueError("eventlist must have same length as cols")
        event_names = [
            pd.NA if _cell_is_na(value) else _r_state_text(value) for value in event_names
        ]

    if id is None:
        person_ids: list[object] = list(range(1, len(seqdata) + 1))
    else:
        id_values = _fce_argument_values(id)
        if not id_values:
            raise ValueError("argument is of length zero")
        if len(id_values) == 1:
            id_column = _fce_column(seqdata, id_values[0])
            person_ids = seqdata[id_column].tolist()
        else:
            person_ids = list(id_values[: len(seqdata)])
            if len(person_ids) < len(seqdata):
                person_ids.extend([pd.NA] * (len(seqdata) - len(person_ids)))

    observed_times = [
        value
        for column in event_columns
        for value in seqdata[column].tolist()
        if not _cell_is_na(value)
    ]
    character_times = any(isinstance(value, str) for value in observed_times)

    def normalized_time(value: object) -> object:
        if character_times:
            return _r_state_text(value)
        if isinstance(value, (bool, np.bool_)):
            return int(value)
        return value

    zero_time: object = "0" if character_times else 0
    first_event = (
        None
        if firstEvent is None
        else (pd.NA if _cell_is_na(firstEvent) else _r_state_text(firstEvent))
    )
    # FCE_to_TSE.R:40 order(timessel, eventsel): a character key sorts by R's
    # collation, not by code points (audit S2-07).
    event_rank = r_character_rank(
        [
            event
            for event in [*event_names, first_event]
            if event is not None and not _cell_is_na(event)
        ]
    )
    time_rank = (
        r_character_rank([zero_time, *(normalized_time(value) for value in observed_times)])
        if character_times
        else None
    )
    records: list[tuple[object, object, object]] = []
    for row_position, person_id in enumerate(person_ids):
        local: list[tuple[object, object, int]] = []
        if firstEvent is not None:
            local.append((zero_time, first_event, -1))
        for event_position, (column, event_name) in enumerate(
            zip(event_columns, event_names, strict=True)
        ):
            time = seqdata[column].iloc[row_position]
            if not _cell_is_na(time):
                local.append((normalized_time(time), event_name, event_position))

        def order_key(item: tuple[object, object, int]) -> tuple[object, bool, int, int]:
            time, event, position = item
            missing_event = _cell_is_na(event)
            return (
                time if time_rank is None else time_rank[str(time)],
                missing_event,
                -1 if missing_event else event_rank[str(event)],
                position,
            )

        local.sort(key=order_key)
        records.extend((person_id, time, event) for time, event, _position in local)

    return pd.DataFrame(records, columns=["id", time_col, "event"])


def _rows_by_person(id_values: Sequence[object]) -> dict[str, list[int]]:
    """Each person's row positions, in table order, found in one pass.

    ``TSE_to_STS`` scanned every row for each person, which made it quadratic
    (2026-10-25 study, finding F-02).
    """
    rows: dict[str, list[int]] = {}
    for position, value in enumerate(id_values):
        if not _cell_is_na(value):
            rows.setdefault(str(value), []).append(position)
    return rows


def TSE_to_STS(
    df: pd.DataFrame,
    *,
    id_col: int | str = "id",
    time_col: int | str = "time",
    state_col: int | str | None = "state",
    event_col: int | str = "event",
    stm: pd.DataFrame | None = None,
    tmin: int = 1,
    tmax: int | None = None,
    first_state: str = "None",
    **kwargs,
) -> STSList:
    """Convert a time-stamped event or state table to state sequences.

    A table with the selected ``state_col`` and no event-only options retains
    the original Python state-at-time pivot. Otherwise, the selected
    ``id_col``/``time_col``/``event_col`` columns follow
    ``TraMineRextras::TSE_to_STS``: events are ordered by time and then name
    (in R's collation order, ``_x, a, b, B, C``, not by code points) within
    each person, unknown explicit-matrix events are ignored, and the
    completed ``1..tmax`` grid is sliced to ``tmin..tmax``. Numeric column
    selectors are one-based in both modes. In both modes the rows keep the
    persons' order of first appearance, as R's event mode and R's route for state
    data (``reshape``, then ``seqdef``) do (TX-AUD-20261019-039); the state mode
    refuses a missing id, as R's ``seqdef`` refuses it as a row name.

    Event output is an :class:`STSList` re-expression of R's data frame. Its
    IDs are R-style character row names, its time labels are ``a<tmin>`` through
    ``a<tmax>``, and its state values are the transition-matrix strings. The NA
    cells an NA transition leaves are decided as ``seqdef`` decides R's NA
    cells: trailing ones become void, leading and internal ones nonresponse.
    In state mode, too, an unobserved period follows those ``seqdef`` rules,
    and a state spelled like a marker stays a real state, as in
    :func:`seqdef_from_long` (where the difference from R's reshape-then-seqdef
    route is described).

    Three upstream pre-release defects are deliberately corrected: automatic
    matrices honor ``first_state`` and retain their dimensions for a one-event
    alphabet, and non-contiguous rows for one ID are grouped exactly once in
    stable first-appearance order.

    The alphabet is R's ``seqdef`` order of the states (``seqstatl``: numbers
    by value, R's collation otherwise), and ``labels`` and ``cpal`` are read in
    that order, as :func:`seqdef` reads them (in R they are given to ``seqdef``
    after ``TSE_to_STS``): one entry per real state or per codebook entry, by
    default the states and ``seqdef``'s palette, the markers taking their own
    label and colour (TX-AUD-20261019-035, TX-REP-20261019-013). There is no ``alphabet=``: ``seqdef(x, alphabet=[...])``
    on the result sets another order and keeps every cell and marker.
    """
    validate_seq_kwargs(kwargs, func_name="TSE_to_STS")
    if not isinstance(df, pd.DataFrame):
        raise ValueError("'seqdata' must be a data frame")

    event_options = (
        stm is not None
        or tmax is not None
        or tmin != 1
        or first_state != "None"
        or event_col != "event"
        or state_col is None
    )
    state_column: object | None = None
    if state_col is not None:
        try:
            state_column = _fce_column(df, state_col)
        except ValueError:
            state_column = None
    if state_column is not None and not event_options:
        id_column = _fce_column(df, id_col)
        time_column = _fce_column(df, time_col)
        long = df[[id_column, time_column, state_column]].copy()
        long.columns = ["id", "time", "state"]
        return seqdef_from_long(
            long,
            id_col="id",
            time_col="time",
            state_col="state",
            **kwargs,
        )

    # R resolves selectors in event/id/time order before validating the horizon.
    event_column = _fce_column(df, event_col)
    id_column = _fce_column(df, id_col)
    time_column = _fce_column(df, time_col)
    if tmin < 1:
        raise ValueError(" [!] tmin must be a valid numeric column index! (cannot be < 1)")
    if tmax is None:
        raise ValueError(" [!] tmax required to determine the length of the sequence")
    if isinstance(tmin, bool) or int(tmin) != tmin:
        raise ValueError("tmin must be an integer")
    if isinstance(tmax, bool) or int(tmax) != tmax or tmax < 1:
        raise ValueError("tmax must be a positive integer")
    tmin_value = int(tmin)
    tmax_value = int(tmax)
    if tmin_value > tmax_value:
        raise ValueError("subscript out of bounds")

    event_values = [
        pd.NA if _cell_is_na(value) else _r_state_text(value) for value in df[event_column].tolist()
    ]
    id_values = [
        pd.NA if _cell_is_na(value) else _r_state_text(value) for value in df[id_column].tolist()
    ]
    try:
        time_values = pd.to_numeric(df[time_column], errors="raise").to_numpy(
            dtype=float,
            na_value=np.nan,
        )
    except (TypeError, ValueError):
        raise ValueError("NAs introduced by coercion") from None

    def stable_nonmissing(values: Sequence[object]) -> list[str]:
        result: list[str] = []
        seen: set[str] = set()
        for value in values:
            if _cell_is_na(value):
                continue
            text = str(value)
            if text not in seen:
                result.append(text)
                seen.add(text)
        return result

    person_ids = stable_nonmissing(id_values)
    event_names = stable_nonmissing(event_values)
    initial_state = _r_state_text(first_state)

    if stm is None:
        if not event_names:
            raise ValueError("subscript out of bounds")
        from .events._object import seqe2stm

        transition = seqe2stm(event_names, first_state=initial_state)
    else:
        if not isinstance(stm, pd.DataFrame):
            raise TypeError("stm must be a pandas DataFrame with labelled rows and columns")
        transition = stm.copy()
        transition.index = [_r_state_text(value) for value in transition.index]
        transition.columns = [_r_state_text(value) for value in transition.columns]
        if transition.index.has_duplicates or transition.columns.has_duplicates:
            raise ValueError("stm row and column labels must be unique")
        for index, value in np.ndenumerate(transition.to_numpy(dtype=object)):
            transition.iat[index] = pd.NA if _cell_is_na(value) else _r_state_text(value)

    states: np.ndarray = np.full((len(person_ids), tmax_value), initial_state, dtype=object)
    transition_events = set(str(value) for value in transition.columns)
    # TSE_TO_STS.R:35 order(mtime, mevent) sorts tied events by R's collation,
    # not by code points (audit S2-07).
    event_rank = r_character_rank(event_names)
    rows_of_person = _rows_by_person(id_values)
    for person_position, person_id in enumerate(person_ids):
        positions = rows_of_person[person_id]

        def event_order(position: int) -> tuple[bool, float, bool, int, int]:
            time = float(time_values[position])
            event = event_values[position]
            return (
                bool(np.isnan(time)),
                time if not np.isnan(time) else math.inf,
                _cell_is_na(event),
                -1 if _cell_is_na(event) else event_rank[str(event)],
                position,
            )

        positions.sort(key=event_order)
        previous_state = initial_state
        previous_time = 1
        for position in positions:
            event = event_values[position]
            if _cell_is_na(event) or str(event) not in transition_events:
                continue
            time = float(time_values[position])
            if np.isnan(time) or np.isneginf(time):
                raise ValueError("missing value where TRUE/FALSE needed")
            transition_time = (
                tmax_value if np.isposinf(time) else min(tmax_value, math.floor(time + 1.0))
            )
            if transition_time < 0:
                raise ValueError("subscript out of bounds")
            if transition_time >= previous_time:
                assignment_start = max(previous_time, 1) - 1
                states[person_position, assignment_start:transition_time] = previous_state
            elif transition_time == 0 and previous_time == 1:
                # R's 1:0 index writes the first cell (index zero is ignored).
                states[person_position, 0] = previous_state
            previous_time = transition_time
            if previous_state not in transition.index:
                raise ValueError("subscript out of bounds")
            next_state = transition.loc[previous_state, str(event)]
            if isinstance(next_state, pd.Series):
                raise ValueError("stm row and column labels must be unique")
            previous_state = pd.NA if _cell_is_na(next_state) else _r_state_text(next_state)
        final_start = max(previous_time, 1) - 1
        states[person_position, final_start:tmax_value] = previous_state

    output = states[:, tmin_value - 1 : tmax_value]
    derived = {"ids", "xtlab", "start"}.intersection(kwargs)
    if derived:
        raise TypeError(
            "TSE_to_STS event mode derives ids, xtlab, and start; "
            f"remove explicit values for {sorted(derived)}"
        )
    event_kwargs = dict(kwargs)
    observed_states = [
        _r_state_text(value) for value in output.ravel(order="C") if not _cell_is_na(value)
    ]
    requested_nr = event_kwargs.get("nr", "*")
    requested_void = event_kwargs.get("void", "%")
    effective_nr = (
        None
        if requested_nr is None
        else _collision_free_marker(_r_state_text(requested_nr), observed_states, "NR")
    )
    effective_void = (
        None
        if requested_void is None
        else _collision_free_marker(
            _r_state_text(requested_void),
            [*observed_states, *([] if effective_nr is None else [effective_nr])],
            "VOID",
        )
    )
    event_kwargs.update(
        ids=person_ids,
        xtlab=[f"a{position}" for position in range(tmin_value, tmax_value + 1)],
        start=tmin_value,
        nr=effective_nr,
        void=effective_void,
    )
    return build_sts_from_state_matrix(output, kwargs=event_kwargs)


def _resolve_cols(df: pd.DataFrame, cols: Sequence[str | int]) -> list[str]:
    """Column names for zero-based (or negative) positions and names, as R refuses others."""
    width = len(df.columns)
    resolved = []
    for column in cols:
        if isinstance(column, str):
            if column not in df.columns:
                raise ValueError("undefined columns selected")
            resolved.append(column)
            continue
        position = int(column)
        if not -width <= position < width:
            raise ValueError("undefined columns selected")
        resolved.append(df.columns[position])
    return resolved


def _hspell_horizontal_sts(
    df: pd.DataFrame,
    *,
    id_col: str | None,
    begin: Sequence[str | int],
    end: Sequence[str | int],
    status: Sequence[str | int] | None,
    fixed_status: str | Sequence[str] | None,
    pvar: str | None,
    overwrite: bool,
    overlap_priority: str,
    filler: str | None,
    tmin: float | None,
    tmax: float | None,
    endObs: str | float | Sequence[float] | None,
    time_grid: Optional[Sequence[int]],
    interval_closed: str,
    tmin_col: str | None,
    tmax_col: str | None,
    seq_kwargs: dict,
) -> STSList:
    """Port of ``TraMineRextras::HSPELL_to_STS`` (``HSPELL_to_STS.R``).

    One output row per input row, in input order. Each spell column is painted,
    in declared order, onto a matrix of NA: a cell ``t`` belongs to spell ``d``
    when ``begin <= t <= end``; a missing begin skips the spell and a missing
    end runs it to ``tmax``. ``overwrite=False`` paints only cells still NA.
    ``fixed_status[d]`` (NA beyond its length) replaces the status columns when
    given. Blank cells stay NA unless ``fillblanks``; then ``endObs``, the first
    unobserved period, blanks every cell from time ``endObs`` on. Times are read
    after subtracting ``pvar``; the grid is ``tmin:tmax``, inferred from every
    begin and end value, and columns are named ``a<t>``.
    The NA cells then go through ``seqdef``'s R defaults: leading and internal
    ones become missing, trailing ones void.

    Corrected R defect: R shifts the spell times by ``1 - tmin`` but ``endObs``
    only by ``- tmin`` (``HSPELL_to_STS.R:25,67-68``), so it also blanks the
    period ``endObs - 1``, before the end of observation; traminex blanks from
    ``endObs`` on (TX-REP-20261017-005).

    Deliberate divergence: R sets a missing end to the *unshifted* ``tmax``
    after shifting all times by ``1 - tmin`` (``HSPELL_to_STS.R:37``), so with
    ``tmin < 1`` an open spell stops ``1 - tmin`` periods before the end of the
    grid. Here it ends at ``tmax`` (``docs/API_PARITY.md``).

    ``overlap_priority``, ``interval_closed``, ``time_grid`` and the per-row
    ``tmin_col``/``tmax_col`` bounds are traminex extensions whose defaults
    leave R's behaviour unchanged.
    """
    if overlap_priority not in {"later", "earlier", "longer", "shorter"}:
        raise ValueError("overlap_priority must be later, earlier, longer, or shorter")
    begin_cols = _resolve_cols(df, begin)
    end_cols = _resolve_cols(df, end)
    if len(begin_cols) != len(end_cols):
        raise ValueError("begin and end must name the same number of spell columns")
    status_cols = None if status is None else _resolve_cols(df, status)
    if status_cols is not None and len(status_cols) != len(begin_cols):
        raise ValueError("status must name one column per spell")
    named = (id_col, pvar, endObs if isinstance(endObs, str) else None, tmin_col, tmax_col)
    _resolve_cols(df, [column for column in named if column is not None])
    n_rows = len(df)

    def numeric(column: str | int) -> np.ndarray:
        try:
            return pd.to_numeric(df[column]).to_numpy(dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"column {column!r} must hold numeric times (begin, end, endObs or pvar)"
            ) from error

    begins = np.column_stack([numeric(c) for c in begin_cols]) if n_rows else np.empty((0, 0))
    ends = np.column_stack([numeric(c) for c in end_cols]) if n_rows else np.empty((0, 0))
    origin = None
    if pvar is not None:
        origin = numeric(pvar)
        begins = begins - origin[:, None]
        ends = ends - origin[:, None]

    if time_grid is not None:
        grid = np.asarray(list(time_grid), dtype=float)
    else:
        observed = np.concatenate([begins.ravel(), ends.ravel()])
        observed = observed[~np.isnan(observed)]
        if (tmin is None or tmax is None) and observed.size == 0:
            raise ValueError("no spell begin or end value to infer tmin/tmax from")
        low = float(tmin) if tmin is not None else float(observed.min())
        high = float(tmax) if tmax is not None else float(observed.max())
        grid = np.arange(low, high + 1)
    open_end = float(grid.max()) if grid.size else 0.0

    if fixed_status is not None:
        fixed = [fixed_status] if isinstance(fixed_status, str) else list(fixed_status)
    status_values = None if status_cols is None else [df[column].tolist() for column in status_cols]
    cells = np.full((n_rows, grid.size), None, dtype=object)
    durations = np.zeros((n_rows, grid.size), dtype=float)
    for spell in range(len(begin_cols)):
        for row in range(n_rows):
            start, stop = begins[row, spell], ends[row, spell]
            if np.isnan(start):
                continue
            if np.isnan(stop):
                stop = open_end
            if fixed_status is not None:
                label = r_as_character(fixed[spell]) if spell < len(fixed) else None
            elif status_values is not None:
                label = r_as_character(status_values[spell][row])
            else:
                label = f"S{spell + 1}"
            lower_ok = grid > start if interval_closed == "right" else grid >= start
            upper_ok = grid < stop if interval_closed == "left" else grid <= stop
            duration = (
                (stop - (interval_closed == "left")) - (start + (interval_closed == "right")) + 1
            )
            for column in np.flatnonzero(lower_ok & upper_ok):
                unset = cells[row, column] is None
                if not overwrite:
                    paint = unset
                elif overlap_priority == "later":
                    paint = True
                elif overlap_priority == "earlier":
                    paint = unset
                elif overlap_priority == "longer":
                    paint = unset or duration >= durations[row, column]
                else:
                    paint = unset or duration < durations[row, column]
                if paint:
                    cells[row, column] = label
                    durations[row, column] = duration

    if filler is not None:
        cells[pd.isna(cells)] = filler
    # The per-row bounds mark each person's observation window: as R's endObs,
    # the periods outside it are unobserved, so they are blanked after the fill
    # and stay missing (before) or void (after) (TX-AUD-20261019-037).
    for bound_col, outside in ((tmin_col, np.less), (tmax_col, np.greater)):
        if bound_col is not None:
            bounds = numeric(bound_col)
            cells[outside(grid[None, :], bounds[:, None])] = None
    if endObs is not None:
        if isinstance(endObs, str):
            observation_end = numeric(endObs)
        elif np.ndim(endObs) == 0:
            observation_end = np.full(n_rows, float(endObs))  # type: ignore[arg-type]
        else:
            observation_end = np.asarray(endObs, dtype=float)
            if observation_end.shape != (n_rows,):
                raise ValueError("endObs must be a column name, one number, or one per row")
        if origin is not None:
            observation_end = observation_end - origin
        with np.errstate(invalid="ignore"):
            cells[grid[None, :] >= observation_end[:, None]] = None

    if id_col is None:
        ids = [str(index) for index in range(1, n_rows + 1)]
    else:
        ids = [str(r_as_character(value)) for value in df[id_col].tolist()]
        if len(set(ids)) != len(ids):
            raise ValueError(
                "id values must be unique: HSPELL_to_STS keeps one row per input row "
                "and R refuses duplicate row names"
            )
    nr = seq_kwargs.get("nr", "*")
    void = seq_kwargs.get("void", "%")
    # R's seqdef sorts the alphabet as seqstatl does: numbers by value when every
    # state is an R number, R's collation otherwise (TX-AUD-20261010-031).
    states = r_seqstatl_order(
        [
            str(cell)
            for cell in cells.ravel(order="F")
            if cell is not None and str(cell) not in {str(nr), str(void)}
        ]
    )
    seqdef_kwargs = {key: value for key, value in seq_kwargs.items() if key not in {"ids", "xtlab"}}
    return seqdef(
        pd.DataFrame(cells),
        alphabet=states or None,
        id=seq_kwargs.get("ids", ids),
        cnames=seq_kwargs.get("xtlab", [f"a{r_number_text(t)}" for t in grid]),
        **seqdef_kwargs,
    )


def HSPELL_to_STS(
    df: pd.DataFrame,
    *,
    id_col: str | None = None,
    # long (vertical spells) arguments
    start_col: str | None = None,
    end_col: str | None = None,
    state_col: str | None = None,
    # horizontal (HSPELL) arguments
    begin: Sequence[str | int] | None = None,
    end: Sequence[str | int] | None = None,
    status: Sequence[str | int] | None = None,
    fixed_status: str | Sequence[str] | None = None,
    pvar: str | None = None,
    overwrite: bool = True,
    overlap_priority: str = "later",  # "later" | "earlier" | "longer" | "shorter"
    fillblanks: str | None = None,
    tmin: int | None = None,
    tmax: int | None = None,
    endObs: str | float | Sequence[float] | None = None,
    time_grid: Optional[Sequence[int]] = None,
    fill: str | None = None,
    interval_closed: str = "both",  # "both" | "left" | "right"
    tmin_col: str | None = None,
    tmax_col: str | None = None,
    **kwargs,
) -> STSList:
    """Convert spell data to STS (supports both vertical and horizontal HSPELL inputs).

    - Vertical: provide start_col, end_col, state_col (optionally id_col) — one row per spell.
      This layout is a traminex extension: the spells are regrouped to one row
      per person (first-appearance order, the person's spells in input order)
      and converted as horizontal input, so both layouts give the same object.
      endObs, pvar and per-person bounds must hold one value per person.
    - Every column argument must name a column of ``df``; otherwise R's
      ``undefined columns selected`` is raised.
    - Horizontal: provide lists begin, end, status (and id_col) — one row per subject.
    - overwrite: if False, a spell paints only cells no earlier spell covered.
    - fillblanks: code for uncovered periods; as in R they stay missing when
      both ``fillblanks`` and ``fill`` are None.
    - tmin/tmax or time_grid define the common time grid; otherwise inferred from data.
    - endObs: end-of-observation time, the first unobserved period -- a column
      name, one number, or one value per row; cells from time ``endObs`` on are
      blanked (R also blanks the period before, a corrected R defect).
    - fixed_status: one status per spell column (R ``fixed.status``), used
      instead of the status columns.
    - tmin_col / tmax_col: columns giving each person's first and last observed
      period (a traminex extension). As ``endObs``, the periods outside that
      window are unobserved: they are blanked after ``fillblanks``, so they are
      never filled, and become missing (before) or void (after) by ``seqdef``'s
      defaults (TX-AUD-20261019-037).

    Horizontal input returns one row per input row in input order, with ids
    from ``id_col`` (which must be unique) and time labels ``a<t>``, exactly as
    ``TraMineRextras::HSPELL_to_STS`` followed by ``seqdef``; see
    ``_hspell_horizontal_sts`` for the full contract and the one deliberate
    divergence (an open spell with ``tmin < 1``).

    Notes
    -----
    Horizontal spell columns are applied in their declared order, matching
    ``TraMineRextras::HSPELL_to_STS``. ``overlap_priority="later"`` is the R
    overwrite policy. The ``"earlier"``, ``"longer"``, and ``"shorter"`` modes
    are deterministic Python extensions; priority changes overlap cells only.

    Internals are split into helper units for: horizontal spell parsing, overlap
    resolution, bounds/grid construction, and matrix rendering/encoding.
    """
    validate_seq_kwargs(kwargs, func_name="HSPELL_to_STS")
    if overlap_priority not in {"later", "earlier", "longer", "shorter"}:
        raise ValueError("overlap_priority must be later, earlier, longer, or shorter")
    if begin is not None or end is not None:
        if not (begin and end):
            raise ValueError("begin and end must both be provided for horizontal HSPELL")
        return _hspell_horizontal_sts(
            df,
            id_col=id_col,
            begin=begin,
            end=end,
            status=status,
            fixed_status=fixed_status,
            pvar=pvar,
            overwrite=overwrite,
            overlap_priority=overlap_priority,
            filler=fillblanks if fillblanks is not None else fill,
            tmin=tmin,
            tmax=tmax,
            endObs=endObs,
            time_grid=time_grid,
            interval_closed=interval_closed,
            tmin_col=tmin_col,
            tmax_col=tmax_col,
            seq_kwargs=kwargs,
        )
    # Vertical spells (a traminex extension): one row per spell. They are
    # reshaped to the horizontal layout, one row per person in first-appearance
    # order with that person's spells in input order, so both layouts give the
    # same object (TX-AUD-20261011-019).
    if not (start_col and end_col and state_col):
        raise ValueError("Provide start_col, end_col, state_col for vertical spells")
    if endObs is not None and not isinstance(endObs, str):
        raise ValueError("vertical spells take endObs as a column name")
    # The person-level columns travel under the argument names, which cannot
    # collide with the reshaped id, b<k>, e<k> and s<k> columns.
    person_cols = {
        argument: column
        for argument, column in (
            ("endObs", endObs),
            ("pvar", pvar),
            ("tmin_col", tmin_col),
            ("tmax_col", tmax_col),
        )
        if column is not None
    }
    _resolve_cols(
        df,
        [column for column in (id_col, start_col, end_col, state_col) if column is not None]
        + list(person_cols.values()),
    )
    wide, begin_cols, end_cols, status_cols, wide_id = _vertical_spells_to_horizontal(
        df,
        id_col=id_col,
        start_col=start_col,
        end_col=end_col,
        state_col=state_col,
        person_cols=person_cols,
    )
    return _hspell_horizontal_sts(
        wide,
        id_col=wide_id,
        begin=begin_cols,
        end=end_cols,
        status=status_cols,
        fixed_status=None,
        pvar="pvar" if "pvar" in person_cols else None,
        overwrite=overwrite,
        overlap_priority=overlap_priority,
        filler=fillblanks if fillblanks is not None else fill,
        tmin=tmin,
        tmax=tmax,
        endObs="endObs" if "endObs" in person_cols else None,
        time_grid=time_grid,
        interval_closed=interval_closed,
        tmin_col="tmin_col" if "tmin_col" in person_cols else None,
        tmax_col="tmax_col" if "tmax_col" in person_cols else None,
        seq_kwargs=kwargs,
    )


def _vertical_spells_to_horizontal(
    df: pd.DataFrame,
    *,
    id_col: str | None,
    start_col: str,
    end_col: str,
    state_col: str,
    person_cols: Mapping[str, str],
) -> tuple[pd.DataFrame, list[str], list[str], list[str], str | None]:
    """One row per person, with spell columns ``b<k>``/``e<k>``/``s<k>``.

    Persons are the distinct ids in first-appearance order (without an id
    column and without an ``id`` column, each row is a person); a person's
    spells keep the input order and shorter lists are padded with NaN begins,
    which the horizontal route skips. Person-level columns (endObs, pvar,
    per-person bounds) must hold one value per person; ``person_cols`` maps
    the name each gets in the result to the input column.
    """
    if id_col is None and "id" in df.columns:
        id_col = "id"
    if id_col is None:
        keys = [str(index) for index in range(1, len(df) + 1)]
        first_values: list[object] = list(range(1, len(df) + 1))
    else:
        keys = [str(r_as_character(value)) for value in df[id_col].tolist()]
        first_values = df[id_col].tolist()
    persons: dict[str, list[int]] = {}
    for row, key in enumerate(keys):
        persons.setdefault(key, []).append(row)
    width = max((len(rows) for rows in persons.values()), default=0)
    starts = df[start_col].tolist()
    ends = df[end_col].tolist()
    states = df[state_col].tolist()
    columns: dict[str, list[object]] = {"id": [first_values[rows[0]] for rows in persons.values()]}
    for k in range(width):
        columns[f"b{k + 1}"] = [
            starts[rows[k]] if k < len(rows) else np.nan for rows in persons.values()
        ]
        columns[f"e{k + 1}"] = [
            ends[rows[k]] if k < len(rows) else np.nan for rows in persons.values()
        ]
        columns[f"s{k + 1}"] = [
            states[rows[k]] if k < len(rows) else None for rows in persons.values()
        ]
    for name, column in person_cols.items():
        values = df[column].tolist()
        per_person = []
        for key, rows in persons.items():
            distinct = {str(r_as_character(values[row])) for row in rows}
            if len(distinct) > 1:
                raise ValueError(
                    f"column {column!r} must hold one value per person; id {key!r} has {sorted(distinct)}"
                )
            per_person.append(values[rows[0]])
        columns[name] = per_person
    wide = pd.DataFrame(columns)
    return (
        wide,
        [f"b{k + 1}" for k in range(width)],
        [f"e{k + 1}" for k in range(width)],
        [f"s{k + 1}" for k in range(width)],
        "id" if id_col is not None else None,
    )


# ---------------------------------------------------------------------------
# Faithful R seqformat STS->SPELL / STS->SRS hub (RP6b; closes T18)
#
# Provenance: ports of TraMineR 2.2-14 ``seqformat-STS_to_SPELL.R`` /
# ``seqformat-STS_to_SRS.R`` (+ ``seqshift.R``). traminex is GPL-3.0-or-later, so a
# GPL->GPL port of TraMineR's own source is clean (no BSD-3 header / NOTICE append).
# These are the spell/shifted-replication converters R's ``seqformat`` dispatches to;
# unlike the wide<->long convenience ``seqformat()`` they are R-faithful. The
# downstream ``seqsurv`` (STS->SPELL) and ``seqauto``/``seqcta`` (STS->SRS) build on
# them. Golden-pinned vs ``seqformat(from="STS", to="SPELL"/"SRS")``.
# ---------------------------------------------------------------------------


def _real_alphabet_codes(x: STSList) -> list[int]:
    """Alphabet codes of real states (excluding nr / void), in alphabet order."""
    return [i for i in range(len(x.alphabet)) if i != x.void_code and i != x.nr_code]


def _spell_pdata(
    pdata: pd.DataFrame,
    pvar: Sequence[int | str] | None,
    nseq: int,
) -> tuple[list[object], np.ndarray]:
    if pvar is None:
        raise ValueError("'pvar' required when 'pdata' is a data frame")
    selectors = list(pvar)
    if len(selectors) < 2:
        raise ValueError("'pvar' must select ID and birth columns")
    columns: list[object] = []
    for selector in selectors[:2]:
        if isinstance(selector, str):
            if selector not in pdata.columns:
                raise ValueError("undefined columns selected")
            columns.append(selector)
        elif isinstance(selector, Integral) and not isinstance(selector, bool):
            position = int(selector) - 1
            if position < 0 or position >= pdata.shape[1]:
                raise ValueError("undefined columns selected")
            columns.append(pdata.columns[position])
        else:
            raise ValueError("'pvar' must contain positive integers or strings")
    selected = pdata.loc[:, columns]
    if len(selected) != nseq:
        raise ValueError("'pdata' id column must contain one entry per sequence")
    births = selected.iloc[:, 1]
    if (
        pd.api.types.is_datetime64_any_dtype(births.dtype)
        or pd.api.types.is_object_dtype(births.dtype)
        or isinstance(births.dtype, pd.CategoricalDtype)
        or pd.api.types.is_string_dtype(births.dtype)
    ):
        raise ValueError(
            f" [!] Column {columns[1]} of pdata should contain integer (birth year) values"
        )
    try:
        offsets = pd.to_numeric(births, errors="raise").to_numpy(dtype=float) - 1.0
    except (TypeError, ValueError):
        raise ValueError(
            f" [!] Column {columns[1]} of pdata should contain integer (birth year) values"
        ) from None
    return list(selected.iloc[:, 0]), offsets


def _prepare_spell_source(x: STSList, right: object) -> STSList:
    # Both marker representations give the same spells (TX-REP-20261023-003): a
    # marker declared in the alphabet beside an external code is the marker.
    x = x._canonicalized()
    marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
    raw = np.empty(x.data.shape, dtype=object)
    for index, code in np.ndenumerate(x.data):
        numeric_code = int(code)
        raw[index] = None if numeric_code in marker_codes else x.alphabet[numeric_code]
    real_states = [x.alphabet[code] for code in _real_alphabet_codes(x)]
    # The decoded text is re-read by seqdef, which would take a real state
    # spelled like a marker ("*") for that marker: give the markers spellings
    # no real state uses (audit S2-09).
    nr = _collision_free_marker("*" if x.nr is None else str(x.nr), real_states, "NR")
    void = _collision_free_marker(
        "%" if x.void is None else str(x.void), [*real_states, nr], "VOID"
    )
    return seqdef(
        raw,
        alphabet=real_states,
        right=right,
        start=x.start,
        nr=nr,
        void=void,
    )


def _sts_to_spell_core(
    x: STSList,
    *,
    with_missing: bool,
    right: object,
    pdata: pd.DataFrame | None,
    pvar: Sequence[int | str] | None,
    ids: Sequence[object] | None,
    rerun: bool = True,
) -> pd.DataFrame:
    # seqformat re-runs seqdef with its ``right`` before extracting spells
    # (seqformat.R:228-248); R's internal STS_to_SPELL reads the object as it is.
    prepared = _prepare_spell_source(x, right) if rerun else x._canonicalized()
    if pdata is not None:
        person_ids, offsets = _spell_pdata(pdata, pvar, prepared.nseq)
    else:
        person_ids = list(ids) if ids is not None else list(range(1, prepared.nseq + 1))
        if len(person_ids) != prepared.nseq:
            raise ValueError("ids must contain one entry per sequence")
        offsets = np.zeros(prepared.nseq, dtype=float)

    real_codes = _real_alphabet_codes(prepared)
    statl = [prepared.alphabet[code] for code in real_codes]
    code_to_level = {code: level for level, code in enumerate(real_codes, start=1)}
    if with_missing:
        statl.append(prepared.nr)
        if prepared.nr_code is not None:
            code_to_level[prepared.nr_code] = len(statl)

    output_ids: list[object] = []
    begins: list[float] = []
    ends: list[float] = []
    state_codes: list[int] = []
    for row_index, row in enumerate(prepared.data):
        stop = len(row)
        if prepared.void_code is not None:
            while stop > 0 and int(row[stop - 1]) == prepared.void_code:
                stop -= 1
        position = 0
        while position < stop:
            code = int(row[position])
            end = position + 1
            while end < stop and int(row[end]) == code:
                end += 1
            keep = code != prepared.void_code and (with_missing or code != prepared.nr_code)
            if keep:
                output_ids.append(person_ids[row_index])
                begins.append(float(offsets[row_index] + position + 1))
                ends.append(float(offsets[row_index] + end))
                state_codes.append(code_to_level[code])
            position = end

    states = [statl[code - 1] for code in state_codes]
    out = pd.DataFrame({"id": output_ids, "begin": begins, "end": ends, "states": states})
    if all(float(value).is_integer() for value in begins + ends):
        out["begin"] = out["begin"].astype(int)
        out["end"] = out["end"].astype(int)
    out.attrs["statl"] = statl
    # A tuple, not an array: pandas compares attrs when it combines frames
    # (astype with a mapping, concat), and an array comparison is ambiguous.
    out.attrs["state_code"] = tuple(int(code) for code in state_codes)
    return out


def _STS_to_SPELL_hub(
    x: STSList,
    *,
    with_missing: bool = True,
    right: object = "DEL",
    pdata: pd.DataFrame | None = None,
    pvar: Sequence[int | str] | None = None,
    ids: Sequence[object] | None = None,
) -> pd.DataFrame:
    return _sts_to_spell_core(
        x,
        with_missing=with_missing,
        right=right,
        pdata=pdata,
        pvar=pvar,
        ids=ids,
    )


@_require_stslist
def STS_to_SPELL(
    x: STSList,
    *,
    with_missing: bool = True,
    right: object = None,
    pdata: pd.DataFrame | None = None,
    pvar: Sequence[int | str] | None = None,
) -> pd.DataFrame:
    """Convert STS sequences to SPELL form (TraMineR ``seqformat(from='STS', to='SPELL')``).

    Returns one row per spell — a maximal run of an identical state — with columns
    ``id`` (the object's ids, as R's row names; 1..n when it has none),
    ``begin`` / ``end`` (1-based positions of the run) and ``states``
    (the state label). Void cells are dropped and break runs. With ``with_missing``
    (default, as in R's hub) the ``nr`` code is kept as an ordinary spell state;
    with ``with_missing=False`` the ``nr`` code is mapped to void (so nr cells are
    dropped and runs do not merge across them). The markers are told apart by
    their codes, so a real state spelled like a marker (a literal ``"*"``)
    stays a real state.

    By default the object is read as it is, as R's internal ``STS_to_SPELL`` does,
    so trailing nonresponse cells (``seqdef(..., right=None)``) form a spell. A
    ``right`` policy (``"DEL"``, as ``seqformat(to="SPELL")`` uses) first re-runs
    ``seqdef`` with it on the decoded cells, as R's ``seqformat`` does.

    Provenance: faithful port of TraMineR 2.2-14 ``seqformat-STS_to_SPELL.R``.
    """
    return _sts_to_spell_core(
        x,
        with_missing=with_missing,
        right=right,
        pdata=pdata,
        pvar=pvar,
        ids=x.ids,
        rerun=right is not None,
    )


@_require_stslist
def STS_to_SRS(x: STSList, *, nrep: int = 0) -> pd.DataFrame:
    """Convert STS sequences to SRS (shifted replication) form.

    For every position ``idx`` in ``(1+nrep)..L`` and every sequence, emits a row
    ``(id, idx, T-(L-1), ..., T-1, T)`` where the state at ``idx`` sits in column
    ``T`` and the preceding states fill ``T-1, T-2, ...`` right-aligned, with leading
    ``NaN`` for the unavailable history (R ``seqshift``). Rows are ordered by
    ``(id, idx)``. State cells carry the state labels (``NaN`` in the leading region).

    Provenance: faithful port of TraMineR 2.2-14 ``seqformat-STS_to_SRS.R`` /
    ``seqshift.R``.
    """
    L = x.seqlen
    n = x.nseq
    if nrep == -1:
        raise ValueError("numbers of columns of arguments do not match")
    if nrep < -1:
        raise ValueError("only 0's may be mixed with negative subscripts")
    if 1 + nrep >= L:
        raise ValueError("invalid 'times' argument")
    labels = np.array(x.alphabet, dtype=object)
    label_mat = np.full(x.data.shape, None, dtype=object)
    # R unifies nr + void into NA before shifting (seqformat: mseqdata[... %in% missing] <- NA).
    miss_cell = np.zeros_like(x.data, dtype=bool)
    if x.void_code is not None:
        miss_cell |= x.data == x.void_code
    if x.nr_code is not None:
        miss_cell |= x.data == x.nr_code
    real_cells = ~miss_cell
    label_mat[real_cells] = labels[x.data[real_cells]]  # (n, L) symbols
    label_mat[miss_cell] = np.nan

    tcols = [f"T-{k}" for k in range(L - 1, 0, -1)] + ["T"]
    rows: list[list] = []
    for idx in range(1 + nrep, L + 1):
        for i in range(n):
            window = [np.nan] * (L - idx) + list(label_mat[i, 0:idx])
            rows.append([i + 1, idx] + window)
    rows.sort(key=lambda r: (r[0], r[1]))
    df = pd.DataFrame(rows, columns=["id", "idx"] + tcols)
    df["id"] = df["id"].astype(int)
    df["idx"] = df["idx"].astype(int)
    return df
