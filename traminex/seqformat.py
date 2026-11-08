from __future__ import annotations

import re
import warnings
from collections.abc import Callable, Mapping
from numbers import Real
from typing import Any, Sequence, cast

import numpy as np
import pandas as pd

from ._categorical_order import r_character_rank, r_seqstatl_order
from ._deprecated import _caller_stacklevel
from ._rformat import r_as_character
from ._sequence_builders import (
    build_sts_from_state_matrix,
    seqdef_labels_and_palette,
    validate_seq_kwargs,
)
from .sequences import STSList, _require_stslist


class _OmittedId:
    """Sentinel whose repr preserves the public ``id=1`` signature spelling."""

    def __repr__(self) -> str:
        return "1"


_OMITTED_ID = _OmittedId()


def seqformat(
    data: Any,
    *,
    from_: str | None = None,
    to: str | None = None,
    id_col: str | None = None,
    time_col: str | None = None,
    state_col: str | None = None,
    var: Sequence[int | str] | int | str | None = None,
    compress: bool = False,
    nrep: int | None = None,
    tevent: Any = None,
    stsep: str | None = None,
    covar: Sequence[int | str] | None = None,
    sps_in: Mapping[str, str] | None = None,
    sps_out: Mapping[str, str] | None = None,
    id: int | str | Sequence[object] | None = _OMITTED_ID,  # type: ignore[assignment]
    begin: int | str = 2,
    end: int | str = 3,
    status: int | str = 4,
    process: bool = True,
    pdata: pd.DataFrame | str | None = None,
    pvar: Sequence[int | str] | None = None,
    limit: int | float = 100,
    overwrite: bool = True,
    fillblanks: str | None = None,
    tmin: int | None = None,
    tmax: int | None = None,
    missing: str | Sequence[str] = "*",
    with_missing: bool = True,
    right: str | None = "DEL",
    compressed: bool | None = None,
    nr: str | None = None,
) -> pd.DataFrame | STSList | list[str]:
    """Convert sequence representations while retaining legacy wide/long routes.

    The lower-case ``wide``/``long`` API remains the original Python pivot.
    Upper-case ``STS``/``SPS``/``SPELL`` inputs follow TraMineR 2.2-14. R's
    raw ``to="STS"`` matrix is re-expressed as an :class:`STSList`; trailing
    missing cells become structural void while leading/interior missing cells
    use the missing-state code. If literal states already use ``*`` or ``%``,
    collision-free ``<NR>``/``<VOID>`` markers preserve both meanings. The
    public ``right`` option is intentionally ignored for this target, exactly
    as in R.

    STS and DSS outputs use :class:`STSList`, SPS uses one spell string per
    sequence, and SRS/TSE/SPELL preserve the R data-frame column schemas.
    """
    if from_ in ("wide", "long") or to in ("wide", "long"):
        return _seqformat_wide_long(
            data,
            from_=from_,
            to=to,
            id_col=id_col,
            time_col=time_col,
            state_col=state_col,
        )

    if compressed is not None:
        if compress is not False:
            raise ValueError(
                "In seqformat() : compress and compressed cannot be specified together. "
                "compressed is deprecated, use compress instead."
            )
        warnings.warn(
            "In seqformat() : compressed is deprecated, use compress instead.",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        compress = bool(compressed)
    if nr is not None:
        if missing != "*":
            raise ValueError(
                "In seqformat() : missing and nr cannot be specified together. "
                "nr is deprecated, use missing instead."
            )
        warnings.warn(
            "In seqformat() : nr is deprecated, use missing instead.",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        missing = nr

    id_was_omitted = id is _OMITTED_ID
    effective_id: int | str | Sequence[object] | None = 1 if id_was_omitted else id

    is_stslist = isinstance(data, STSList)
    if from_ is None:
        if is_stslist:
            from_ = "STS"
        else:
            raise ValueError("'from' must be specified")
    if to is None:
        raise ValueError("'to' must be specified")

    input_formats = ("STS", "SPS", "SPELL")
    output_formats = ("STS", "DSS", "SPS", "SRS", "TSE", "SPELL")
    if from_ not in input_formats:
        raise ValueError("'from' must be one of: STS, SPS, SPELL")
    if to not in output_formats:
        raise ValueError("'to' must be one of: STS, DSS, SPS, SRS, TSE, SPELL")

    missing_tokens = _validate_missing(missing)
    _validate_stsep(stsep, from_)

    # R validates these options whenever SPELL appears on either side. The
    # SPELL-input normalizer handles its own branch; output-only routes need
    # the same warning/error boundary before dispatch.
    if to == "SPELL" and from_ != "SPELL":
        if pvar is not None:
            _selector_list(pvar, name="pvar")
            if not isinstance(pdata, pd.DataFrame):
                warnings.warn(
                    "'pvar' ignored because 'pdata' is not a data frame",
                    UserWarning,
                    stacklevel=2,
                )
        elif isinstance(pdata, pd.DataFrame):
            raise ValueError("'pvar' required when 'pdata' is a data frame")

    if covar is not None and to != "SRS":
        warnings.warn(
            'covar ignored when to is not "SRS"',
            UserWarning,
            stacklevel=2,
        )
        covar = None
    if to == "TSE" and tevent is None:
        raise ValueError("'tevent' must be specified")
    if compress is True and to not in {"STS", "DSS", "SPS"}:
        raise ValueError(f"{to} is not a compressible format")

    if to == "TSE":
        raw, normalized_ids = _normalize_tse_input(
            data,
            from_=from_,
            var=var,
            stsep=stsep,
            missing=missing_tokens,
            sps_in=sps_in,
            id=effective_id,
            begin=begin,
            end=end,
            status=status,
            process=process,
            pdata=pdata,
            pvar=pvar,
            limit=limit,
            overwrite=overwrite,
            fillblanks=fillblanks,
            tmin=tmin,
            tmax=tmax,
        )
        output_ids = _resolve_tse_ids(
            data,
            from_=from_,
            id=effective_id,
            id_was_omitted=id_was_omitted,
            normalized_ids=normalized_ids,
            nseq=raw.shape[0],
        )
        assert tevent is not None
        return _sts_to_tse_events(raw, output_ids, tevent)

    normalized = _normalize_input_to_sts(
        data,
        from_=from_,
        var=var,
        stsep=stsep,
        missing=missing_tokens,
        sps_in=sps_in,
        id=effective_id,
        begin=begin,
        end=end,
        status=status,
        process=process,
        pdata=pdata,
        pvar=pvar,
        limit=limit,
        overwrite=overwrite,
        fillblanks=fillblanks,
        tmin=tmin,
        tmax=tmax,
    )
    if to == "STS":
        return _compress_sts(normalized) if compress is True else normalized
    if to == "DSS":
        converted_dss = _sts_to_dss(normalized)
        return _compress_dss(converted_dss) if compress is True else converted_dss
    if to == "SPS":
        return _sts_to_sps(normalized, sps_out=sps_out, right=right)
    if to == "SRS":
        from .format import STS_to_SRS

        converted_srs = STS_to_SRS(normalized, nrep=0 if nrep is None else nrep)
        return _append_srs_covariates(converted_srs, data, covar)

    from .format import _STS_to_SPELL_hub

    spell_source = normalized
    if from_ == "SPS" and spell_source.ids is None:
        spell_source = spell_source.copy()
        spell_source.ids = [f"[{index}]" for index in range(1, spell_source.nseq + 1)]
    return _STS_to_SPELL_hub(
        spell_source,
        with_missing=with_missing,
        right=right,
        pdata=pdata if isinstance(pdata, pd.DataFrame) else None,
        pvar=pvar,
        ids=spell_source.ids,
    )


def _seqformat_wide_long(
    data: pd.DataFrame,
    *,
    from_: str | None,
    to: str | None,
    id_col: str | None,
    time_col: str | None,
    state_col: str | None,
) -> pd.DataFrame:
    """The original Python wide/long pivot, kept byte-for-byte compatible."""
    if from_ == to:
        return data.copy()
    if from_ == "wide" and to == "long":
        df = data.copy()
        if id_col is None:
            df = df.copy()
            df.insert(0, "id", np.arange(1, len(df) + 1, dtype=int))
            id_col = "id"
        long = df.melt(id_vars=[id_col], var_name="time", value_name="state")
        # Normalize time to integer order 1..L preserving original column order
        col_order = [c for c in df.columns if c != id_col]
        mapping = {col: i + 1 for i, col in enumerate(col_order)}
        long["time"] = long["time"].map(mapping)
        long = long.rename(columns={id_col: "id"})
        return long[["id", "time", "state"]].sort_values(["id", "time"]).reset_index(drop=True)
    if from_ == "long" and to == "wide":
        if not (id_col and time_col and state_col):
            raise ValueError("id_col, time_col, and state_col are required for long->wide")
        df = data[[id_col, time_col, state_col]].copy()
        wide = df.pivot(index=id_col, columns=time_col, values=state_col)
        # Sort columns by the time column's own order: numbers and dates in
        # value order, plain strings as text, and an ordered categorical by
        # its categories (sorted() ignored those; audit S2-06).
        wide = wide.sort_index(axis=1)
        wide.columns.name = None
        wide = wide.reset_index()
        return wide
    raise ValueError("Unsupported format conversion")


def _normalize_input_to_sts(
    data: Any,
    *,
    from_: str,
    var: Sequence[int | str] | int | str | None,
    stsep: str | None,
    missing: tuple[str, ...],
    sps_in: Mapping[str, str] | None,
    id: int | str | Sequence[object] | None,
    begin: int | str,
    end: int | str,
    status: int | str,
    process: bool,
    pdata: pd.DataFrame | str | None,
    pvar: Sequence[int | str] | None,
    limit: int | float,
    overwrite: bool,
    fillblanks: str | None,
    tmin: int | None,
    tmax: int | None,
) -> STSList:
    if from_ == "STS":
        return _normalize_sts_input(data, var=var, stsep=stsep, missing=missing)
    if from_ == "SPS":
        return _normalize_sps_input(
            data,
            var=var,
            stsep=stsep,
            missing=missing,
            sps_in=sps_in,
        )
    return _normalize_spell_input(
        data,
        var=var,
        id=id,
        begin=begin,
        end=end,
        status=status,
        process=process,
        pdata=pdata,
        pvar=pvar,
        limit=limit,
        overwrite=overwrite,
        fillblanks=fillblanks,
        tmin=tmin,
        tmax=tmax,
        missing=missing,
    )


def _real_states(x: STSList) -> list[str]:
    marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
    return [state for code, state in enumerate(x.alphabet) if code not in marker_codes]


def _decoded_sts_matrix(x: STSList) -> np.ndarray:
    marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
    out = np.empty(x.data.shape, dtype=object)
    for index, code in np.ndenumerate(x.data):
        numeric_code = int(code)
        out[index] = None if numeric_code in marker_codes else x.alphabet[numeric_code]
    return out


def _joined_nonmissing_rows(x: STSList) -> np.ndarray:
    marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
    out = np.empty((x.nseq, 1), dtype=object)
    for row_index, row in enumerate(x.data):
        values = [x.alphabet[int(code)] for code in row if int(code) not in marker_codes]
        out[row_index, 0] = "-".join(values)
    return out


def _compress_sts(x: STSList) -> STSList:
    return _materialize_sts(
        _joined_nonmissing_rows(x),
        ids=None if x.ids is None else list(x.ids),
        columns=["Sequence"],
        source=x,
        start=x.start,
    )


def _sts_to_dss(x: STSList) -> STSList:
    raw = np.full(x.data.shape, None, dtype=object)
    for row_index, row in enumerate(x.data):
        stop = len(row)
        if x.void_code is not None:
            while stop > 0 and int(row[stop - 1]) == x.void_code:
                stop -= 1
        collapsed: list[int] = []
        for value in row[:stop]:
            code = int(value)
            if not collapsed or code != collapsed[-1]:
                collapsed.append(code)
        for position, code in enumerate(collapsed):
            if code in {x.nr_code, x.void_code}:
                raw[row_index, position] = None
            else:
                raw[row_index, position] = x.alphabet[code]

    return _materialize_sts(
        raw,
        ids=[f"[{index}]" for index in range(1, x.nseq + 1)],
        columns=[f"[{index}]" for index in range(1, x.seqlen + 1)],
        nr="*",
        void="%",
        real_alphabet=_real_states(x),
        source=x,
        start=1,
    )


def _compress_dss(x: STSList) -> STSList:
    out = np.empty((x.nseq, 1), dtype=object)
    for row_index, row in enumerate(x.data):
        values = [
            x._state_for_code(int(code))
            for code in row
            if x.void_code is None or int(code) != x.void_code
        ]
        out[row_index, 0] = "-".join(values)
    return _materialize_sts(
        out,
        ids=None if x.ids is None else list(x.ids),
        columns=["Sequence"],
        source=x,
        start=1,
    )


def _right_replacement(right: str | None) -> str | None:
    if right is None or right == "DEL":
        return right
    if right == "NEUTRAL":
        return "#"
    return _r_state_text(right)


def _sts_to_sps(
    x: STSList,
    *,
    sps_out: Mapping[str, str] | None,
    right: str | None,
) -> list[str]:
    spec: Mapping[str, object] = {"xfix": "()", "sdsep": ","} if sps_out is None else sps_out
    if not isinstance(spec, Mapping):
        raise ValueError("'sps_out' must be a mapping")
    xfix_value = spec.get("xfix", "")
    separator_value = spec.get("sdsep", "")
    xfix = "" if xfix_value is None else str(xfix_value)
    duration_separator = "" if separator_value is None else str(separator_value)
    prefix = xfix[:1]
    suffix = xfix[1:2]
    replacement = _right_replacement(right)
    special_codes = {code for code in (x.nr_code, x.void_code) if code is not None}

    output: list[str] = []
    for row in x.data:
        real_positions = [
            position for position, code in enumerate(row) if int(code) not in special_codes
        ]
        last_real = max(real_positions, default=-1)
        states: list[str] = []
        for position, value in enumerate(row):
            code = int(value)
            if code not in special_codes:
                states.append(x.alphabet[code])
            elif position <= last_real:
                states.append("*")
            elif replacement == "DEL":
                continue
            elif replacement is None:
                states.append("*")
            else:
                states.append(replacement)

        tokens: list[str] = []
        position = 0
        while position < len(states):
            state = states[position]
            end = position + 1
            while end < len(states) and states[end] == state:
                end += 1
            tokens.append(f"{prefix}{state}{duration_separator}{end - position}{suffix}")
            position = end
        output.append("-".join(tokens))
    return output


def _append_srs_covariates(
    converted: pd.DataFrame,
    data: Any,
    covar: Sequence[int | str] | None,
) -> pd.DataFrame:
    if covar is None:
        return converted

    if isinstance(data, pd.DataFrame):
        selected = _select_dataframe_columns(data, covar, name="covar")
    elif isinstance(data, np.ndarray):
        matrix = np.asarray(data, dtype=object)
        if matrix.ndim == 1:
            matrix = matrix.reshape(-1, 1)
        if matrix.ndim != 2:
            raise ValueError("undefined columns selected")
        selectors = _selector_list(covar, name="covar")
        assert selectors is not None
        if not isinstance(selectors[0], int):
            raise ValueError("undefined columns selected")
        positions = [int(selector) - 1 for selector in selectors]
        if any(position < 0 or position >= matrix.shape[1] for position in positions):
            raise ValueError("undefined columns selected")
        names = (
            ["subset.data....covar."]
            if len(positions) == 1
            else [f"X{index}" for index in range(1, len(positions) + 1)]
        )
        selected = pd.DataFrame(matrix[:, positions], columns=names)
    elif isinstance(data, STSList):
        positions = _stslist_column_positions(data, covar)
        names = (
            [f"T{position + 1}" for position in positions]
            if data.xtlab is None
            else [str(data.xtlab[position]) for position in positions]
        )
        decoded = np.empty((data.nseq, len(positions)), dtype=object)
        for index, code in np.ndenumerate(data.data[:, positions]):
            decoded[index] = data._state_for_code(int(code))
        categories = list(
            dict.fromkeys(data._state_for_code(code) for code in data._codebook_codes())
        )
        selected = pd.DataFrame(
            {
                name: pd.Categorical(decoded[:, index], categories=categories)
                for index, name in enumerate(names)
            }
        )
    else:
        raise ValueError("undefined columns selected")

    row_positions = converted["id"].to_numpy(dtype=int) - 1
    if row_positions.size and len(selected) <= int(row_positions.max()):
        raise ValueError("undefined columns selected")
    selected_rows = selected.iloc[row_positions].reset_index(drop=True)

    used = {str(column) for column in converted.columns}
    repaired: list[str] = []
    for column in selected_rows.columns:
        base = str(column)
        candidate = base
        suffix = 1
        while candidate in used:
            candidate = f"{base}.{suffix}"
            suffix += 1
        repaired.append(candidate)
        used.add(candidate)
    selected_rows.columns = repaired
    return pd.concat([converted.reset_index(drop=True), selected_rows], axis=1)


def _raw_sts_for_tse(
    data: Any,
    *,
    var: Sequence[int | str] | int | str | None,
    stsep: str | None,
    missing: tuple[str, ...],
) -> tuple[np.ndarray, list[object] | None]:
    if isinstance(data, STSList):
        positions = _stslist_column_positions(data, var)
        warnings.warn(
            f"'missing' set as \"c('{data.nr}','{data.void}')\", the 'nr' and 'void' "
            "code from the 'data' state sequence object",
            UserWarning,
            stacklevel=3,
        )
        raw = np.empty((data.nseq, len(positions)), dtype=object)
        for row_index in range(data.nseq):
            for output_position, source_position in enumerate(positions):
                raw[row_index, output_position] = data._state_for_code(
                    int(data.data[row_index, source_position])
                )
        return raw, None if data.ids is None else list(data.ids)

    matrix, ids, _columns = _table_to_matrix(data, var=var)
    compressed = matrix.shape[1] == 1
    if compressed:
        matrix, _columns = _decompress_one_column(matrix, stsep=stsep)
        return _consume_missing_tokens(matrix, missing), ids
    raw = np.empty(matrix.shape, dtype=object)
    for index, value in np.ndenumerate(matrix):
        raw[index] = None if _cell_is_na(value) else _r_state_text(value)
    return raw, ids


def _normalize_tse_input(
    data: Any,
    *,
    from_: str,
    var: Sequence[int | str] | int | str | None,
    stsep: str | None,
    missing: tuple[str, ...],
    sps_in: Mapping[str, str] | None,
    id: int | str | Sequence[object] | None,
    begin: int | str,
    end: int | str,
    status: int | str,
    process: bool,
    pdata: pd.DataFrame | str | None,
    pvar: Sequence[int | str] | None,
    limit: int | float,
    overwrite: bool,
    fillblanks: str | None,
    tmin: int | None,
    tmax: int | None,
) -> tuple[np.ndarray, list[object] | None]:
    if from_ == "STS":
        return _raw_sts_for_tse(data, var=var, stsep=stsep, missing=missing)
    normalized = _normalize_input_to_sts(
        data,
        from_=from_,
        var=var,
        stsep=stsep,
        missing=missing,
        sps_in=sps_in,
        id=id,
        begin=begin,
        end=end,
        status=status,
        process=process,
        pdata=pdata,
        pvar=pvar,
        limit=limit,
        overwrite=overwrite,
        fillblanks=fillblanks,
        tmin=tmin,
        tmax=tmax,
    )
    return _decoded_sts_matrix(normalized), normalized.ids


def _selector_values(data: Any, selector: int | str) -> list[object]:
    if isinstance(data, pd.DataFrame):
        column = _resolve_column(data, selector, name="id")
        return list(pd.unique(data[column]))
    if isinstance(data, np.ndarray):
        if not isinstance(selector, int) or isinstance(selector, bool) or selector < 1:
            raise ValueError("'id' invalid type or unknown value")
        matrix = np.asarray(data, dtype=object)
        if matrix.ndim != 2 or selector > matrix.shape[1]:
            raise ValueError("undefined columns selected")
        return list(pd.unique(matrix[:, selector - 1]))
    if isinstance(data, STSList):
        if not isinstance(selector, int) or isinstance(selector, bool) or selector < 1:
            raise ValueError("'id' invalid type or unknown value")
        if selector > data.seqlen:
            raise ValueError("undefined columns selected")
        return list(
            pd.unique([data._state_for_code(int(code)) for code in data.data[:, selector - 1]])
        )
    raise ValueError("'id' invalid type or unknown value")


def _resolve_tse_ids(
    data: Any,
    *,
    from_: str,
    id: int | str | Sequence[object] | None,
    id_was_omitted: bool,
    normalized_ids: list[object] | None,
    nseq: int,
) -> list[object]:
    if from_ == "SPELL":
        return list(normalized_ids) if normalized_ids is not None else list(range(1, nseq + 1))
    if id_was_omitted:
        warnings.warn(
            "'id' set to NULL as it is not specified (backward compatibility with TraMineR 1.8)",
            UserWarning,
            stacklevel=3,
        )
        warnings.warn(
            "replacing original IDs in the output by the sequence indexes",
            UserWarning,
            stacklevel=3,
        )
        return list(range(1, nseq + 1))
    if id is None:
        warnings.warn(
            "replacing original IDs in the output by the sequence indexes",
            UserWarning,
            stacklevel=3,
        )
        return list(range(1, nseq + 1))
    if isinstance(id, Sequence) and not isinstance(id, (str, bytes, bytearray)):
        values = list(id)
        if not values:
            raise ValueError("'id' invalid type or unknown value")
        if len(values) == 1:
            only = values[0]
            if not isinstance(only, (int, str)) or isinstance(only, bool):
                raise ValueError("'id' invalid type or unknown value")
            values = _selector_values(data, only)
        elif pd.Index(values).has_duplicates:
            raise ValueError("'id' must contain unique IDs")
    elif isinstance(id, (int, str)) and not isinstance(id, bool):
        values = _selector_values(data, id)
    else:
        raise ValueError("'id' invalid type or unknown value")
    if len(values) != nseq:
        raise ValueError(f"'id' must contain {nseq} IDs")
    return values


def _tevent_values(tevent: pd.DataFrame, source: str, target: str) -> list[str] | None:
    if source not in tevent.index or target not in tevent.columns:
        return None
    value = tevent.loc[source, target]
    if _cell_is_na(value):
        return []
    if isinstance(value, str):
        return _r_strsplit(value, ",")
    return [_r_state_text(value)]


def _sts_to_tse_events(
    raw: np.ndarray,
    ids: Sequence[object],
    tevent: Any,
) -> pd.DataFrame:
    if not isinstance(tevent, pd.DataFrame):
        raise ValueError("'tevent' must be a data frame with state row and column names")
    rows: list[dict[str, object]] = []
    for row_index, sequence in enumerate(raw):
        if len(sequence) == 0:
            continue
        first = sequence[0]
        if first is not None:
            entrance = _tevent_values(tevent, str(first), str(first))
            if entrance:
                rows.extend({"id": ids[row_index], "time": 0, "event": event} for event in entrance)
        for position in range(len(sequence) - 1):
            source = sequence[position]
            target = sequence[position + 1]
            if source is None or target is None or source == target or source == "" or target == "":
                continue
            events = _tevent_values(tevent, str(source), str(target))
            if events is None:
                raise ValueError("argument is of length zero")
            rows.extend(
                {"id": ids[row_index], "time": position + 1, "event": event} for event in events
            )
    return pd.DataFrame(rows, columns=["id", "time", "event"])


def _validate_missing(missing: str | Sequence[str]) -> tuple[str, ...]:
    if isinstance(missing, str):
        return (missing,)
    if isinstance(missing, Sequence) and not isinstance(missing, (bytes, bytearray)):
        values = tuple(missing)
        if values and all(isinstance(value, str) for value in values):
            return values
    raise ValueError("'missing' invalid type or unknown value")


def _validate_stsep(stsep: str | None, from_: str) -> None:
    if from_ in ("STS", "SPS") and stsep is not None:
        if not isinstance(stsep, str) or len(stsep) != 1:
            raise ValueError("'stsep' must be a character")


def _selector_list(
    selector: Sequence[int | str] | int | str | None,
    *,
    name: str,
) -> list[int | str] | None:
    if selector is None:
        return None
    if isinstance(selector, (str, int)) and not isinstance(selector, bool):
        values: list[int | str] = [selector]
    elif isinstance(selector, Sequence) and not isinstance(selector, (str, bytes, bytearray)):
        values = list(selector)
    else:
        raise ValueError(f"'{name}' must be a vector of positive integers or strings")
    all_int = bool(values) and all(
        isinstance(value, int) and not isinstance(value, bool) and value > 0 for value in values
    )
    all_str = bool(values) and all(isinstance(value, str) for value in values)
    if not (all_int or all_str):
        raise ValueError(f"'{name}' must be a vector of positive integers or strings")
    return values


def _select_dataframe_columns(
    data: pd.DataFrame,
    selector: Sequence[int | str] | int | str | None,
    *,
    name: str = "var",
) -> pd.DataFrame:
    values = _selector_list(selector, name=name)
    if values is None:
        return data.copy()
    if isinstance(values[0], int):
        positions = [int(value) - 1 for value in values]
        if any(position < 0 or position >= data.shape[1] for position in positions):
            raise ValueError("undefined columns selected")
        return data.iloc[:, positions].copy()
    missing_columns = [value for value in values if value not in data.columns]
    if missing_columns:
        raise ValueError("undefined columns selected")
    return data.loc[:, values].copy()


def _select_array_columns(
    data: np.ndarray,
    selector: Sequence[int | str] | int | str | None,
) -> np.ndarray:
    values = _selector_list(selector, name="var")
    if values is None:
        return data.copy()
    if not isinstance(values[0], int):
        raise ValueError("'var' must be a vector of positive integers or strings")
    positions = [int(value) - 1 for value in values]
    if any(position < 0 or position >= data.shape[1] for position in positions):
        raise ValueError("undefined columns selected")
    return data[:, positions].copy()


def _has_default_index(data: pd.DataFrame | pd.Series) -> bool:
    return isinstance(data.index, pd.RangeIndex) and data.index.equals(pd.RangeIndex(len(data)))


def _dataframe_to_r_matrix(data: pd.DataFrame) -> np.ndarray:
    """Apply R's homogeneous logical/integer/double matrix promotion."""
    dtypes = list(data.dtypes)
    numeric = all(
        pd.api.types.is_numeric_dtype(dtype) or pd.api.types.is_bool_dtype(dtype)
        for dtype in dtypes
    )
    # pandas 3 may return a read-only view under Copy-on-Write.  This helper
    # normalizes cells in place, so request an owned writable array explicitly.
    matrix = data.to_numpy(dtype=object, copy=True)
    if not numeric:
        for column, dtype in enumerate(dtypes):
            if isinstance(dtype, pd.CategoricalDtype):
                for row in range(matrix.shape[0]):
                    if not _cell_is_na(matrix[row, column]):
                        matrix[row, column] = _r_state_text(matrix[row, column])
        return matrix

    convert: Callable[[Any], object]
    if any(pd.api.types.is_complex_dtype(dtype) for dtype in dtypes):
        convert = complex
    elif any(pd.api.types.is_float_dtype(dtype) for dtype in dtypes):
        convert = float
    elif any(pd.api.types.is_integer_dtype(dtype) for dtype in dtypes):
        convert = int
    else:
        convert = bool

    for index, value in np.ndenumerate(matrix):
        if not _cell_is_na(value):
            matrix[index] = convert(value)
    return matrix


def _table_to_matrix(
    data: Any,
    *,
    var: Sequence[int | str] | int | str | None,
    allow_empty_columns: bool = False,
) -> tuple[np.ndarray, list[object] | None, list[str]]:
    ids: list[object] | None = None
    if isinstance(data, pd.DataFrame):
        selected = _select_dataframe_columns(data, var)
        if selected.shape[1] == 0 and not allow_empty_columns:
            raise ValueError("'data' must contain at least one column after 'var' filtering")
        if not _has_default_index(selected):
            ids = list(selected.index)
        matrix = _dataframe_to_r_matrix(selected)
        return matrix, ids, [str(column) for column in selected.columns]
    if isinstance(data, pd.Series):
        if var is not None:
            raise ValueError("'var' must be a vector of positive integers or strings")
        if not _has_default_index(data):
            ids = list(data.index)
        name = "[1]" if data.name is None else str(data.name)
        return data.to_numpy(dtype=object).reshape(-1, 1), ids, [name]
    if isinstance(data, np.ndarray):
        matrix = np.asarray(data, dtype=object)
        if matrix.ndim == 1:
            matrix = matrix.reshape(-1, 1)
        if matrix.ndim != 2:
            raise ValueError("'data' must be a state sequence object, a data frame, or a matrix")
        matrix = _select_array_columns(matrix, var)
        if matrix.shape[1] == 0 and not allow_empty_columns:
            raise ValueError("'data' must contain at least one column after 'var' filtering")
        return matrix, None, [f"[{i}]" for i in range(1, matrix.shape[1] + 1)]
    if isinstance(data, Sequence) and not isinstance(data, (str, bytes, bytearray)):
        values = list(data)
        if all(isinstance(value, str) for value in values):
            if var is not None:
                raise ValueError("'var' must be a vector of positive integers or strings")
            return np.asarray(values, dtype=object).reshape(-1, 1), None, ["[1]"]
    raise ValueError("'data' must be a state sequence object, a data frame, or a matrix")


def _cell_is_na(value: object) -> bool:
    if value is None:
        return True
    try:
        result = pd.isna(value)
    except (TypeError, ValueError):
        return False
    return bool(result) if isinstance(result, (bool, np.bool_)) else False


def _r_state_text(value: object) -> str:
    """Render a scalar as R does when a state matrix becomes character."""
    if isinstance(value, (bool, np.bool_)):
        return "TRUE" if bool(value) else "FALSE"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if np.isposinf(numeric):
            return "Inf"
        if np.isneginf(numeric):
            return "-Inf"
        if numeric == 0:
            return "0"
        fixed = format(numeric, ".15g")
        scientific = np.format_float_scientific(
            numeric,
            precision=14,
            unique=True,
            trim="-",
            exp_digits=2,
        )
        return scientific if len(scientific) < len(fixed) else fixed
    return str(value)


def _r_strsplit(value: str, separator: str) -> list[str]:
    """Approximate base R ``strsplit`` for seqformat's scalar regex routes."""
    if separator == "":
        return list(value)
    try:
        parts = re.split(separator, value)
    except re.error as exc:
        raise ValueError(str(exc)) from exc
    # Base R discards exactly one empty component at the right boundary.
    if parts and parts[-1] == "":
        parts.pop()
    return parts


def _decompress_one_column(
    matrix: np.ndarray,
    *,
    stsep: str | None,
) -> tuple[np.ndarray, list[str]]:
    source_values = list(matrix[:, 0])
    observed_source = [value for value in source_values if not _cell_is_na(value)]
    separator = stsep
    if separator is None:
        if not observed_source:
            raise ValueError("missing value where TRUE/FALSE needed")
        all_strings = all(isinstance(value, str) for value in observed_source)
        if all_strings and any("-" in value for value in observed_source):
            separator = "-"
        elif all_strings and any(":" in value for value in observed_source):
            separator = ":"
        elif len(observed_source) != len(source_values):
            # seqformat evaluates max(nchar(seqdata) > 1) only when separator
            # detection failed; an NA makes that condition indeterminate.
            raise ValueError("missing value where TRUE/FALSE needed")
        elif not all_strings:
            rendered = [_r_state_text(value) for value in observed_source]
            if any(len(value) > 1 for value in rendered):
                raise ValueError("'stsep' must be specified as it is neither '-' nor ':'")
            raise ValueError("non-character argument")
        elif any(len(value) > 1 for value in observed_source):
            raise ValueError("'stsep' must be specified as it is neither '-' nor ':'")

    if any(not isinstance(value, str) for value in observed_source):
        raise ValueError("non-character argument")

    values = [None if _cell_is_na(value) else str(value) for value in source_values]

    rows: list[list[object]] = []
    for value in values:
        if value is None:
            rows.append([None])
        elif value == "":
            rows.append([])
        elif separator is None:
            rows.append([value])
        elif separator == "":
            rows.append(list(value))
        else:
            rows.append(cast(list[object], _r_strsplit(value, separator)))
    width = max((len(row) for row in rows), default=0)
    if width == 0:
        raise ValueError("length of 'dimnames' [2] not equal to array extent")
    out: np.ndarray = np.full((len(rows), width), None, dtype=object)
    for i, row in enumerate(rows):
        out[i, : len(row)] = row
    return out, [f"[{i}]" for i in range(1, width + 1)]


def _consume_missing_tokens(matrix: np.ndarray, missing: tuple[str, ...]) -> np.ndarray:
    out = np.empty(matrix.shape, dtype=object)
    for index, value in np.ndenumerate(matrix):
        text = None if _cell_is_na(value) else _r_state_text(value)
        if text is None or text in missing:
            out[index] = None
        else:
            out[index] = text
    return out


def _stable_real_states(matrix: np.ndarray) -> list[str]:
    states: list[str] = []
    seen: set[str] = set()
    for value in matrix.ravel(order="C"):
        if _cell_is_na(value):
            continue
        state = _r_state_text(value)
        if state not in seen:
            states.append(state)
            seen.add(state)
    return states


def _collision_free_marker(preferred: str, occupied: Sequence[str], role: str) -> str:
    if preferred not in occupied:
        return preferred
    candidate = f"<{role}>"
    suffix = 1
    while candidate in occupied:
        candidate = f"<{role}:{suffix}>"
        suffix += 1
    return candidate


def _materialize_sts(
    matrix: np.ndarray,
    *,
    ids: list[object] | None,
    columns: list[str],
    nr: str = "*",
    void: str = "%",
    real_alphabet: Sequence[str] | None = None,
    source: STSList | None = None,
    start: int = 1,
    issues: dict[str, list[int]] | None = None,
) -> STSList:
    raw = np.asarray(matrix, dtype=object)
    if raw.ndim != 2:
        raise ValueError("normalized STS input must be two-dimensional")

    # States the input does not declare take R's seqdef order (seqstatl):
    # owner decision D2 of 2026-10-19. A categorical SPELL status keeps its
    # declared categories, as R's factor status keeps its levels. An STSList
    # input keeps its alphabet's order, a traminex choice: R returns a plain
    # character matrix there, and seqdef of either result re-sorts the states.
    observed_states = r_seqstatl_order(_stable_real_states(raw))
    alphabet = [] if real_alphabet is None else [str(state) for state in real_alphabet]
    for state in observed_states:
        if state not in alphabet:
            alphabet.append(state)

    nr_used = False
    void_used = False
    last_observed_by_row: list[int] = []
    for i in range(raw.shape[0]):
        observed = [j for j in range(raw.shape[1]) if not _cell_is_na(raw[i, j])]
        last_observed = max(observed) if observed else -1
        last_observed_by_row.append(last_observed)
        for j in range(raw.shape[1]):
            if _cell_is_na(raw[i, j]) and j <= last_observed:
                nr_used = True
            elif _cell_is_na(raw[i, j]):
                void_used = True

    if nr_used:
        nr = _collision_free_marker(nr, alphabet, "NR")
        alphabet.append(nr)
    if void_used:
        void = _collision_free_marker(void, alphabet, "VOID")
        alphabet.append(void)

    classified = raw.copy()
    for i, last_observed in enumerate(last_observed_by_row):
        for j in range(raw.shape[1]):
            if not _cell_is_na(raw[i, j]):
                classified[i, j] = _r_state_text(raw[i, j])
            elif j <= last_observed:
                classified[i, j] = nr
            else:
                classified[i, j] = void

    code = {state: i for i, state in enumerate(alphabet)}
    encoded = np.empty(classified.shape, dtype=int)
    for index, value in np.ndenumerate(classified):
        encoded[index] = code[_r_state_text(value)]

    if source is None:
        # From a table, the labels and colours seqdef gives the same cells
        # (TX-AUD-20261019-035): R's seqformat returns a table that goes to seqdef.
        labels, cpal = seqdef_labels_and_palette(
            alphabet, nr=nr if nr_used else None, void=void if void_used else None
        )
    else:
        # Key the metadata repair by codebook code, not by alphabet text, so a
        # declared external marker resolves through _state_for_code exactly
        # like the data itself (cycle-8 audit F-35). Canonicalizing folds the
        # markers into the alphabet and extends labels/cpal with the marker
        # text / missing_color, which is what the fallbacks below hand out.
        labels = None
        cpal = None
        canonical_source = source._canonicalized()
        if canonical_source.labels is not None:
            label_map = dict(zip(canonical_source.alphabet, canonical_source.labels))
            labels = [label_map.get(state, state) for state in alphabet]
        if canonical_source.cpal is not None:
            color_map = dict(zip(canonical_source.alphabet, canonical_source.cpal))
            cpal = [
                color_map.get(
                    state,
                    "white" if state == void else canonical_source.missing_color,
                )
                for state in alphabet
            ]

    return STSList(
        data=encoded,
        alphabet=alphabet,
        labels=labels,
        cpal=cpal,
        xtlab=list(columns),
        xtstep=1 if source is None else source.xtstep,
        tick_last=False if source is None else source.tick_last,
        nr=nr,
        void=void,
        nr_code=code.get(nr) if nr_used else None,
        void_code=code.get(void) if void_used else None,
        weights=(
            None
            if source is None or source.weights is None
            else np.asarray(source.weights, dtype=float).copy()
        ),
        start=start,
        ids=None if ids is None else list(ids),
        missing_color="darkgrey" if source is None else source.missing_color,
        issues=issues if source is None or issues is not None else source.issues,
    )


def _stslist_column_positions(
    x: STSList,
    var: Sequence[int | str] | int | str | None,
) -> list[int]:
    values = _selector_list(var, name="var")
    if values is None:
        return list(range(x.seqlen))
    if isinstance(values[0], int):
        positions = [int(value) - 1 for value in values]
    else:
        if x.xtlab is None:
            raise ValueError("undefined columns selected")
        lookup = {str(label): i for i, label in enumerate(x.xtlab)}
        if any(str(value) not in lookup for value in values):
            raise ValueError("undefined columns selected")
        positions = [lookup[str(value)] for value in values]
    if any(position < 0 or position >= x.seqlen for position in positions):
        raise ValueError("undefined columns selected")
    return positions


def _normalize_sts_input(
    data: Any,
    *,
    var: Sequence[int | str] | int | str | None,
    stsep: str | None,
    missing: tuple[str, ...],
) -> STSList:
    if isinstance(data, STSList):
        # Both marker representations give the same output (TX-REP-20261023-003):
        # a marker declared in the alphabet beside an external code is the marker.
        data = data._canonicalized()
        positions = _stslist_column_positions(data, var)
        warnings.warn(
            f"'missing' set as \"c('{data.nr}','{data.void}')\", the 'nr' and 'void' "
            "code from the 'data' state sequence object",
            UserWarning,
            stacklevel=2,
        )
        raw = np.empty((data.nseq, len(positions)), dtype=object)
        special = {code for code in (data.nr_code, data.void_code) if code is not None}
        for i in range(data.nseq):
            for j, position in enumerate(positions):
                value = int(data.data[i, position])
                raw[i, j] = None if value in special else data.alphabet[value]
        real_alphabet = [
            state
            for i, state in enumerate(data.alphabet)
            if i not in (data.nr_code, data.void_code)
        ]
        if data.xtlab is None:
            columns = [f"T{position + 1}" for position in positions]
        else:
            columns = [str(data.xtlab[position]) for position in positions]
        return _materialize_sts(
            raw,
            ids=None if data.ids is None else list(data.ids),
            columns=columns,
            nr=data.nr,
            void=data.void,
            real_alphabet=real_alphabet,
            source=data,
            start=data.start,
        )

    matrix, ids, columns = _table_to_matrix(data, var=var)
    if matrix.shape[1] == 1:
        matrix, columns = _decompress_one_column(matrix, stsep=stsep)
    raw = _consume_missing_tokens(matrix, missing)
    return _materialize_sts(raw, ids=ids, columns=columns)


def _normalize_sps_input(
    data: Any,
    *,
    var: Sequence[int | str] | int | str | None,
    stsep: str | None,
    missing: tuple[str, ...],
    sps_in: Mapping[str, str] | None,
) -> STSList:
    table_types = (pd.DataFrame, pd.Series, np.ndarray)
    is_python_sequence = isinstance(data, Sequence) and not isinstance(
        data,
        (str, bytes, bytearray),
    )
    if isinstance(data, STSList):
        raise ValueError("'data' must be a data frame or a matrix")
    if is_python_sequence and not isinstance(data, table_types):
        if not all(isinstance(value, str) for value in data):
            raise ValueError("'data' SPS sequence must be one-dimensional and contain only strings")
    elif not isinstance(data, table_types):
        raise ValueError("'data' must be a data frame or a matrix")
    matrix, _source_ids, _columns = _table_to_matrix(data, var=var)
    if matrix.shape[1] == 1:
        matrix, _columns = _decompress_one_column(matrix, stsep=stsep)
    if not any(not _cell_is_na(value) for value in matrix.ravel(order="C")):
        raise ValueError("length of 'dimnames' [2] not equal to array extent")

    if sps_in is None:
        spec = {"xfix": "()", "sdsep": ","}
    else:
        if not isinstance(sps_in, Mapping):
            raise ValueError("'sps_in' must be a mapping")
        spec = dict(sps_in)
        if "xfix" not in spec:
            raise ValueError("invalid regular expression '[]', reason 'Missing ']''")
        if "sdsep" not in spec:
            raise ValueError("missing value where TRUE/FALSE needed")
    xfix = spec["xfix"]
    sdsep = spec["sdsep"]
    if not isinstance(xfix, str) or not isinstance(sdsep, str):
        raise ValueError("'sps_in' values must be strings")

    rows: list[list[object]] = []
    for source_row in matrix:
        trans: list[str] = []
        for cell in source_row:
            if _cell_is_na(cell):
                continue
            token = str(cell)
            if xfix:
                try:
                    token = re.sub(f"[{xfix}]", "", token)
                except re.error as exc:
                    raise ValueError(str(exc)) from exc
            parts = _r_strsplit(token, sdsep)
            if len(parts) < 2:
                raise ValueError("missing value where TRUE/FALSE needed")
            state = None if parts[0] in missing or parts[0] == "NA" else parts[0]
            try:
                duration_value = float(parts[1])
                if not np.isfinite(duration_value):
                    raise ValueError
                duration = int(duration_value)
            except (TypeError, ValueError):
                raise ValueError("missing value where TRUE/FALSE needed") from None
            trans.extend(["NA" if state is None else state] * max(duration, 1))
        expanded = _r_strsplit("-".join(trans), "-")
        rows.append([None if part == "NA" else part for part in expanded])

    width = max((len(row) for row in rows), default=0)
    if width == 0:
        raise ValueError("length of 'dimnames' [2] not equal to array extent")
    raw: np.ndarray = np.full((len(rows), width), None, dtype=object)
    for i, row in enumerate(rows):
        raw[i, : len(row)] = row
    columns = [f"[{i}]" for i in range(1, width + 1)]
    return _materialize_sts(raw, ids=None, columns=columns)


def _resolve_column(data: pd.DataFrame, selector: object, *, name: str) -> object:
    if isinstance(selector, str):
        if selector not in data.columns:
            raise ValueError("undefined columns selected")
        return selector
    if isinstance(selector, int) and not isinstance(selector, bool) and selector > 0:
        position = selector - 1
        if position >= data.shape[1]:
            raise ValueError("undefined columns selected")
        return data.columns[position]
    raise ValueError(f"'{name}' must be a positive integer or a string")


def _numeric_time_values(data: pd.Series) -> np.ndarray:
    if (
        pd.api.types.is_datetime64_any_dtype(data.dtype)
        or pd.api.types.is_object_dtype(data.dtype)
        or isinstance(data.dtype, pd.CategoricalDtype)
        or pd.api.types.is_string_dtype(data.dtype)
    ):
        raise ValueError(" [!] 'begin' and 'end' columns of seqdata should contain integer values")
    try:
        return pd.to_numeric(data, errors="raise").to_numpy(dtype=float, na_value=np.nan)
    except (TypeError, ValueError):
        raise ValueError(
            " [!] 'begin' and 'end' columns of seqdata should contain integer values"
        ) from None


def _positive_numeric(values: np.ndarray) -> np.ndarray:
    return values[np.isfinite(values) & (values > 0)]


def _resolve_pdata(
    pdata: pd.DataFrame,
    pvar: Sequence[int | str],
) -> tuple[object, object, dict[object, list[float | None]]]:
    selected = _select_dataframe_columns(pdata, pvar, name="pvar")
    if selected.shape[1] < 2:
        raise ValueError("'pvar' must select ID and birth columns")
    id_column, birth_column = selected.columns[:2]
    birth_series = selected.iloc[:, 1]
    if (
        pd.api.types.is_datetime64_any_dtype(birth_series.dtype)
        or pd.api.types.is_object_dtype(birth_series.dtype)
        or isinstance(birth_series.dtype, pd.CategoricalDtype)
        or pd.api.types.is_string_dtype(birth_series.dtype)
    ):
        raise ValueError(
            f" [!] Column '{birth_column}' of pdata should contain integer (birth year) values"
        )
    try:
        births = pd.to_numeric(birth_series, errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError):
        raise ValueError(
            f" [!] Column '{birth_column}' of pdata should contain integer (birth year) values"
        ) from None
    mapping: dict[object, list[float | None]] = {}
    for i in range(len(selected)):
        mapping.setdefault(selected.iloc[i, 0], []).append(
            float(births[i]) if np.isfinite(births[i]) else None
        )
    return id_column, birth_column, mapping


def _spell_position(
    value: float,
    *,
    mode: str,
    birth: float | None,
    axis_start: int,
) -> int:
    if mode == "age2age":
        return int(value)
    if mode == "age2year":
        assert birth is not None
        return int((birth - axis_start) + value + 1)
    if mode == "year2year":
        return int(value - axis_start + 1)
    assert birth is not None
    return int(value - birth + 1)


def _normalize_spell_input(
    data: Any,
    *,
    var: Sequence[int | str] | int | str | None,
    id: int | str | Sequence[object] | None,
    begin: int | str,
    end: int | str,
    status: int | str,
    process: bool,
    pdata: pd.DataFrame | str | None,
    pvar: Sequence[int | str] | None,
    limit: int | float,
    overwrite: bool,
    fillblanks: str | None,
    tmin: int | None,
    tmax: int | None,
    missing: tuple[str, ...],
) -> STSList:
    if not isinstance(data, pd.DataFrame) or data.shape[1] < 3:
        raise ValueError("'data' must be a data frame with at least three columns")
    if not isinstance(limit, Real) or isinstance(limit, bool):
        raise ValueError("'limit' must be a number")
    if fillblanks is not None and (not isinstance(fillblanks, str) or len(fillblanks) != 1):
        raise ValueError("'fillblanks' must be a character")
    if isinstance(pdata, pd.DataFrame) and pvar is None:
        raise ValueError("'pvar' required when 'pdata' is a data frame")
    if pvar is not None and not isinstance(pdata, pd.DataFrame):
        warnings.warn(
            "'pvar' ignored because 'pdata' is not a data frame",
            UserWarning,
            stacklevel=2,
        )
    if process and isinstance(pdata, str) and pdata != "auto":
        raise ValueError("incorrect number of dimensions")

    seqdata = _select_dataframe_columns(data, var)
    if seqdata.shape[1] == 0:
        raise ValueError("'data' must contain at least one column after 'var' filtering")
    id_column = _resolve_column(seqdata, id, name="id")
    begin_column = _resolve_column(seqdata, begin, name="begin")
    end_column = _resolve_column(seqdata, end, name="end")
    status_column = _resolve_column(seqdata, status, name="status")

    status_series = seqdata[status_column]
    real_alphabet: list[str] | None = None
    if isinstance(status_series.dtype, pd.CategoricalDtype):
        real_alphabet = [_r_state_text(value) for value in status_series.cat.categories]
        if fillblanks is not None and fillblanks not in real_alphabet:
            real_alphabet.append(fillblanks)

    begins = _numeric_time_values(seqdata[begin_column])
    ends = _numeric_time_values(seqdata[end_column])
    finite_times = np.concatenate([begins[np.isfinite(begins)], ends[np.isfinite(ends)]])
    if np.any(finite_times != np.trunc(finite_times)):
        raise ValueError(
            " [!] Found non-integer values in colunms referred to by begin and/or end argument!"
        )
    if np.any(begins[np.isfinite(begins)] < 1):
        raise ValueError(" [!] Found one or more spell with starting time < 1")
    reversed_rows = np.flatnonzero(np.isfinite(begins) & np.isfinite(ends) & (ends < begins)) + 1
    if reversed_rows.size:
        indexes = " ".join(str(int(value)) for value in reversed_rows)
        raise ValueError(
            " [!] Found one or more spell with ending time < starting time\n" + indexes
        )

    pdata_mapping: dict[object, list[float | None]] | None = None
    if isinstance(pdata, pd.DataFrame):
        assert pvar is not None
        _, _, pdata_mapping = _resolve_pdata(pdata, pvar)

    if process:
        if _positive_numeric(ends).size == 0:
            warnings.warn(
                "No positive value in 'end' column!",
                UserWarning,
                stacklevel=2,
            )
        if isinstance(pdata, pd.DataFrame) or pdata == "auto":
            mode = "year2age"
        else:
            mode = "age2age"
        axis_start = 1
        width = int(limit)
        columns = [f"a{i}" for i in range(1, width + 1)]
    else:
        if not isinstance(pdata, pd.DataFrame):
            pdata = None
            pdata_mapping = None
        mode = "age2year" if isinstance(pdata, pd.DataFrame) else "year2year"
        if tmin is None or tmax is None:
            if mode == "year2year":
                positive_begin = _positive_numeric(begins)
                positive_end = _positive_numeric(ends)
                if positive_begin.size == 0 or positive_end.size == 0:
                    raise ValueError(
                        " [!] Could not find the minimum or maximum year automatically, "
                        "please use tmin/tmax options"
                    )
                tmin = int(positive_begin.min())
                tmax = int(positive_end.max())
            else:
                assert pdata_mapping is not None
                births = np.asarray(
                    [
                        np.nan if value is None else value
                        for values in pdata_mapping.values()
                        for value in values
                    ],
                    dtype=float,
                )
                positive_births = _positive_numeric(births)
                positive_end = _positive_numeric(ends)
                if positive_births.size == 0 or positive_end.size == 0:
                    raise ValueError(
                        " [!] Could not find the minimum or maximum year automatically, "
                        "please use tmin/tmax options"
                    )
                tmin = int(positive_births.min())
                tmax = int(positive_end.max() + positive_births.max())
        assert tmin is not None and tmax is not None
        axis_start = int(tmin)
        width = int(tmax - tmin + 1)
        columns = [f"y{year}" for year in range(int(tmin), int(tmax) + 1)]

    if any(_cell_is_na(value) for value in seqdata[id_column].to_numpy(dtype=object)):
        raise ValueError("missing values in 'row.names' are not allowed")
    ids = list(pd.unique(seqdata[id_column]))
    raw = np.full((len(ids), width), None, dtype=object)
    truncated: list[int] = []
    empty_missing_start: list[int] = []
    start_before_birth: list[int] = []
    spell_before_birth: list[int] = []
    status_values = status_series.to_numpy(dtype=object)

    for sequence_index, sequence_id in enumerate(ids):
        source_positions = np.flatnonzero(seqdata[id_column].to_numpy() == sequence_id)
        local_begins = begins[source_positions].copy()
        local_ends = ends[source_positions].copy()
        local_states: list[object] = []
        for position in source_positions:
            state = status_values[position]
            text = None if _cell_is_na(state) else _r_state_text(state)
            local_states.append(None if text is None or text in missing else text)

        first_begin = local_begins[0]
        birth: float | None = None
        if mode == "year2age":
            if isinstance(pdata, str) and pdata == "auto":
                if np.isfinite(first_begin):
                    birth = float(first_begin)
                    age1 = 0.0
                else:
                    age1 = np.nan
            else:
                assert pdata_mapping is not None
                if sequence_id not in pdata_mapping:
                    raise ValueError(
                        ' [>] pdata must either contain a birth year per individual or be set as "auto"'
                    )
                birth_values = pdata_mapping[sequence_id]
                if len(birth_values) != 1:
                    raise ValueError("the condition has length > 1")
                birth = birth_values[0]
                if birth is None:
                    age1 = np.nan
                else:
                    age1 = first_begin - birth
                    if np.isfinite(age1) and age1 < 0:
                        start_before_birth.append(sequence_index + 1)
                        if not np.isfinite(local_ends[0]):
                            raise ValueError("missing value where TRUE/FALSE needed")
                        if local_ends[0] >= birth:
                            local_begins[0] = float(birth)
                            age1 = 0.0
                        else:
                            spell_before_birth.append(sequence_index + 1)
        elif mode == "age2year":
            assert pdata_mapping is not None
            if sequence_id not in pdata_mapping:
                raise ValueError("argument is of length zero")
            birth_values = pdata_mapping[sequence_id]
            if len(birth_values) != 1:
                raise ValueError("the condition has length > 1")
            birth = birth_values[0]
            age1 = first_begin
        else:
            age1 = first_begin

        if not np.isfinite(age1):
            empty_missing_start.append(sequence_index + 1)
            continue
        if age1 < 0:
            continue

        spell_order = range(len(source_positions))
        if not overwrite:
            spell_order = reversed(range(len(source_positions)))
        was_truncated = False
        for local_index in spell_order:
            begin_value = local_begins[local_index]
            end_value = local_ends[local_index]
            if not np.isfinite(begin_value) or not np.isfinite(end_value):
                warnings.warn(
                    f"skipping episode {local_index + 1} for case {sequence_index + 1} "
                    "due to missing start/end time",
                    UserWarning,
                    stacklevel=2,
                )
                continue
            if mode == "age2year" and birth is None:
                warnings.warn(
                    f"skipping episode {local_index + 1} for case {sequence_index + 1} "
                    "due to missing start/end time",
                    UserWarning,
                    stacklevel=2,
                )
                continue
            sstart = _spell_position(
                begin_value,
                mode=mode,
                birth=birth,
                axis_start=axis_start,
            )
            sstop = _spell_position(
                end_value,
                mode=mode,
                birth=birth,
                axis_start=axis_start,
            )

            if fillblanks is not None and local_index > 0:
                previous_end = local_ends[local_index - 1]
                if not np.isfinite(previous_end):
                    raise ValueError("missing value where TRUE/FALSE needed")
                if mode in ("age2age", "age2year"):
                    previous_position = int(previous_end)
                elif mode == "year2age":
                    assert birth is not None
                    previous_position = int(previous_end - birth)
                else:
                    previous_position = int(previous_end - axis_start + 1)
                if sstart not in (previous_position, previous_position + 1):
                    duration = sstart - (previous_position + 1)
                    if duration > 0 and sstart - 1 < width and sstart > 0 and previous_end > 0:
                        lo = max(previous_position, 0)
                        hi = min(sstart - 1, width)
                        raw[sequence_index, lo:hi] = fillblanks

            state = local_states[local_index]
            duration = sstop - sstart + 1
            if duration < 0 or sstop <= 0 or sstart < 0 or state is None:
                continue
            if sstop <= width:
                lo = max(sstart - 1, 0)
                raw[sequence_index, lo:sstop] = state
            else:
                if not was_truncated:
                    truncated.append(sequence_index + 1)
                    was_truncated = True
                if sstart <= width:
                    lo = max(sstart - 1, 0)
                    raw[sequence_index, lo:width] = state

    if start_before_birth:
        warnings.warn(
            f"{len(start_before_birth)} cases with a spell starting before birth year. "
            "Indexes in attribute 'issues'",
            UserWarning,
            stacklevel=2,
        )
    if spell_before_birth:
        warnings.warn(
            f"{len(spell_before_birth)} empty sequences because of spell occurring before "
            "birth year. Indexes in attribute 'issues'",
            UserWarning,
            stacklevel=2,
        )
    if empty_missing_start:
        warnings.warn(
            f"{len(empty_missing_start)} empty sequences because of missing start time. "
            "Indexes in attribute 'issues'",
            UserWarning,
            stacklevel=2,
        )
    if truncated:
        warnings.warn(
            f"{len(truncated)} sequences truncated at limit {width} . "
            "Indexes in attribute 'issues'",
            UserWarning,
            stacklevel=2,
        )
    return _materialize_sts(
        raw,
        ids=ids,
        columns=columns,
        real_alphabet=real_alphabet,
        start=axis_start,
        issues={
            "truncated": truncated,
            "empty.seq": empty_missing_start,
            "start.before.birth": start_before_birth,
            "spell.before.birth": spell_before_birth,
        },
    )


def seqdef_from_long(
    data: pd.DataFrame,
    *,
    id_col: str = "id",
    time_col: str = "time",
    state_col: str = "state",
    time_grid: Sequence[object] | None = None,
    **kwargs,
) -> STSList:
    """Create STSList from long format DataFrame (id, time, state).

    One row per person, in order of first appearance, as R's route for long data
    (``reshape`` to wide, then ``seqdef``) gives them (TX-AUD-20261019-039).
    A missing id is refused, as R's ``seqdef`` refuses it as a row name.
    There is one position per observed time value, and the time values become
    ``xtlab`` (an explicit ``xtlab=`` wins). A period no person reports is not
    invented: wave-coded times such as 1, 3, 5 have no period 2. Give
    ``time_grid`` (e.g. ``range(2001, 2006)``) to place the sequences on a
    fixed grid; a period of the grid nobody reports then follows the rules
    below like any unobserved period, and a time outside the grid, or a grid
    entry given twice, is refused (TX-AUD-20261010-018). The labels are R's
    ``as.character`` of the times: a date, or a date-time at midnight, reads
    as the date alone.

    A period with no row for a person, or with a missing state, is decided by
    ``seqdef``'s default rules: leading and internal ones become the ``nr``
    marker, trailing ones the ``void`` marker. The states keep their own text,
    a number read as R prints it (an integer state stays ``"1"`` even when a
    period is unobserved, and ``1.0`` is ``"1"``, as in ``seqdef``;
    TX-AUD-20261010-073). A state
    spelled like a marker (``"*"`` by default) stays a real state, and the
    marker takes another spelling when a cell needs it (TX-AUD-20261019-005).
    R's own route for long data, ``reshape`` to wide and then ``seqdef``, lists
    such a state yet counts its cells as missing: for persons A, *, B and A, A,
    B, R's ``seqlength(with.missing=FALSE)`` is 2 and 3, here 3 and 3.

    Positions follow the order of the time column: numbers and dates in value
    order, an ordered categorical in its category order, and plain strings in
    text order -- so labels such as ``"t1"`` ... ``"t12"`` need an ordered
    categorical (or zero padding) to keep ``"t10"`` after ``"t9"``.

    As in :func:`seqdef`, the alphabet takes R's ``seqstatl`` order (numeric
    when every state is an R number, R's collation otherwise), and ``labels``
    and ``cpal`` are read in that order, so the same data gets the same
    alphabet, labels and colours as in ``seqdef``: one entry per real state
    (R's length) or per codebook entry, by default the states and ``seqdef``'s
    palette, the markers taking their own label and colour (TX-AUD-20261019-035,
    TX-REP-20261019-013). There is no ``alphabet=``:
    ``seqdef(x, alphabet=[...])`` on the result sets another order and keeps
    every cell and marker.
    """
    validate_seq_kwargs(kwargs, func_name="seqdef_from_long")
    if isinstance(data, pd.DataFrame) and state_col in data.columns:
        # The pivot inserts NaN for unobserved periods; an object column keeps
        # integer states from being upcast to float text ("1.0", audit S2-05).
        data = data.astype({state_col: object})
    wide = cast(
        pd.DataFrame,
        seqformat(
            data,
            from_="long",
            to="wide",
            id_col=id_col,
            time_col=time_col,
            state_col=state_col,
        ),
    )
    # R's route (reshape to wide, then seqdef with the ids as row names) keeps a
    # person with no id until seqdef, which stops.
    if bool(wide[id_col].isna().any()):
        raise ValueError(
            "missing values in 'row.names' are not allowed: the id column has missing values"
        )
    # The pivot sorts the persons; R's route keeps them in order of first
    # appearance (TX-AUD-20261019-039). One positional pass: linear in persons.
    order = pd.Index(wide[id_col]).get_indexer(pd.unique(data[id_col]))
    wide = wide.iloc[order[order >= 0]].reset_index(drop=True)
    state_wide = wide.drop(columns=[id_col], errors="raise")
    if time_grid is not None:
        grid = list(time_grid)
        index = pd.Index(grid)
        if index.has_duplicates:
            repeated = index[index.duplicated()].unique().tolist()
            raise ValueError(f"time_grid has duplicate entries: {repeated}")
        outside = [time for time in state_wide.columns if time not in grid]
        if outside:
            raise ValueError(f"time values outside time_grid: {outside}")
        state_wide = state_wide.reindex(columns=grid)
    output_kwargs = dict(kwargs)
    output_kwargs.setdefault("ids", wide[id_col].tolist())
    output_kwargs.setdefault("xtlab", [str(r_as_character(time)) for time in state_wide.columns])
    return build_sts_from_state_matrix(state_wide.to_numpy(), kwargs=output_kwargs)


@_require_stslist
def seqstart(
    x: STSList,
    data_start: int | Sequence[int],
    new_start: int | Sequence[int],
    *,
    tmin: int | None = None,
    tmax: int | None = None,
    missing: str | None = None,
) -> STSList:
    """Align sequences with the ``TraMineRextras::seqstart`` placement algorithm.

    Scalar start values are recycled by row. The public R algorithm computes
    one-based ``new.index = data.start - new.start + 1``, applies its negative
    index correction, adjusts or derives ``tmin``/``tmax``, and clips copied
    source cells to the resulting width.

    R returns an unlabeled character matrix. Traminex deliberately re-expresses
    those cells as an :class:`STSList`: state, identity, weight, display, and
    issue metadata are copied by value. ``start`` and ``xtlab`` use the
    zero-based placement axis mechanically derived as
    ``corrected_tmin - correction - 1``. R's default ``NA`` fill is represented
    by ``x.nr`` unless an explicit ``missing`` state is supplied; when a real
    state is spelled like ``x.nr`` (a literal ``"*"``), the fill marker takes a
    spelling no state uses (``"<NR>"``), as ``seqdef`` does, so those cells stay
    real. R's ``seqdata`` may be any data frame or matrix; traminex takes an
    :class:`STSList`.

    As R, it stops when the frame ends before position 0 ("There are no data
    in the specified new time frame"), when ``x`` has no sequences and
    ``tmax`` is not given (R's maximum of no rows is ``-Inf``) or ``tmin`` is
    not given, and when ``tmin`` is more than one position after ``tmax``
    (R's ``matrix()`` refuses a negative number of columns); ``tmin`` one
    after ``tmax`` gives sequences of length 0, as in R.
    """
    # Both marker representations are equivalent inputs (identity for
    # in-alphabet objects; cycle-9 audit N-06) — the placement below copies
    # source codes verbatim while reassigning nr_code, which the result's
    # constructor rejects for external codes.
    x = x._canonicalized()
    n, L = x.nseq, x.seqlen

    def start_vector(value: int | Sequence[int]) -> np.ndarray:
        values = np.asarray(value)
        if values.ndim == 0:
            return np.full(n, int(values), dtype=int)
        if values.ndim != 1:
            raise ValueError("start values must be scalars or one-dimensional row vectors")
        if values.size == 1:
            return np.full(n, int(values[0]), dtype=int)
        if values.size != n:
            raise ValueError("number of rows and data_start and/or new_start values mismatch")
        return values.astype(int)

    data_start_values = start_vector(data_start)
    new_start_values = start_vector(new_start)
    new_index = data_start_values - new_start_values + 1
    if np.any(new_index < 1):
        correction = -int(np.min(new_index)) + 1
        new_index = new_index + correction
    else:
        correction = 0

    # In R's order (seqstart.R:18-31). With no rows R's max of the empty index is
    # -Inf, so the default tmax stops with "no data"; its min is Inf, so a
    # default tmin leaves a frame matrix() refuses (TX-AUD-20261011-062).
    no_data = "There are no data in the specified new time frame"
    if n == 0 and tmax is None:
        raise ValueError(no_data)
    corrected_tmax = (int(np.max(new_index)) + L - 1 if tmax is None else int(tmax)) + correction
    if corrected_tmax < 0:
        raise ValueError(no_data)
    if n == 0 and tmin is None:
        raise ValueError(
            "invalid 'ncol' value (too large or NA): with no sequences tmin has no default"
        )
    corrected_tmin = (max(int(np.min(new_index)), 1) if tmin is None else int(tmin)) + correction

    width = corrected_tmax - corrected_tmin + 1
    if width < 0:
        raise ValueError(
            "invalid 'ncol' value (< 0): the new time frame from tmin "
            f"{corrected_tmin - correction} to tmax {corrected_tmax - correction} "
            "has a negative length"
        )
    alphabet = list(x.alphabet)
    nr = x.nr
    fill_code: int | None
    if missing is None and x.nr_code is not None:
        fill, fill_code = x.nr, int(x.nr_code)
    elif missing is None:
        # R fills with NA, which no state can equal: a real state spelled like
        # the nr marker keeps its cells, and the fill gets a spelling no state
        # uses, as seqdef gives it (audit S2-09).
        occupied = [state for code, state in enumerate(alphabet) if code != x.void_code]
        fill = _collision_free_marker("*" if x.nr is None else str(x.nr), occupied, "NR")
        fill_code = None
    else:
        fill = missing
        fill_code = next((i for i, state in enumerate(alphabet) if state == fill), None)
    labels = None if x.labels is None else list(x.labels)
    cpal = None if x.cpal is None else list(x.cpal)
    data = np.full((n, width), -1 if fill_code is None else fill_code, dtype=int)
    target_starts = new_index - corrected_tmin
    for row in range(n):
        for source_column in range(L):
            target_column = int(target_starts[row]) + source_column
            if 0 <= target_column < width:
                data[row, target_column] = x.data[row, source_column]

    if fill_code is None and np.any(data == -1):
        fill_code = len(alphabet)
        alphabet.append(fill)
        data[data == -1] = fill_code
        if missing is None:
            nr = fill
        if labels is not None:
            labels.append(fill)
        if cpal is not None:
            fill_color = (
                x.missing_color if fill == nr else ("white" if fill == x.void else "#808080")
            )
            cpal.append(fill_color)

    nr_code = fill_code if fill == nr and fill_code is not None else x.nr_code
    void_code = fill_code if fill == x.void and fill_code is not None else x.void_code
    axis_start = corrected_tmin - correction - 1

    return STSList(
        data=data,
        alphabet=alphabet,
        labels=labels,
        cpal=cpal,
        xtlab=[str(axis_start + position) for position in range(width)],
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        nr=nr,
        void=x.void,
        nr_code=nr_code,
        void_code=void_code,
        weights=None if x.weights is None else x.weights.copy(),
        start=axis_start,
        ids=None if x.ids is None else list(x.ids),
        missing_color=x.missing_color,
        issues=(
            None
            if x.issues is None
            else {name: list(indexes) for name, indexes in x.issues.items()}
        ),
    )


@_require_stslist
def toPersonPeriod(x: STSList, *, time_col: str = "timestamp") -> pd.DataFrame:
    """Convert an STSList to TraMineRextras person-period form.

    The returned columns are ``id``, ``state``, and ``time_col``. Sequence
    identities come from ``x.ids`` when present, time is always zero-based,
    structural void cells are omitted, and nonresponse cells become true
    missing values. Rows follow R's ID/time/state ordering, the ids as R's
    collation orders character row names (``_z, a, b, C``, not code points).
    R's factor state is re-expressed as object values plus ``pd.NA``; an
    STSList without IDs uses one-based integer IDs while retaining R's
    character-row-name sort order.
    """
    if time_col in {"id", "state"}:
        raise ValueError("time_col must differ from 'id' and 'state'")

    custom_ids = x.ids is not None
    if custom_ids:
        id_values = np.empty(x.nseq, dtype=object)
        id_values[:] = x.ids
    else:
        id_values = np.arange(1, x.nseq + 1, dtype=int)

    codes = x.data.ravel(order="C")
    states = np.asarray([x._state_for_code(int(code)) for code in codes], dtype=object)
    ids: np.ndarray = np.repeat(id_values, x.seqlen)
    times = np.tile(np.arange(x.seqlen, dtype=int), x.nseq)

    keep = np.ones(codes.shape, dtype=bool)
    if x.void_code is not None:
        keep &= codes != x.void_code
    missing = np.zeros(codes.shape, dtype=bool)
    if x.nr_code is not None:
        missing = codes == x.nr_code

    states = states[keep]
    states[missing[keep]] = pd.NA
    state_order = codes[keep].astype(float)
    state_order[missing[keep]] = np.nan
    ids = ids[keep]
    times = times[keep]
    result = pd.DataFrame(
        {
            "id": ids,
            "state": states,
            time_col: times,
        },
        columns=["id", "state", time_col],
    )

    # toPersonPeriod.R:15 orders the character row names by R's collation,
    # not by code points (audit S2-07).
    id_text = result["id"].map(str)
    order = pd.DataFrame(
        {
            "id": id_text.map(r_character_rank(id_text.unique().tolist())),
            "time": times,
            "state": state_order,
        }
    )
    row_order = order.sort_values(
        ["id", "time", "state"],
        kind="stable",
        na_position="last",
    ).index
    return result.iloc[row_order].reset_index(drop=True)


@_require_stslist
def STS_to_TSE(
    x: STSList, *, state_changes_only: bool = True, time_col: str = "time"
) -> pd.DataFrame:
    """Convert STSList to a time-stamped event/state table (id, time, state).

    - If state_changes_only=True (default), only emit records when the state changes
      (including the first position). If False, emit one record per time step.
    - time_col controls the name of the time column ('time' by default).

    The ``id`` column carries the object's ids, and 1..n when it has none
    (TX-AUD-20261010-033). R's STS_to_TSE numbers the sequences unless an
    ``id`` is passed, and ``seqformat(to="TSE")`` then warns that it replaces
    the original ids; this function has no ``id`` argument and keeps them.

    The time column is ``x.start + position``: 1-based by default, as
    :func:`TSE_to_STS` reads it. A void cell (the end of a shorter sequence) is
    no observed state and gets no record, so a sequence of void cells only has
    none, and a table without records still has its three columns. With
    ``state_changes_only=False`` (one record per cell), :func:`TSE_to_STS` reads
    the table back with its void cells as void, except that it has no column
    for a position no sequence reaches and no row for a sequence without
    records; a change table does not read back. A missing (nr) cell keeps its
    record with the nr marker as its state, which ``TSE_to_STS`` reads as a
    state of that name. With ``state_changes_only=True`` the state after an
    internal void starts a new run (TX-AUD-20261019-082).

    This remains the Python state-at-time convenience API. For TraMineR's
    event-driven conversion, use ``seqformat(..., to="TSE", tevent=...)``,
    which records transitions (into void and missing too) at 0-based times.
    """
    import pandas as pd

    ids: list[object] = list(x.ids) if x.ids is not None else list(range(1, x.nseq + 1))
    rows = []
    for i in range(x.nseq):
        row = x.data[i]
        prev = None
        for t in range(x.seqlen):
            cur = int(row[t])
            if x.void_code is not None and cur == x.void_code:
                prev = None
                continue
            is_change = prev is None or cur != prev
            if (state_changes_only and is_change) or (not state_changes_only):
                rows.append(
                    {
                        "id": ids[i],
                        time_col: int(x.start + t),
                        "state": x._state_for_code(cur),
                    }
                )
            prev = cur
    # Rows are built sequence by sequence and time by time, R's order; sorting
    # by id would reorder sequences by their id text. The columns are named so
    # a table without records keeps them (batch-11 review M-8).
    return pd.DataFrame(rows, columns=["id", time_col, "state"])


def to_SPS(
    x: STSList,
    *,
    state_sep: str = "-",
    sdsep: str = "/",
    compress: bool = True,
) -> list[str]:
    """Write sequences to a simple SPS-like string format.

    - If compress=True, runs are collapsed and printed as state + sdsep + duration.
      Example: A/2-B/1-C/1 for A,A,B,C.
    - If compress=False, prints all states with state_sep.

    Void cells are left out wherever they are, so this is R's SPS form of the
    sequence without its void cells; a missing cell is a state
    (TX-AUD-20261010-033). R's ``seqformat`` (and :func:`seqformat`) writes an
    internal void, which only a recode can create, as a missing spell instead.
    """
    out: list[str] = []
    for row in x.data:
        states = [x._state_for_code(int(code)) for code in row if int(code) != x.void_code]
        if not compress:
            out.append(state_sep.join([str(s) for s in states]))
            continue
        # run-length encode
        tokens: list[str] = []
        if len(states) == 0:
            out.append("")
            continue
        cur = states[0]
        dur = 1
        for s in states[1:]:
            if s == cur:
                dur += 1
            else:
                tokens.append(f"{cur}{sdsep}{dur}")
                cur = s
                dur = 1
        tokens.append(f"{cur}{sdsep}{dur}")
        out.append(state_sep.join(tokens))
    return out


def to_DSS(x: STSList, *, state_sep: str = "-") -> list[str]:
    """Write sequences to a simple DSS-like format (distinct successive states).

    Void cells are left out wherever they are, so this is R's DSS form of the
    sequence without its void cells; a missing cell is a state
    (TX-AUD-20261010-033). R's ``seqformat`` (and :func:`seqformat`) writes an
    internal void, which only a recode can create, as a missing state instead.
    """
    out: list[str] = []
    for row in x.data:
        states = [x._state_for_code(int(code)) for code in row if int(code) != x.void_code]
        tokens: list[str] = []
        prev = None
        for s in states:
            if prev is None or s != prev:
                tokens.append(str(s))
                prev = s
        out.append(state_sep.join(tokens))
    return out
