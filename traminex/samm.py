"""Sequence Analysis Multistate Model person-period frames."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from numbers import Integral

import numpy as np
import pandas as pd

from ._categorical_order import stable_sorted_unique
from ._rformat import r_as_character
from .sequences import STSList
from .sequences import alphabet as _real_state_alphabet


@dataclass
class SAMMResult:
    """Explicit Python representation of R's attributed ``SAMM`` data frame.

    TraMineRextras stores the person-period frame, the extracted sequence data,
    and SAMM metadata in one data frame plus R attributes.  Keeping those parts
    explicit avoids hidden pandas attributes while preserving their values.
    """

    frame: pd.DataFrame
    sequences: STSList
    spell: tuple[str, ...]
    typology: np.ndarray
    state_columns: tuple[str, ...]
    sublength: int


def _state_categories(x: STSList) -> list[str]:
    categories = list(x.alphabet)
    for marker in (x.nr, x.void):
        if marker not in categories:
            categories.append(marker)
    return categories


def _state_text_lookup(x: STSList) -> dict[int, str]:
    lookup = {code: state for code, state in enumerate(x.alphabet)}
    if x.nr_code is not None:
        lookup[x.nr_code] = x.nr
    if x.void_code is not None:
        lookup[x.void_code] = x.void
    unknown = sorted(set(np.asarray(x.data, dtype=int).ravel()) - set(lookup))
    if unknown:
        raise ValueError(f"sequence data contains unknown state codes: {unknown}")
    return lookup


def seqsamm(
    x: STSList,
    sublength: int,
    covar: pd.DataFrame | None = None,
) -> SAMMResult:
    """Build the TraMineRextras SAMM person-period frame.

    ``time`` and ``begin`` intentionally use R's 1-based positions relative to
    the sequence matrix, independently of ``x.start``.  As in the R function,
    windows begin only through ``seqlen - sublength`` and any window containing
    the declared void code is omitted.
    """
    if not isinstance(x, STSList):
        raise TypeError("x must be an STSList")
    if isinstance(sublength, (bool, np.bool_)) or not isinstance(sublength, Integral):
        raise TypeError("sublength must be an integer")
    window_width = int(sublength)
    if window_width < 2:
        raise ValueError("sublength must be at least 2")
    if window_width >= x.seqlen:
        raise ValueError("sublength must be smaller than the sequence width")
    if covar is not None:
        if not isinstance(covar, pd.DataFrame):
            raise TypeError("covar must be a pandas DataFrame")
        if len(covar) != x.nseq:
            raise ValueError("covar must have one row per sequence")
        if not covar.columns.is_unique:
            raise ValueError("covar columns must be unique")

    state_columns = tuple(f"s.{position}" for position in range(1, window_width + 1))
    reserved = {"id", "time", "begin", "spell.time", "transition", *state_columns}
    if covar is not None:
        overlap = [str(column) for column in covar.columns if column in reserved]
        if overlap:
            raise ValueError(f"covar overlaps reserved columns: {', '.join(overlap)}")

    rows: list[tuple[int, int, float, float, bool, np.ndarray]] = []
    source_rows: list[int] = []
    for sequence_index, sequence in enumerate(x.data):
        spell_begin = 1
        for start in range(x.seqlen - window_width):
            if start > 0 and sequence[start - 1] != sequence[start]:
                spell_begin = start + 1
            window = np.asarray(sequence[start : start + window_width], dtype=int)
            if x.void_code is not None and np.any(window == x.void_code):
                continue
            rows.append(
                (
                    sequence_index + 1,
                    start + 1,
                    float(spell_begin),
                    float(start + 1 - spell_begin),
                    bool(window[0] != window[1]),
                    window.copy(),
                )
            )
            source_rows.append(sequence_index)

    lookup = _state_text_lookup(x)
    categories = _state_categories(x)
    frame = pd.DataFrame(
        {
            "id": pd.Series([row[0] for row in rows], dtype="int64"),
            "time": pd.Series([row[1] for row in rows], dtype="int64"),
            "begin": pd.Series([row[2] for row in rows], dtype="float64"),
            "spell.time": pd.Series([row[3] for row in rows], dtype="float64"),
            "transition": pd.Series([row[4] for row in rows], dtype="bool"),
        }
    )
    for offset, column in enumerate(state_columns):
        frame[column] = pd.Categorical(
            [lookup[int(row[5][offset])] for row in rows], categories=categories
        )
    if covar is not None:
        replicated = covar.iloc[source_rows].reset_index(drop=True)
        frame = pd.concat([frame, replicated], axis=1)

    sequence_data = (
        np.vstack([row[5] for row in rows]).astype(int, copy=False)
        if rows
        else np.empty((0, window_width), dtype=int)
    )
    sequences = STSList(
        data=sequence_data,
        alphabet=list(x.alphabet),
        labels=None if x.labels is None else list(x.labels),
        cpal=None if x.cpal is None else list(x.cpal),
        xtlab=None,
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        nr=x.nr,
        void=x.void,
        nr_code=x.nr_code,
        void_code=x.void_code,
        weights=None,
        start=x.start,
        ids=None,
        missing_color=x.missing_color,
        issues=None,
    )
    return SAMMResult(
        frame=frame,
        sequences=sequences,
        # R: attr(ret, "spell") <- alphabet(seqdata) (seqsamm.R:47) — the real
        # states only, never the missing/void markers (cycle-9 audit N-08).
        spell=tuple(_real_state_alphabet(x)),
        typology=np.full(len(frame), "None", dtype=object),
        state_columns=state_columns,
        sublength=window_width,
    )


def _validate_samm(samm: SAMMResult) -> None:
    if not isinstance(samm, SAMMResult):
        raise TypeError("samm must be a SAMMResult")
    if not samm.state_columns or samm.state_columns[0] not in samm.frame.columns:
        raise ValueError("samm is missing its leading state column")
    if len(samm.frame) != samm.sequences.nseq or len(samm.typology) != len(samm.frame):
        raise ValueError("samm frame, sequences, and typology must have equal lengths")


def _subset_sequences(x: STSList, positions: np.ndarray) -> STSList:
    return STSList(
        data=x.data[positions].copy(),
        alphabet=list(x.alphabet),
        labels=None if x.labels is None else list(x.labels),
        cpal=None if x.cpal is None else list(x.cpal),
        xtlab=None if x.xtlab is None else list(x.xtlab),
        xtstep=x.xtstep,
        tick_last=x.tick_last,
        nr=x.nr,
        void=x.void,
        nr_code=x.nr_code,
        void_code=x.void_code,
        weights=None if x.weights is None else x.weights[positions].copy(),
        start=x.start,
        ids=None if x.ids is None else [x.ids[int(position)] for position in positions],
        missing_color=x.missing_color,
        issues=None,
    )


def seqsammseq(samm: SAMMResult, spell: str) -> STSList:
    """Return transition windows that exit ``spell`` from a SAMM result."""
    _validate_samm(samm)
    if not isinstance(spell, str):
        raise TypeError("spell must be a string")
    leading_state = samm.frame[samm.state_columns[0]]
    selected = (leading_state == spell) & samm.frame["transition"].astype(bool)
    positions = np.flatnonzero(selected.to_numpy(dtype=bool))
    return _subset_sequences(samm.sequences, positions)


def _typology_parts(typology: Sequence[object]) -> tuple[list[str], list[str]]:
    categorical: pd.Categorical | None = None
    if isinstance(typology, pd.Categorical):
        categorical = typology
    elif isinstance(typology, pd.Series) and isinstance(typology.dtype, pd.CategoricalDtype):
        categorical = typology.array

    # R stores as.character(typology) and names the indicator columns from it,
    # so a typology given as the numbers 1, 2 is "1", "2" -- never "1.0"
    # (2026-09-22 audit S5-12).
    if categorical is not None:
        raw_values = list(categorical)
        types = [_typology_text(value) for value in categorical.categories]
    elif isinstance(typology, str):
        raw_values = [typology]
        types = [typology]
    else:
        try:
            raw_values = list(typology)
        except TypeError as error:
            raise TypeError("typology must be a sequence") from error
        types = list(dict.fromkeys(_typology_text(value) for value in raw_values))
    if any(pd.isna(value) for value in raw_values):
        raise ValueError("typology must not contain missing values")
    return [_typology_text(value) for value in raw_values], types


def _typology_text(value: object) -> str:
    """One typology value as R's ``as.character`` writes it."""
    return str(r_as_character(value))


def seqsammeha(
    samm: SAMMResult,
    spell: str,
    typology: Sequence[object],
    persper: bool = True,
) -> pd.DataFrame:
    """Prepare competing-risk event-history data for exits from ``spell``."""
    _validate_samm(samm)
    if not isinstance(spell, str):
        raise TypeError("spell must be a string")
    if not isinstance(persper, (bool, np.bool_)):
        raise TypeError("persper must be a boolean")

    leading_state = samm.frame[samm.state_columns[0]]
    spell_mask = (leading_state == spell).to_numpy(dtype=bool)
    if not spell_mask.any():
        raise ValueError(f"spell {spell!r} has no observations in samm")
    transition_mask = spell_mask & samm.frame["transition"].to_numpy(dtype=bool)
    values, types = _typology_parts(typology)
    transition_count = int(transition_mask.sum())
    if len(values) != transition_count:
        raise ValueError(
            "typology must have one value per transition out of the selected spell "
            f"({transition_count})"
        )

    assigned = np.asarray(samm.typology, dtype=object).copy()
    assigned[transition_mask] = values
    period = samm.frame.loc[spell_mask].copy().reset_index(drop=True)
    period_values = [str(value) for value in assigned[spell_mask]]
    levels = [
        "None",
        *(str(value) for value in stable_sorted_unique(period_values) if value != "None"),
    ]
    period["SAMMtypology"] = pd.Categorical(period_values, categories=levels)
    maxima = period.groupby(["id", "begin"], sort=False, dropna=False)["spell.time"].transform(
        "max"
    )
    period["lastobs"] = (maxima == period["spell.time"]).astype(bool)
    for type_name in types:
        period[f"SAMM{type_name}"] = (period["SAMMtypology"].astype(object) == type_name).astype(
            float
        )
    if not bool(persper):
        period = period.loc[period["lastobs"]].reset_index(drop=True)
    return period
