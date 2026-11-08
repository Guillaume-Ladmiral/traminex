"""Public-R goldens for ``TraMineRextras::seqstart``.

The R function returns a character matrix. Traminex deliberately re-expresses
those cells as an ``STSList`` while preserving Python identity/state metadata.

Generator: ``scripts/r_goldens/make_goldens_seqstart.R`` against
TraMineRextras 0.6.9.
"""

from __future__ import annotations

import numpy as np
import pytest
from r_goldens_utils import load_golden_lines

from traminex import STSList, seqstart


def _parse_vector(value: str) -> int | list[int]:
    values = [int(part) for part in value.split(",")]
    return values[0] if len(values) == 1 else values


def _parse_bound(value: str) -> int | None:
    return None if value == "NULL" else int(value)


def _success_cases() -> dict[str, dict[str, object]]:
    cases: dict[str, dict[str, object]] = {}
    for line in load_golden_lines("seqstart_success_cases"):
        if line.startswith("#") or not line:
            continue
        case, data_start, new_start, tmin, tmax, missing, shape, cells = line.split("\t")
        cases[case] = {
            "data_start": _parse_vector(data_start),
            "new_start": _parse_vector(new_start),
            "tmin": _parse_bound(tmin),
            "tmax": _parse_bound(tmax),
            "missing": missing,
            "shape": tuple(int(part) for part in shape.split("x")),
            "cells": [row.split(",") for row in cells.split(";")],
        }
    return cases


def _error_cases() -> dict[str, dict[str, object]]:
    cases: dict[str, dict[str, object]] = {}
    for line in load_golden_lines("seqstart_error_cases"):
        if line.startswith("#") or not line:
            continue
        case, data_start, new_start, tmin, tmax, message = line.split("\t")
        cases[case] = {
            "data_start": _parse_vector(data_start),
            "new_start": _parse_vector(new_start),
            "tmin": _parse_bound(tmin),
            "tmax": _parse_bound(tmax),
            "message": message,
        }
    return cases


def _source() -> STSList:
    return STSList(
        data=np.arange(6, dtype=int).reshape(2, 3),
        alphabet=["A", "B", "C", "D", "E", "F"],
        ids=["row-1", "row-2"],
    )


def _decoded(x: STSList) -> list[list[object]]:
    alphabet = np.empty(len(x.alphabet), dtype=object)
    alphabet[:] = x.alphabet
    return alphabet[x.data].tolist()


@pytest.mark.parametrize("case", list(_success_cases()))
def test_seqstart_cells_shapes_and_axis_match_public_r(case):
    oracle = _success_cases()[case]
    source = _source()
    kwargs = {
        "tmin": oracle["tmin"],
        "tmax": oracle["tmax"],
    }
    if oracle["missing"] != "<DEFAULT_NA>":
        kwargs["missing"] = oracle["missing"]

    result = seqstart(
        source,
        data_start=oracle["data_start"],
        new_start=oracle["new_start"],
        **kwargs,
    )
    expected = [[source.nr if cell == "<NA>" else cell for cell in row] for row in oracle["cells"]]
    expected_start = -2 if case == "clip_left" else 0

    assert result.data.shape == oracle["shape"]
    assert _decoded(result) == expected
    assert result.start == expected_start
    assert result.xtlab == [
        str(position) for position in range(expected_start, expected_start + result.seqlen)
    ]
    assert len(result.xtlab) == result.seqlen


@pytest.mark.parametrize(
    ("case", "message"),
    [("length_mismatch", "mismatch"), ("negative_tmax", "no data")],
)
def test_seqstart_rejects_public_r_error_branches(case, message):
    source = _source()
    oracle = _error_cases()[case]
    with pytest.raises(ValueError, match=message):
        seqstart(
            source,
            data_start=oracle["data_start"],
            new_start=oracle["new_start"],
            tmin=oracle["tmin"],
            tmax=oracle["tmax"],
        )


def test_seqstart_preserves_stslist_identity_state_metadata_and_isolation():
    source = STSList(
        data=np.array([[0, 1, 2], [2, 1, 0]], dtype=int),
        alphabet=["A", "B", "%"],
        labels=["Alpha", "Beta", "Void"],
        cpal=["red", "blue", "white"],
        xtlab=["1998", "1999", "2000"],
        xtstep=2,
        tick_last=True,
        nr="?",
        void="%",
        nr_code=None,
        void_code=2,
        weights=np.array([0.5, 1.5]),
        start=1998,
        ids=["first", ("second", 2)],
        missing_color="#aabbcc",
        issues={"truncated": [1], "empty.seq": [2]},
    )

    result = seqstart(source, data_start=[1998, 2000], new_start=2000)

    assert result.ids == source.ids
    assert result.ids is not source.ids
    np.testing.assert_array_equal(result.weights, source.weights)
    assert result.weights is not source.weights
    assert result.nr == "?" and result.void == "%"
    assert result.nr_code == 3 and result.void_code == 2
    assert result.alphabet == ["A", "B", "%", "?"]
    assert result.labels == ["Alpha", "Beta", "Void", "?"]
    assert result.cpal == ["red", "blue", "white", "#aabbcc"]
    assert result.missing_color == "#aabbcc"
    assert result.xtstep == 2 and result.tick_last is True
    assert result.start == 0
    assert result.xtlab == ["0", "1", "2", "3", "4"]
    assert result.issues == source.issues
    assert result.issues is not source.issues
    assert result.issues["truncated"] is not source.issues["truncated"]

    result.issues["truncated"][0] = 2
    assert source.issues["truncated"] == [1]


def test_seqstart_preserves_non_string_state_values_without_string_coercion():
    source = STSList(
        data=np.array([[0, 1]], dtype=int),
        alphabet=[10, 20],  # type: ignore[list-item]
    )

    result = seqstart(source, data_start=1, new_start=1)

    assert result.alphabet == [10, 20]
    assert _decoded(result) == [[10, 20]]
