"""Public-R goldens for horizontal ``TraMineRextras::HSPELL_to_STS``.

Generator: ``scripts/r_goldens/make_goldens_hspell_horizontal.R`` against
TraMineRextras 0.6.9.
"""

from __future__ import annotations

import pandas as pd
import pytest
from r_goldens_utils import load_golden_lines

from traminex import HSPELL_to_STS

_CASES = {
    "right_tail": {
        "data": {
            "id": ["right"],
            "b1": [1],
            "e1": [3],
            "s1": ["A"],
            "b2": [2],
            "e2": [4],
            "s2": ["B"],
        },
        "begin": ["b1", "b2"],
        "end": ["e1", "e2"],
        "status": ["s1", "s2"],
        "overwrite": True,
    },
    "right_tail_no_overwrite": {
        "data": {
            "id": ["first"],
            "b1": [1],
            "e1": [3],
            "s1": ["A"],
            "b2": [2],
            "e2": [4],
            "s2": ["B"],
        },
        "begin": ["b1", "b2"],
        "end": ["e1", "e2"],
        "status": ["s1", "s2"],
        "overwrite": False,
    },
    "out_of_chronological_order": {
        "data": {
            "id": ["ordered"],
            "b1": [3],
            "e1": [4],
            "s1": ["B"],
            "b2": [1],
            "e2": [3],
            "s2": ["A"],
        },
        "begin": ["b1", "b2"],
        "end": ["e1", "e2"],
        "status": ["s1", "s2"],
        "overwrite": True,
    },
    "nested": {
        "data": {
            "id": ["nested"],
            "b1": [1],
            "e1": [6],
            "s1": ["A"],
            "b2": [2],
            "e2": [4],
            "s2": ["B"],
        },
        "begin": ["b1", "b2"],
        "end": ["e1", "e2"],
        "status": ["s1", "s2"],
        "overwrite": True,
    },
    "chained": {
        "data": {
            "id": ["chain"],
            "b1": [1],
            "e1": [5],
            "s1": ["A"],
            "b2": [2],
            "e2": [3],
            "s2": ["B"],
            "b3": [3],
            "e3": [6],
            "s3": ["C"],
        },
        "begin": ["b1", "b2", "b3"],
        "end": ["e1", "e2", "e3"],
        "status": ["s1", "s2", "s3"],
        "overwrite": True,
    },
    "same_bounds": {
        "data": {
            "id": ["tie"],
            "b1": [1],
            "e1": [4],
            "s1": ["A"],
            "b2": [1],
            "e2": [4],
            "s2": ["B"],
        },
        "begin": ["b1", "b2"],
        "end": ["e1", "e2"],
        "status": ["s1", "s2"],
        "overwrite": True,
    },
    "adjacent": {
        "data": {
            "id": ["adjacent"],
            "b1": [1],
            "e1": [2],
            "s1": ["A"],
            "b2": [3],
            "e2": [4],
            "s2": ["B"],
        },
        "begin": ["b1", "b2"],
        "end": ["e1", "e2"],
        "status": ["s1", "s2"],
        "overwrite": True,
    },
}


def _r_cases() -> dict[str, dict[str, str]]:
    cases = {}
    for line in load_golden_lines("hspell_horizontal_cases"):
        if line.startswith("#") or not line:
            continue
        case, row_id, spell_order, overwrite, columns, cells = line.split("\t")
        cases[case] = {
            "id": row_id,
            "spell_order": spell_order,
            "overwrite": overwrite,
            "columns": columns,
            "cells": cells,
        }
    return cases


@pytest.mark.parametrize("case", list(_CASES))
def test_hspell_horizontal_cells_match_public_r_by_id_and_position(case):
    spec = _CASES[case]
    oracle = _r_cases()[case]
    result = HSPELL_to_STS(
        pd.DataFrame(spec["data"]),
        id_col="id",
        begin=spec["begin"],
        end=spec["end"],
        status=spec["status"],
        overwrite=spec["overwrite"],
    )

    assert oracle["id"] in result.ids
    row_index = result.ids.index(oracle["id"])
    decoded = [str(result.alphabet[int(value)]) for value in result.data[row_index]]
    r_columns = oracle["columns"].split(",")
    expected = dict(zip(r_columns, oracle["cells"].split(","), strict=True))
    actual = dict(zip(r_columns, decoded, strict=True))

    assert actual == expected
    assert oracle["spell_order"] == ">".join(
        f"{begin}:{end}:{status}"
        for begin, end, status in zip(spec["begin"], spec["end"], spec["status"], strict=True)
    )
    assert oracle["overwrite"] == str(spec["overwrite"]).upper()


def test_hspell_pvar_alignment_and_fill_match_public_r():
    rows = []
    oracle: dict[str, tuple[list[str], list[str]]] = {}
    for line in load_golden_lines("hspell_pvar_case"):
        if line.startswith("#") or not line:
            continue
        row_id, pvar, begin, end, status, columns, cells = line.split("\t")
        rows.append(
            {
                "id": row_id,
                "origin": int(pvar),
                "begin1": int(begin),
                "end1": int(end),
                "status1": status,
            }
        )
        oracle[row_id] = (columns.split(","), cells.split(","))

    result = HSPELL_to_STS(
        pd.DataFrame(rows),
        id_col="id",
        begin=["begin1"],
        end=["end1"],
        status=["status1"],
        pvar="origin",
        fillblanks="_",
    )
    for row_id, (columns, expected_cells) in oracle.items():
        row_index = result.ids.index(row_id)
        decoded = [str(result.alphabet[int(value)]) for value in result.data[row_index]]
        assert dict(zip(columns, decoded, strict=True)) == dict(
            zip(columns, expected_cells, strict=True)
        )
