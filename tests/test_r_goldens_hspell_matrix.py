"""Whole-matrix R goldens for horizontal ``HSPELL_to_STS``.

Generator: ``scripts/r_goldens/make_goldens_hspell_horizontal.R`` (the
``hspell_horizontal_matrix_cases.txt`` section), live TraMineRextras 0.6.9.

The earlier goldens were single-row cases with complete spells, so they could
not see that traminex sorted rows by id text (``1, 10, 11, 12, 2, ...``),
dropped a spell whose end is missing instead of running it to ``tmax``, filled
uncovered cells with ``"_"`` instead of leaving them NA, applied ``endObs`` one
period late and wrote the filler instead of NA, applied one ``fixed_status``
to every spell, and rendered integer statuses as ``"1.0"`` (2026-09-22 audit
S2-01..S2-04). R returns NA cells in a data frame; the traminex result is that
frame after ``seqdef``, so an NA cell is the missing or void marker here.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "hspell_horizontal_matrix_cases.txt"


def _load() -> dict[str, dict]:
    cases: dict[str, dict] = defaultdict(lambda: {"row_names": [], "rows": []})
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, _, row_name, columns, cells = line.split("\t")
        cases[case]["columns"] = columns.split(",")
        cases[case]["row_names"].append(row_name)
        cases[case]["rows"].append(cells.split(","))
    return dict(cases)


R = _load()

TWO_SPELLS = pd.DataFrame(
    {
        "id": ["p1", "p2"],
        "begin1": [1, 1],
        "end1": [3, 2],
        "status1": ["A", "A"],
        "begin2": [2, 3],
        "end2": [4, 4],
        "status2": ["B", "B"],
    }
)
SPELLS = {"begin": ["begin1", "begin2"], "end": ["end1", "end2"], "status": ["status1", "status2"]}
NA_ENDS = pd.DataFrame(
    {
        "id": [1, 2, 3],
        "b1": [1, 2, 1],
        "e1": [2, np.nan, 3],
        "s1": ["A", "B", "A"],
        "b2": [3, np.nan, 4],
        "e2": [np.nan, np.nan, 5],
        "s2": ["C", "D", "C"],
    }
)
NA_SPELLS = {"begin": ["b1", "b2"], "end": ["e1", "e2"]}

CASES = {
    "fillblanks": lambda: tx.HSPELL_to_STS(TWO_SPELLS, id_col="id", fillblanks="_", **SPELLS),
    "no_id": lambda: tx.HSPELL_to_STS(TWO_SPELLS, **SPELLS),
    "overwrite_false": lambda: tx.HSPELL_to_STS(TWO_SPELLS, id_col="id", overwrite=False, **SPELLS),
    "explicit_bounds": lambda: tx.HSPELL_to_STS(TWO_SPELLS, id_col="id", tmin=0, tmax=6, **SPELLS),
    "endobs_scalar": lambda: tx.HSPELL_to_STS(TWO_SPELLS, id_col="id", endObs=3, **SPELLS),
    "endobs_column": lambda: tx.HSPELL_to_STS(
        TWO_SPELLS.assign(obs=[3, 4]), id_col="id", fillblanks="_", endObs="obs", **SPELLS
    ),
    "row_order": lambda: tx.HSPELL_to_STS(
        pd.DataFrame({"id": range(1, 13), "b": 1, "e": 3, "s": list("ABCDEFGHIJKL")}),
        id_col="id",
        begin=["b"],
        end=["e"],
        status=["s"],
    ),
    "missing_ends": lambda: tx.HSPELL_to_STS(
        NA_ENDS, id_col="id", status=["s1", "s2"], **NA_SPELLS
    ),
    "fixed_status_vector": lambda: tx.HSPELL_to_STS(
        NA_ENDS, id_col="id", fixed_status=["X", "Y"], **NA_SPELLS
    ),
    "fixed_status_scalar": lambda: tx.HSPELL_to_STS(
        NA_ENDS, id_col="id", fixed_status="X", **NA_SPELLS
    ),
    "numeric_status": lambda: tx.HSPELL_to_STS(
        pd.DataFrame(
            {
                "id": [1, 2],
                "b1": [0, 1],
                "e1": [2, 2],
                "s1": [1, 2],
                "b2": [3, 3],
                "e2": [5, 5],
                "s2": [2.5, np.nan],
            }
        ),
        id_col="id",
        begin=["b1", "b2"],
        end=["e1", "e2"],
        status=["s1", "s2"],
    ),
    "help_example": lambda: tx.HSPELL_to_STS(
        pd.DataFrame(
            {
                "begin1": [1] * 5,
                "end1": [2, 3, 4, 5, np.nan],
                "status1": [1, 2, 3, 4, 5],
                "begin2": [3, 4, 5, 6, np.nan],
                "end2": [np.nan] * 5,
                "status2": [5, 4, 3, 2, 1],
            }
        ),
        begin=["begin1", "begin2"],
        end=["end1", "end2"],
        status=["status1", "status2"],
    ),
}


def _as_r_frame(x: tx.STSList) -> list[list[str]]:
    """Decode to R's view: the missing and void markers are both NA."""
    markers = {x.nr_code, x.void_code} - {None}
    return [["NA" if code in markers else x.alphabet[code] for code in row] for row in x.data]


# R blanks one period early under endObs (HSPELL_to_STS.R:25,67-68); traminex
# blanks from endObs on (TX-REP-20261017-005), so these two cases keep R's ids and
# columns but assert the corrected rows and that they differ from R's.
ENDOBS_CORRECTED = {
    "endobs_scalar": [["A", "B", "NA", "NA"], ["A", "A", "NA", "NA"]],
    "endobs_column": [["A", "B", "NA", "NA"], ["A", "A", "B", "NA"]],
}


@pytest.mark.parametrize("case", sorted(set(CASES) - set(ENDOBS_CORRECTED)))
def test_horizontal_spells_reproduce_the_r_matrix(case):
    result = CASES[case]()

    assert result.ids == R[case]["row_names"]
    assert list(result.xtlab) == R[case]["columns"]
    assert _as_r_frame(result) == R[case]["rows"]


@pytest.mark.parametrize("case", sorted(ENDOBS_CORRECTED))
def test_endobs_cases_blank_from_endobs_on_where_r_blanks_one_period_early(case):
    result = CASES[case]()

    assert result.ids == R[case]["row_names"]
    assert list(result.xtlab) == R[case]["columns"]
    assert _as_r_frame(result) == ENDOBS_CORRECTED[case]
    assert ENDOBS_CORRECTED[case] != R[case]["rows"]


def test_the_open_spell_with_tmin_below_one_is_a_deliberate_divergence():
    """R stops an open spell ``1 - tmin`` periods early; traminex runs it to tmax.

    R sets the missing end to the unshifted ``tmax`` after shifting every time by
    ``1 - tmin`` (``HSPELL_to_STS.R:37``). With ``tmin = 0`` the open spell ``B``
    from time 1 therefore ends at time 2 in R, one period before ``tmax = 3``.
    """
    frame = pd.DataFrame({"id": [1, 2], "b1": [0, 1], "e1": [3, np.nan], "s1": ["A", "B"]})
    assert R["open_end_tmin_zero"]["rows"] == [["A", "A", "A", "A"], ["NA", "B", "B", "NA"]]

    result = tx.HSPELL_to_STS(frame, id_col="id", begin=["b1"], end=["e1"], status=["s1"])

    assert _as_r_frame(result) == [["A", "A", "A", "A"], ["NA", "B", "B", "B"]]


def test_duplicate_ids_are_refused_rather_than_merged():
    frame = pd.DataFrame({"id": ["a", "a"], "b": [1, 2], "e": [2, 3], "s": ["X", "Y"]})
    with pytest.raises(ValueError, match="id values must be unique"):
        tx.HSPELL_to_STS(frame, id_col="id", begin=["b"], end=["e"], status=["s"])


def test_a_float_column_does_not_turn_vertical_ids_into_float_text():
    """``iterrows`` coerced every row to float when any column was float (S2-03)."""
    frame = pd.DataFrame(
        {
            "id": [7, 7, 8],
            "start": [1, 3, 1],
            "end": [2, 4, 4],
            "state": [1, 3, 1],
            "w": [0.5, 0.5, 1.0],
        }
    )
    result = tx.HSPELL_to_STS(
        frame, id_col="id", start_col="start", end_col="end", state_col="state"
    )

    assert result.ids == ["7", "8"]
    assert sorted(result.alphabet) == ["1", "3"]
