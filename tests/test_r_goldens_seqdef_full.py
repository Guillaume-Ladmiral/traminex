"""TraMineR 2.2-14 goldens for the complete ``seqdef`` contract.

TraMineR stores real states in its ``alphabet`` attribute and keeps the
nonresponse/void factor levels separately.  ``STSList`` needs integer codes for
every cell, so the Python representation appends only the marker levels that
are actually used.  The metadata golden records both views explicitly.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import seqdef
from traminex.sequences import STSList
from traminex.sequences import alphabet as real_alphabet

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"
CELLS = pd.read_csv(
    GOLDEN_DIR / "seqdef_full_cells.csv",
    dtype=str,
    keep_default_na=False,
)

EXPECTED_CASES = {
    "sts_default_policies",
    "sts_all_del",
    "sts_replacement_full_meta",
    "sts_compressed",
    "states_mapping",
    "sps_default",
    "spell_limit4",
}
EXPECTED_METADATA_FIELDS = {
    "r_alphabet",
    "python_alphabet",
    "labels",
    "cpal",
    "ids",
    "weights",
    "start",
    "missing",
    "void",
    "nr",
    "xtlab",
    "xtstep",
    "tick_last",
    "missing_color",
}
EXPECTED_ERROR_CASES = {
    "duplicate_alphabet",
    "unknown_state",
    "weights_length",
    "weights_na",
    "weights_negative",
    "palette_length",
    "label_length",
    "void_type",
    "void_equals_left",
    "void_equals_gaps",
    "void_equals_right",
    "id_length",
    "id_duplicate",
    "id_missing",
}


def _pipe_records(path: Path, *, columns: int) -> list[list[str]]:
    records: list[list[str]] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("|", columns - 1)
        assert len(fields) == columns, f"invalid record in {path.name}: {line!r}"
        records.append(fields)
    return records


METADATA = {
    (case, field): value
    for case, field, value in _pipe_records(
        GOLDEN_DIR / "seqdef_full_metadata.txt",
        columns=3,
    )
}
ERRORS = {
    name: message.replace(r"\n", "\n")
    for kind, name, message in _pipe_records(
        GOLDEN_DIR / "seqdef_full_errors.txt",
        columns=3,
    )
    if kind == "ERROR"
}


def _vector(case: str, field: str) -> list[str] | None:
    value = METADATA[(case, field)]
    if value == "<NULL>":
        return None
    return value.split(";")


def _scalar(case: str, field: str) -> str:
    values = _vector(case, field)
    assert values is not None and len(values) == 1
    return values[0]


def _case(name: str) -> pd.DataFrame:
    case = CELLS.loc[CELLS["case"] == name].copy()
    assert not case.empty, f"missing cell golden {name!r}"
    case["row"] = case["row"].astype(int)
    case["position"] = case["position"].astype(int)
    return case.sort_values(["row", "position"], kind="stable")


def _decoded(x: STSList) -> list[list[str]]:
    return [[x.alphabet[int(code)] for code in row] for row in x.data]


def _expected_rows(name: str) -> list[list[str]]:
    case = _case(name)
    return [
        group.sort_values("position", kind="stable")["value"].tolist()
        for _, group in case.groupby("row", sort=True)
    ]


def _expected_marker_code(name: str, marker_field: str) -> int | None:
    marker = _scalar(name, marker_field)
    rows = _expected_rows(name)
    if not any(marker in row for row in rows):
        return None
    expected_alphabet = _vector(name, "python_alphabet")
    assert expected_alphabet is not None
    return expected_alphabet.index(marker)


def _expected_labels(name: str) -> list[str]:
    expected_alphabet = _vector(name, "python_alphabet")
    r_alphabet = _vector(name, "r_alphabet")
    r_labels = _vector(name, "labels")
    assert expected_alphabet is not None
    assert r_alphabet is not None
    assert r_labels is not None

    if name in {"sts_replacement_full_meta", "states_mapping"}:
        labels = list(r_labels)
    else:
        # Stable first-observed order is the documented Python adaptation.
        labels = list(expected_alphabet[: len(r_alphabet)])
    labels.extend(expected_alphabet[len(r_alphabet) :])
    return labels


def _expected_palette(name: str) -> list[str]:
    expected_alphabet = _vector(name, "python_alphabet")
    r_alphabet = _vector(name, "r_alphabet")
    r_palette = _vector(name, "cpal")
    assert expected_alphabet is not None
    assert r_alphabet is not None
    assert r_palette is not None

    # R's cpal is aligned with R's alphabet. traminex stores R's order (owner
    # decision D2), and each real state must carry the colour R gives that state.
    palette = [r_palette[r_alphabet.index(state)] for state in expected_alphabet[: len(r_alphabet)]]
    for state in expected_alphabet[len(r_alphabet) :]:
        if state == _scalar(name, "nr"):
            palette.append(_scalar(name, "missing_color"))
        elif state == _scalar(name, "void"):
            palette.append("white")
        else:  # pragma: no cover - malformed golden metadata
            raise AssertionError(f"unexpected internal marker {state!r}")
    return palette


def assert_seqdef_case(actual: STSList, name: str) -> None:
    expected_alphabet = _vector(name, "python_alphabet")
    r_alphabet = _vector(name, "r_alphabet")
    assert expected_alphabet is not None
    assert r_alphabet is not None

    assert _decoded(actual) == _expected_rows(name)
    assert actual.alphabet == expected_alphabet
    assert real_alphabet(actual) == expected_alphabet[: len(r_alphabet)]
    assert actual.nr == _scalar(name, "nr")
    assert actual.void == _scalar(name, "void")
    assert actual.nr_code == _expected_marker_code(name, "nr")
    assert actual.void_code == _expected_marker_code(name, "void")
    assert actual.labels == _expected_labels(name)
    assert actual.cpal == _expected_palette(name)
    assert actual.start == int(_scalar(name, "start"))
    assert actual.xtstep == int(_scalar(name, "xtstep"))
    assert actual.tick_last is (_scalar(name, "tick_last") == "TRUE")
    assert actual.missing_color == _scalar(name, "missing_color")

    expected_ids = _vector(name, "ids")
    assert expected_ids is not None
    assert actual.ids is not None
    assert [str(value) for value in actual.ids] == expected_ids

    expected_weights = _vector(name, "weights")
    if expected_weights is None:
        assert actual.weights is None
    else:
        np.testing.assert_allclose(
            actual.weights, [float(value) for value in expected_weights], rtol=0, atol=R_GOLDEN_ATOL
        )

    # ``cnames`` is explicitly re-expressed as ``xtlab``.  Implicit R source
    # column names remain documented in the cell golden but are not forced
    # into Python plotting metadata, whose default axis derives from ``start``.
    if name == "sts_replacement_full_meta":
        assert actual.xtlab == _vector(name, "xtlab")


def _sts_default_input() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "t1": ["B", None, "A"],
            "t2": [None, "A", None],
            "t3": ["A", None, None],
        },
        index=["default-b", "default-a", "default-tail"],
    )


def _replacement_input() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ignored": ["x", "y", "z"],
            "w1": ["B", "MISS", "A"],
            "w2": ["MISS", "A", "MISS"],
            "w3": ["A", "MISS", "C"],
            "w4": ["MISS", "C", "MISS"],
        },
        index=["source-1", "source-2", "source-3"],
    )


def test_oracle_artifacts_have_the_exact_task13_case_matrix():
    assert set(CELLS["case"]) == EXPECTED_CASES
    for name in EXPECTED_CASES:
        assert {field for case, field in METADATA if case == name} == EXPECTED_METADATA_FIELDS
    assert {case for case, _ in METADATA if case.startswith("palette_")} == {
        "palette_1",
        "palette_2",
        "palette_3",
        "palette_8",
        "palette_9",
        "palette_12",
        "palette_13",
    }
    assert set(ERRORS) == EXPECTED_ERROR_CASES
    assert all(message != "<NO ERROR>" for message in ERRORS.values())


def test_sts_default_and_all_del_preprocessing_match_r():
    default = seqdef(_sts_default_input())
    assert_seqdef_case(default, "sts_default_policies")

    all_del = seqdef(
        _sts_default_input(),
        left="DEL",
        gaps="DEL",
        right="DEL",
        void="~",
        nr="?",
    )
    assert_seqdef_case(all_del, "sts_all_del")


def test_explicit_replacements_and_full_metadata_match_r():
    actual = seqdef(
        _replacement_input(),
        var=["w1", "w2", "w3", "w4"],
        alphabet=["B", "A", "C"],
        states=["Beta", "Alpha", "Gamma"],
        id=["case-z", "case-a", "case-m"],
        weights=[0.5, 1.25, 2.25],
        start=1999,
        left="LEFT",
        gaps="GAP",
        right="RIGHT",
        missing="MISS",
        void="~",
        nr="?",
        cnames=["Y0", "Y1", "Y2", "Y3"],
        xtstep=2,
        tick_last=True,
        cpal=["#110000", "#001100", "#000011", "#aaaa00", "#00aaaa", "#aa00aa"],
        missing_color="#c0ffee",
        labels=["State B", "State A", "State C", "Before", "Between", "After"],
    )
    assert_seqdef_case(actual, "sts_replacement_full_meta")


def test_compressed_sts_and_explicit_states_mapping_match_r():
    compressed = pd.DataFrame(
        {"sequence": ["B;A;B", "A;C;C"]},
        index=["compressed-b", "compressed-a"],
    )
    assert_seqdef_case(seqdef(compressed, stsep=";"), "sts_compressed")

    mapping = pd.DataFrame(
        [["b", "a"], ["a", "b"]],
        index=["mapping-1", "mapping-2"],
        columns=["first", "second"],
    )
    mapped = seqdef(
        mapping,
        alphabet=["b", "a", "unused"],
        states=["Bee", "Ay", "Unused"],
    )
    assert_seqdef_case(mapped, "states_mapping")


def test_sps_and_spell_informat_dispatch_match_r():
    sps = pd.DataFrame(
        [["(A,2)", "(B,1)"], ["(C,1)", "(@,2)"]],
        index=["sps-1", "sps-2"],
        columns=["spell1", "spell2"],
    )
    assert_seqdef_case(
        seqdef(sps, informat="SPS", missing="@"),
        "sps_default",
    )

    spell = pd.DataFrame(
        {
            "who": ["b", "b", "a"],
            "begin": [1, 3, 2],
            "end": [2, 6, 5],
            "state": ["A", "B", "C"],
        }
    )
    with pytest.warns(UserWarning, match="truncated"):
        actual = seqdef(
            spell,
            informat="SPELL",
            process=True,
            limit=4,
            overwrite=True,
            fillblanks="_",
        )
    assert_seqdef_case(actual, "spell_limit4")
    assert actual.issues == {
        "truncated": [1, 2],
        "empty.seq": [],
        "start.before.birth": [],
        "spell.before.birth": [],
    }


@pytest.mark.parametrize("size", [1, 2, 3, 8, 9, 12, 13])
def test_default_palette_matches_r_without_plotting_imports(size: int):
    states = [f"s{i}" for i in range(1, size + 1)]
    actual = seqdef(np.asarray([states], dtype=object), alphabet=states)
    assert actual.cpal == _vector(f"palette_{size}", "cpal")


def _validation_input() -> np.ndarray:
    return np.array([["A", "B"], ["B", "A"]], dtype=object)


ERROR_CALLS: dict[str, Callable[[], object]] = {
    "duplicate_alphabet": lambda: seqdef(_validation_input(), alphabet=["A", "A", "B"]),
    "unknown_state": lambda: seqdef(_validation_input(), alphabet=["A"]),
    "weights_length": lambda: seqdef(_validation_input(), weights=[1]),
    "weights_na": lambda: seqdef(_validation_input(), weights=[1, np.nan]),
    "weights_negative": lambda: seqdef(_validation_input(), weights=[1, -1]),
    "palette_length": lambda: seqdef(_validation_input(), cpal=["red"]),
    "label_length": lambda: seqdef(_validation_input(), labels=["Alpha"]),
    "void_type": lambda: seqdef(_validation_input(), void=1),  # type: ignore[arg-type]
    "void_equals_left": lambda: seqdef(_validation_input(), void="X", left="X"),
    "void_equals_gaps": lambda: seqdef(_validation_input(), void="X", gaps="X"),
    "void_equals_right": lambda: seqdef(_validation_input(), void="X", right="X"),
    "id_length": lambda: seqdef(_validation_input(), id=["only-one"]),
    "id_duplicate": lambda: seqdef(_validation_input(), id=["same", "same"]),
    "id_missing": lambda: seqdef(_validation_input(), id=["one", None]),
}

PYTHON_ERROR_PATTERNS = {
    "duplicate_alphabet": r"duplicat.*alphabet",
    "unknown_state": r"not present in provided alphabet|should be in 'alphabet'",
    "weights_length": r"weights.*length|number of weights",
    "weights_na": r"finite|NA weights",
    "weights_negative": r"non-negative|Negative or NA weights",
    "palette_length": r"cpal|colors",
    "label_length": r"labels.*alphabet|number of labels",
    "void_type": r"void.*character|void must be a string",
    "void_equals_left": r"void.*different|void.*left",
    "void_equals_gaps": r"void.*different|void.*gaps",
    "void_equals_right": r"void.*different|void.*right",
    "id_length": r"ids?.*length|row.names.*length",
    "id_duplicate": r"duplicat.*(?:ids?|row.names)",
    "id_missing": r"missing.*(?:ids?|row.names)|ids?.*missing",
}


@pytest.mark.parametrize("name", sorted(EXPECTED_ERROR_CASES))
def test_seqdef_validation_boundaries_match_r_categories(name: str):
    # Keep the exact executable-R text in the golden even where established
    # Python validation messages are intentionally more actionable.
    assert ERRORS[name]
    with pytest.raises(
        (TypeError, ValueError), match=re.compile(PYTHON_ERROR_PATTERNS[name], re.I)
    ):
        ERROR_CALLS[name]()
