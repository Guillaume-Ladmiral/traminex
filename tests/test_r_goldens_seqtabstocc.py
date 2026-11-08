"""Complete public TraMineRextras 0.6.9 evidence for ``seqtabstocc``."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import STSList, seqdef, seqtabstocc

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"
ATOL = 1e-12


def _combo_fixture(*, weighted: bool = True) -> STSList:
    rows = [
        ["A"],
        ["B"],
        ["C"],
        ["D"],
        ["E"],
        ["A", "B"],
        ["A", "C"],
        ["A", "D"],
        ["A", "E"],
        ["B", "C"],
        ["B", "D"],
        ["B", "E"],
        ["A", "B"],
        ["A", "B"],
        ["A", "C"],
    ]
    matrix = np.full((len(rows), 3), None, dtype=object)
    for position, row in enumerate(rows):
        matrix[position, : len(row)] = row
    return seqdef(
        matrix,
        alphabet=["A", "B", "C", "D", "E"],
        weights=np.arange(1, len(rows) + 1, dtype=float) if weighted else None,
        id=[f"combo-{position}" for position in range(1, len(rows) + 1)],
        cnames=["t1", "t2", "t3"],
    )


def _marker_fixture() -> STSList:
    return seqdef(
        np.array(
            [
                ["A", "*", "B", "*"],
                ["A", "B", "*", "*"],
                ["*", "*", "C", "*"],
                ["C", "*", "*", "*"],
                ["A", "*", "B", "*"],
            ],
            dtype=object,
        ),
        alphabet=["A", "B", "C"],
        missing="*",
        id=[f"marker-{position}" for position in range(1, 6)],
    )


def _locale_fixture() -> STSList:
    states = ["a", "B", "A", "b", "á", "ä", "Z"]
    return seqdef(
        np.column_stack([states, ["X"] * len(states)]),
        alphabet=[*states, "X"],
    )


def _internal_order_fixture() -> STSList:
    return seqdef(
        np.array(
            [["Z", "ä", "A"], ["é", "!", "A"], ["ä", "Z", "A"]],
            dtype=object,
        ),
        alphabet=["Z", "ä", "A", "é", "!"],
    )


def _expected(case: str) -> pd.DataFrame:
    rows = pd.read_csv(
        GOLDEN_DIR / "seqtabstocc_public_cases.csv",
        dtype={"case": str, "call": str, "r_pattern": str},
        keep_default_na=False,
    )
    selected = rows.loc[rows["case"] == case].sort_values("position").reset_index(drop=True)
    assert not selected.empty
    return selected


def _assert_result(actual: pd.DataFrame, case: str) -> None:
    expected = _expected(case)
    assert list(actual.columns) == ["pattern", "count", "prop", "cumprop"]
    assert actual["pattern"].tolist() == expected["r_pattern"].tolist()
    np.testing.assert_allclose(actual["count"], expected["Freq"], rtol=0, atol=ATOL)
    np.testing.assert_allclose(
        actual["prop"],
        expected["Percent"] / 100.0,
        rtol=0,
        atol=ATOL,
    )
    np.testing.assert_allclose(
        actual["cumprop"],
        np.cumsum(expected["Percent"].to_numpy(dtype=float) / 100.0),
        rtol=0,
        atol=ATOL,
    )

    freq = actual.attrs["freq"]
    assert isinstance(freq, pd.DataFrame)
    assert list(freq.columns) == ["Freq", "Percent"]
    assert freq.index.tolist() == expected["r_pattern"].tolist()
    np.testing.assert_allclose(freq["Freq"], expected["Freq"], rtol=0, atol=ATOL)
    np.testing.assert_allclose(freq["Percent"], expected["Percent"], rtol=0, atol=ATOL)
    assert actual.attrs["format"] == "STS"
    assert actual.attrs["weighted"] is False


@pytest.mark.parametrize(
    ("case", "kwargs"),
    [
        ("default_top_ten", {}),
        ("all_patterns", {"idxs": 0}),
        ("reordered_duplicate_ranks", {"idxs": [2, 1, 2]}),
        ("weighted_true_lost", {"idxs": 0, "weighted": True}),
        ("weighted_false_lost", {"idxs": 0, "weighted": False}),
    ],
)
def test_seqtabstocc_rank_and_weight_contract_matches_r(
    case: str,
    kwargs: dict[str, object],
) -> None:
    source = _combo_fixture()
    before = source.copy()

    actual = seqtabstocc(source, **kwargs)

    _assert_result(actual, case)
    assert actual.attrs["nbseq"] == source.nseq
    np.testing.assert_array_equal(source.data, before.data)
    np.testing.assert_array_equal(source.weights, before.weights)
    assert source.ids == before.ids
    assert source.xtlab == before.xtlab
    assert source.alphabet == before.alphabet
    assert source.labels == before.labels
    assert source.cpal == before.cpal


@pytest.mark.parametrize(
    ("case", "include_missing"),
    [
        ("missing_false", False),
        ("missing_true", True),
    ],
)
def test_seqtabstocc_missing_and_void_contract_matches_r(
    case: str,
    include_missing: bool,
) -> None:
    source = _marker_fixture()
    before = source.copy()

    actual = seqtabstocc(source, idxs=0, include_missing=include_missing)

    _assert_result(actual, case)
    assert all("%" not in pattern for pattern in actual["pattern"])
    assert any("*" in pattern for pattern in actual["pattern"]) is include_missing
    np.testing.assert_array_equal(source.data, before.data)
    assert source.nr == before.nr
    assert source.void == before.void
    assert source.ids == before.ids


@pytest.mark.parametrize(
    ("case", "factory"),
    [
        ("locale_ties", _locale_fixture),
        ("internal_state_order", _internal_order_fixture),
    ],
)
def test_seqtabstocc_uses_r_collation_for_state_and_pattern_order(
    case: str,
    factory,
) -> None:
    _assert_result(seqtabstocc(factory(), idxs=0), case)


def test_seqtabstocc_tlim_forwards_with_deprecation_warning() -> None:
    warning_rows = pd.read_csv(
        GOLDEN_DIR / "seqtabstocc_public_warnings.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    expected = warning_rows.loc[warning_rows["case"] == "deprecated_tlim", "message"].iloc[0]
    assert "deprecated" in expected
    with pytest.warns(DeprecationWarning, match="tlim.*deprecated"):
        actual = seqtabstocc(_combo_fixture(), tlim=[2, 1])
    _assert_result(actual, "deprecated_tlim")


def test_seqtabstocc_accepts_with_missing_without_warning() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        actual = seqtabstocc(_marker_fixture(), idxs=0, with_missing=True)
    _assert_result(actual, "missing_true")


@pytest.mark.parametrize(
    ("include_missing", "with_missing"),
    [(False, False), (False, True), (True, False), (True, True)],
)
def test_seqtabstocc_rejects_simultaneous_missing_spellings(
    include_missing: bool,
    with_missing: bool,
) -> None:
    with pytest.raises(ValueError, match="include_missing.*with_missing"):
        seqtabstocc(
            _marker_fixture(),
            include_missing=include_missing,
            with_missing=with_missing,
        )


def test_seqtabstocc_rejects_conflicting_or_unknown_forwarded_options() -> None:
    errors = pd.read_csv(
        GOLDEN_DIR / "seqtabstocc_public_errors.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    assert set(errors["case"]) == {
        "format_conflict",
        "unknown_forwarded",
        "idxs_tlim_conflict",
    }
    with pytest.raises((TypeError, ValueError), match="format"):
        seqtabstocc(_combo_fixture(), format="STS")
    with pytest.raises(TypeError, match="unexpected keyword|banana"):
        seqtabstocc(_combo_fixture(), banana=1)
    with pytest.raises(ValueError, match="cannot be specified together"):
        seqtabstocc(_combo_fixture(), idxs=[2, 1], tlim=1)


def test_seqtabstocc_retains_explicit_python_separator_extension() -> None:
    actual = seqtabstocc(_combo_fixture(), idxs=0, sep="|")
    assert "A|B" in actual["pattern"].tolist()
    assert actual.attrs["python_sep"] == "|"
