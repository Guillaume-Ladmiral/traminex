"""R2 contract tests for TraMineR state-order canonicalization."""

from __future__ import annotations

import csv
import importlib
from pathlib import Path

import numpy as np
import pytest

from traminex.format import seqdef
from traminex.stats_complexity import seqibad, seqidegrad, seqprecstart

GOLDENS = Path(__file__).parent / "r_goldens"
SEQS_A = np.asarray(
    [
        ["A", "A", "B", "C", "C", "C"],
        ["A", "B", "B", "B", "C", "A"],
        ["D", "D", "D", "D", "D", "D"],
        ["A", "B", "A", "B", "A", "B"],
    ],
    dtype=object,
)
SEQS_B = np.asarray(
    [
        ["A", "A", "*", "B", "B", "C"],
        ["*", "A", "A", "A", "*", "*"],
        ["A", "*", "A", "B", "B", "B"],
        ["C", "C", "C", "C", "C", "C"],
    ],
    dtype=object,
)


def _fixture(tag: str):
    if tag == "A":
        return seqdef(SEQS_A, alphabet=["A", "B", "C", "D"])
    return seqdef(SEQS_B, alphabet=["A", "B", "C"], nr="*")


def _rows() -> list[dict[str, str]]:
    with (GOLDENS / "state_order_r2_cases.csv").open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _module():
    try:
        return importlib.import_module("traminex._state_order")
    except ModuleNotFoundError:
        pytest.fail("R2 state-order module has not been implemented")


def _split(value: str) -> list[str]:
    return [] if value == "" else value.split("|")


def _numbers(value: str) -> np.ndarray:
    return np.asarray(
        [np.nan if item in {"NA", "NaN"} else float(item) for item in _split(value)],
        dtype=float,
    )


def _groups(value: str) -> tuple[tuple[int, ...], ...]:
    if value == "":
        return ()
    return tuple(tuple(int(item) for item in group.split("+")) for group in value.split(";"))


def _case_kwargs(case: str) -> dict[str, object]:
    return {
        "full_default": {},
        "full_reverse": {"state_order": ["D", "C", "B", "A"]},
        "partial_omit_leading": {"state_order": ["B", "C", "D"]},
        "partial_omit_trailing": {"state_order": ["A", "B", "C"]},
        "partial_omit_middle": {"state_order": ["A", "D"]},
        "partial_singleton": {"state_order": ["A"]},
        "equiv_singleton": {"state_equiv": [["B"]]},
        "equiv_multi": {"state_equiv": [["B", "C"]]},
        "partial_promote_one": {
            "state_order": ["A", "D"],
            "state_equiv": [["D", "B"]],
        },
        "partial_promote_multi": {
            "state_order": ["A", "D"],
            "state_equiv": [["D", "B", "C"]],
        },
        "partial_equiv_noncomp": {
            "state_order": ["A"],
            "state_equiv": [["B", "C"]],
        },
        "explicit_tied_negative": {"stprec": [0.0, 0.5, 0.5, -1.0]},
        "explicit_negative_surviving_equiv": {
            "stprec": [0.0, 0.2, 0.8, -1.0],
            "state_equiv": [["B", "C"]],
        },
        "explicit_constant": {"stprec": [1.0, 1.0, 1.0, 1.0]},
        "explicit_all_negative": {"stprec": [-1.0, -1.0, -1.0, -1.0]},
        "missing_observed": {},
        "missing_observed_equiv": {
            "state_order": ["A", "C"],
            "state_equiv": [["C", "*"]],
        },
        "phantom_missing": {},
        "phantom_missing_equiv": {"state_equiv": [["*", "A"]]},
        "phantom_explicit_negative": {"stprec": [0.0, 0.2, 0.4, 0.6, -1.0]},
    }[case]


def test_resolver_matches_complete_r_canonicalization_manifest():
    module = _module()
    assert len(_rows()) == 20
    for row in _rows():
        kwargs = _case_kwargs(row["case"])
        contract = module._resolve_state_order(
            _fixture(row["fixture"]),
            with_missing=row["with_missing"] == "TRUE",
            **kwargs,
        )
        assert isinstance(contract, module.StateOrderContract)
        assert contract.labels == tuple(_split(row["labels"]))
        assert contract.groups == _groups(row["groups"])
        np.testing.assert_allclose(
            contract.precedence,
            _numbers(row["precedence"]),
            rtol=0,
            atol=1e-14,
            equal_nan=True,
            err_msg=row["case"],
        )
        np.testing.assert_array_equal(contract.code_to_group, _numbers(row["code_to_group"]))


def test_contract_arrays_are_defensive_read_only_copies():
    module = _module()
    state_order = ["A", "D"]
    state_equiv = [["D", "B", "C"]]
    stprec = np.asarray([0.0, 0.2, 0.8, -1.0])
    contract = module._resolve_state_order(
        _fixture("A"),
        state_order=state_order,
        state_equiv=state_equiv,
        stprec=stprec,
    )
    expected_precedence = contract.precedence.copy()
    expected_mapping = contract.code_to_group.copy()

    state_order[:] = ["D"]
    state_equiv[0][:] = ["A"]
    stprec[:] = 99
    np.testing.assert_array_equal(contract.precedence, expected_precedence)
    np.testing.assert_array_equal(contract.code_to_group, expected_mapping)
    with pytest.raises(ValueError, match="read-only"):
        contract.precedence[0] = 1
    with pytest.raises(ValueError, match="read-only"):
        contract.code_to_group[0] = 1


@pytest.mark.parametrize(
    ("kwargs", "error", "match"),
    [
        ({"state_order": ["A", "A"]}, ValueError, "Multiple occurrences"),
        ({"state_order": ["A", "Z"]}, ValueError, "not in the alphabet"),
        ({"state_equiv": ["A", "B"]}, TypeError, "list of equivalence classes"),
        ({"state_equiv": [["B", "B"]]}, ValueError, "Multiple occurrence"),
        (
            {"state_equiv": [["A", "B"], ["B", "C"]]},
            ValueError,
            "Multiple occurrence",
        ),
        ({"state_equiv": [["A", "Z"]]}, ValueError, "not in the alphabet"),
        ({"stprec": [0.0, 1.0]}, ValueError, "length"),
        ({"stprec": [0.0, 1.0, np.nan, 2.0]}, ValueError, "NA/NaN"),
        ({"state_order": []}, ValueError, "wrong sign"),
    ],
)
def test_resolver_validation_matches_r_dispositions(kwargs, error, match):
    module = _module()
    with pytest.raises(error, match=match):
        module._resolve_state_order(_fixture("A"), **kwargs)


def test_direct_seqprecstart_seqibad_and_seqidegrad_match_manifest():
    for row in _rows():
        x = _fixture(row["fixture"])
        kwargs = _case_kwargs(row["case"])
        with_missing = row["with_missing"] == "TRUE"
        np.testing.assert_allclose(
            seqprecstart(x, include_missing=with_missing, **kwargs),
            _numbers(row["precedence"]),
            rtol=0,
            atol=1e-14,
            equal_nan=True,
            err_msg=f"seqprecstart:{row['case']}",
        )
        np.testing.assert_allclose(
            seqibad(x, include_missing=with_missing, **kwargs),
            _numbers(row["seqibad"]),
            rtol=0,
            atol=1e-14,
            equal_nan=True,
            err_msg=f"seqibad:{row['case']}",
        )
        np.testing.assert_allclose(
            seqidegrad(x, with_missing=with_missing, method="RANK", spell_integr=False, **kwargs),
            _numbers(row["seqidegrad"]),
            rtol=0,
            atol=1e-14,
            equal_nan=True,
            err_msg=f"seqidegrad:{row['case']}",
        )


def test_r_error_manifest_is_complete_and_deterministic():
    lines = (GOLDENS / "state_order_r2_errors.txt").read_text(encoding="utf-8").splitlines()
    assert lines[0] == "case=message"
    assert [line.split("=", 1)[0] for line in lines[1:]] == [
        "duplicate_order",
        "unknown_order",
        "equiv_not_list",
        "duplicate_equiv",
        "overlap_equiv",
        "unknown_equiv",
        "short_stprec",
        "nan_stprec",
        "empty_order",
        "null_order",
    ]
