"""Exact TraMineR 2.2-14 parity tests for ``seqsep``."""

from __future__ import annotations

import csv
import json
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

import traminex.stats as stats_facade
from traminex.stats_p2 import seqsep

GOLDEN_DIR = Path(__file__).with_name("r_goldens")


def _tsv_rows(name: str) -> list[list[str]]:
    with (GOLDEN_DIR / name).open(encoding="utf-8", newline="") as stream:
        return list(
            csv.reader(
                (line for line in stream if not line.startswith("#")),
                delimiter="\t",
                escapechar="\\",
            )
        )


def _success_cases():
    cases = OrderedDict()
    for case, sl, sep, position, source, expected in _tsv_rows("seqsep_cases.txt"):
        record = cases.setdefault(
            case,
            {"sl": int(sl), "sep": sep, "source": [], "expected": []},
        )
        assert int(position) == len(record["source"]) + 1
        record["source"].append(source)
        record["expected"].append(expected)
    return cases


def _error_cases():
    cases = []
    for case, sl, sep, failing_index, encoded, condition_class, message in _tsv_rows(
        "seqsep_errors.txt"
    ):
        source = []
        if encoded != "[]":
            for expected_position, field in enumerate(encoded.split(" | "), start=1):
                position, value = field.split("=", maxsplit=1)
                assert int(position) == expected_position
                source.append(json.loads(value))
        cases.append(
            {
                "case": case,
                "sl": int(sl),
                "sep": sep,
                "failing_index": int(failing_index),
                "source": source,
                "condition_class": condition_class,
                "message": message,
            }
        )
    return cases


def test_seqsep_list_results_match_r_goldens():
    for case, record in _success_cases().items():
        source = record["source"].copy()

        result = seqsep(source, sl=record["sl"], sep=record["sep"])

        assert isinstance(result, list), case
        assert result == record["expected"], case
        assert source == record["source"], case


def test_seqsep_errors_match_r_and_leave_input_unchanged():
    for record in _error_cases():
        source = record["source"].copy()

        with pytest.raises(ValueError) as exc_info:
            seqsep(source, sl=record["sl"], sep=record["sep"])

        assert record["condition_class"] == "simpleError"
        assert str(exc_info.value) == record["message"], record["case"]
        assert source == record["source"], record["case"]
        if record["case"].startswith("nondivisible"):
            assert str(exc_info.value).endswith(str(record["failing_index"]))


@pytest.mark.parametrize(
    "source",
    [
        np.asarray([], dtype=str),
        pd.Series([], dtype="string", name="trajectory"),
    ],
)
def test_seqsep_empty_array_categories_match_r_empty_vector_error(source):
    record = next(case for case in _error_cases() if case["case"] == "empty_vector")

    with pytest.raises(ValueError) as exc_info:
        seqsep(source, sl=record["sl"], sep=record["sep"])

    assert str(exc_info.value) == record["message"]


def test_seqsep_preserves_one_dimensional_ndarray_category_and_shape():
    source = np.array(["AABB", "CC"], dtype="<U4")

    result = seqsep(source, sl=2, sep="::")

    assert isinstance(result, np.ndarray)
    assert result.shape == source.shape
    assert result.dtype.kind == "U"
    assert result.tolist() == ["AA::BB", "CC"]
    assert source.tolist() == ["AABB", "CC"]


def test_seqsep_preserves_series_index_name_and_string_dtype():
    source = pd.Series(
        ["AABB", "CC"],
        index=pd.Index([101, 205], name="case_id"),
        name="trajectory",
        dtype="string",
    )

    result = seqsep(source, sl=2, sep="|")

    expected = pd.Series(
        ["AA|BB", "CC"],
        index=source.index,
        name=source.name,
        dtype=source.dtype,
    )
    pd.testing.assert_series_equal(result, expected, check_exact=False, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqsep_is_available_from_stats_facade():
    assert stats_facade.seqsep is seqsep


@pytest.mark.parametrize(
    ("source", "error_type", "message"),
    [
        (np.array([["AB", "CD"]]), ValueError, "one-dimensional"),
        ([["AB"], ["CD"]], ValueError, "one-dimensional"),
        (("AB", "CD"), TypeError, "list, one-dimensional numpy.ndarray, or pandas.Series"),
        (["AB", 12], TypeError, "element 2 must be a string"),
    ],
)
def test_seqsep_rejects_unsupported_inputs_actionably(source, error_type, message):
    with pytest.raises(error_type, match=message):
        seqsep(source)
