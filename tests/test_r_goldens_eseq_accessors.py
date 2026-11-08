"""R-golden tests for ``seqeid`` / ``seqeweight`` (+ the ``seqeweight<-`` setter).

``seqeid`` returns the per-sequence integer id (the id factor code) and
``seqeweight`` the per-sequence weight (default 1); ``eseq.weights = value`` is
the Pythonic ``seqeweight<-``. Pinned against live TraMineR 2.2-14.

Generator: scripts/r_goldens/make_goldens_eseq_accessors.R.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden_lines

from traminex import STSList, seqecreate, seqeid, seqetm, seqeweight


def _eseq():
    # Character ids whose factor codes (a=1, m=2, z=3) differ from appearance
    # order (z, a, m) -- mirrors make_goldens_eseq_accessors.R.
    return seqecreate(
        pd.DataFrame(
            {
                "id": ["z", "z", "a", "a", "m"],
                "timestamp": [0, 2, 1, 3, 0],
                "event": ["A", "B", "A", "C", "B"],
            }
        )
    )


def _golden_floats(name):
    return np.array([float(v) for v in load_golden_lines(name)], dtype=float)


def test_seqeid_matches_r():
    np.testing.assert_array_equal(seqeid(_eseq()), _golden_floats("eseq_seqeid").astype(int))


def test_seqeweight_default_matches_r():
    np.testing.assert_allclose(
        seqeweight(_eseq()), _golden_floats("eseq_weight_default"), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqeweight_setter_matches_r():
    e = _eseq()
    e.weights = [1.5, 2.0, 0.5]
    np.testing.assert_allclose(
        seqeweight(e), _golden_floats("eseq_weight_set"), rtol=0, atol=R_GOLDEN_ATOL
    )


def test_seqeweight_setter_length_check():
    e = _eseq()
    try:
        e.weights = [1.0, 2.0]  # wrong length (3 sequences)
    except ValueError:
        pass
    else:  # pragma: no cover
        raise AssertionError("expected ValueError for mismatched weights length")


def test_seqeid_scalar_for_single_eseq():
    e = _eseq()
    assert seqeid(e[0]) == 3  # first sequence is id "z" -> factor code 3
    assert seqeweight(e[0]) == 1.0


def test_seqeid_integer_ids_keep_literal_value():
    # R: integer ids go through as.integer() directly (no factor) -> the literal
    # value, NOT a sorted-unique rank. Regression for the string-sort id bug.
    e = seqecreate(
        pd.DataFrame(
            {
                "id": [2, 10, 1, 20, 3],
                "timestamp": [0, 0, 0, 0, 0],
                "event": ["A", "B", "A", "B", "A"],
            }
        )
    )
    np.testing.assert_array_equal(seqeid(e), _golden_floats("eseq_seqeid_int").astype(int))


def test_seqeid_factor_ids_follow_declared_levels():
    # R: a factor id -> as.integer(factor) follows the declared level order. A
    # pandas Categorical with non-alphabetical categories must do the same.
    e = seqecreate(
        pd.DataFrame(
            {
                "id": pd.Categorical(["z", "a", "m"], categories=["z", "a", "m"]),
                "timestamp": [0, 1, 0],
                "event": ["A", "B", "A"],
            }
        )
    )
    np.testing.assert_array_equal(seqeid(e), _golden_floats("eseq_seqeid_factor").astype(int))


def _r3_vectors() -> tuple[list[str], list[int], pd.Categorical]:
    return (
        ["p1", "p1", "p1", "p2", "p2"],
        [0, 2, 5, 0, 4],
        pd.Categorical(
            ["A", "B", "END", "B", "END"],
            categories=["B", "A", "END", "UNUSED"],
        ),
    )


def _r3_sts() -> STSList:
    return STSList(
        data=np.array(
            [
                [0, 0, 1, 2],
                [1, 2, 2, 4],
                [0, 3, 1, 1],
            ],
            dtype=int,
        ),
        alphabet=["A", "B", "C", "*", "%"],
        labels=["Alpha", "Beta", "Gamma", "*", "%"],
        nr="*",
        void="%",
        nr_code=3,
        void_code=4,
        weights=np.array([2.0, 1.0, 3.0]),
    )


def _assert_r3_case(actual, case_name: str) -> None:
    cases = pd.read_csv(R_GOLDENS_DIR / "seqecreate_r3_cases.csv")
    metadata = pd.read_csv(
        R_GOLDENS_DIR / "seqecreate_r3_metadata.csv",
        dtype={"value": str},
        keep_default_na=False,
    )
    expected = cases.loc[cases["case"] == case_name]
    expected_dictionary = metadata.loc[
        (metadata["case"] == case_name) & (metadata["key"] == "dictionary"), "value"
    ].tolist()
    assert actual.as_character() == expected["as_character"].tolist()
    assert actual.dictionary == expected_dictionary
    np.testing.assert_array_equal(seqeid(actual), expected["id_code"].to_numpy(dtype=int))
    np.testing.assert_allclose(actual.lengths, expected["length"], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(actual.weights, expected["weight"], rtol=0, atol=R_GOLDEN_ATOL)


def test_seqecreate_direct_vectors_and_end_event_match_r() -> None:
    ids, timestamps, events = _r3_vectors()
    _assert_r3_case(
        seqecreate(
            id=ids,
            timestamp=timestamps,
            event=events,
            end_event="END",
        ),
        "direct_vectors_end",
    )


@pytest.mark.parametrize("time_name", ["timestamp", "time"])
def test_seqecreate_tse_frame_time_column_routes_match_r(time_name: str) -> None:
    ids, timestamps, events = _r3_vectors()
    frame = pd.DataFrame({"id": ids, time_name: timestamps, "event": events})
    _assert_r3_case(seqecreate(frame), f"frame_{time_name}")


def test_seqecreate_factor_dictionary_and_simultaneous_order_match_r() -> None:
    frame = pd.DataFrame(
        {
            "id": ["p1", "p1", "p1", "p2", "p2"],
            "timestamp": [0, 2, 2, 0, 0],
            "event": pd.Categorical(
                ["Beta", "Gamma", "Alpha", "Gamma", "Alpha"],
                categories=["Gamma", "Beta", "Alpha", "UNUSED"],
            ),
        }
    )
    _assert_r3_case(seqecreate(frame), "factor_simultaneous")


def test_seqecreate_character_mixed_case_dictionary_ids_and_order_match_r() -> None:
    frame = pd.DataFrame(
        {
            "id": ["a", "a", "A", "A"],
            "timestamp": [0, 1, 0, 0],
            "event": ["a", "A", "A", "a"],
        }
    )
    _assert_r3_case(seqecreate(frame), "character_mixed_case")


def test_seqecreate_character_mixed_case_is_hash_seed_independent() -> None:
    script = """
import json
from traminex import seqecreate

result = seqecreate(
    id=["a", "a", "A", "A"],
    timestamp=[0, 1, 0, 0],
    event=["a", "A", "A", "a"],
)
print(json.dumps([result.dictionary, result.as_character()]))
"""
    outputs = []
    for seed in ("1", "2"):
        environment = os.environ.copy()
        environment["PYTHONHASHSEED"] = seed
        completed = subprocess.run(
            [sys.executable, "-c", script],
            check=True,
            capture_output=True,
            text=True,
            env=environment,
        )
        outputs.append(json.loads(completed.stdout))
    assert outputs == [[["a", "A"], ["(a)-1-(A)", "(a,A)"]]] * 2


def test_seqecreate_locale_order_matches_r() -> None:
    """Event labels take R's collation (TX-AUD-20261011-017); they used a portable key."""
    labels = ["_", "-", "!", "?", ".", "*", "%", "a", "A", "á", "ä", "b"]
    metadata = pd.read_csv(
        R_GOLDENS_DIR / "seqecreate_r3_metadata.csv",
        dtype=str,
        keep_default_na=False,
    )
    r_dictionary = metadata.loc[metadata["case"] == "character_locale_order", "value"].tolist()
    assert r_dictionary == labels

    result = seqecreate(id=[1] * len(labels), timestamp=[0] * len(labels), event=labels[::-1])

    assert result.dictionary == r_dictionary
    assert result.as_character() == [f"({','.join(r_dictionary)})"]


@pytest.mark.parametrize(
    ("tevent", "use_labels", "case_name"),
    [
        ("transition", True, "sts_transition_labels"),
        ("state", True, "sts_state_labels"),
        ("period", True, "sts_period_labels"),
        ("state", False, "sts_state_alphabet"),
    ],
)
def test_seqecreate_sts_method_routes_metadata_and_nonmutation_match_r(
    tevent: str, use_labels: bool, case_name: str
) -> None:
    source = _r3_sts()
    data_before = source.data.copy()
    weights_before = source.weights.copy()
    alphabet_before = list(source.alphabet)

    actual = seqecreate(source, tevent=tevent, use_labels=use_labels)

    _assert_r3_case(actual, case_name)
    np.testing.assert_array_equal(source.data, data_before)
    np.testing.assert_array_equal(source.weights, weights_before)
    assert source.alphabet == alphabet_before


def test_seqecreate_sts_weight_switch_and_explicit_matrix_match_r() -> None:
    source = _r3_sts()
    _assert_r3_case(
        seqecreate(source, tevent="state", use_labels=False, weighted=False),
        "sts_unweighted",
    )
    custom = seqetm(source, method="state", use_labels=False).to_numpy(dtype=object)
    custom[0, 1] = "custom1,custom2"
    _assert_r3_case(
        seqecreate(source, tevent=custom, use_labels=False),
        "sts_explicit_matrix",
    )


def test_seqecreate_preserves_custom_tse_column_compatibility() -> None:
    frame = pd.DataFrame({"person": ["p1", "p1", "p2"], "when": [0, 1, 0], "what": ["A", "B", "A"]})
    actual = seqecreate(frame, id_col="person", time_col="when", event_col="what")
    assert actual.as_character() == ["(A)-1-(B)", "(A)"]


def test_seqecreate_r3_errors_and_warning_boundaries_match_r() -> None:
    conditions = (R_GOLDENS_DIR / "seqecreate_r3_conditions.txt").read_text(encoding="utf-8")
    assert "ungrouped_ids|simpleError| [!] Events not grouped by id!" in conditions
    assert "missing_value|simpleError| [!] Missing values not supported" in conditions
    assert "length_mismatch|simpleError|Time ,idpers and event vector" in conditions
    assert "absent_columns|simpleError| [!] Could not find an id argument" in conditions
    assert "invalid_event_characters|warning| [!] some of your events contain" in conditions
    assert "unknown_end_event|success|(A)-1-(B)" in conditions

    with pytest.raises(ValueError, match="Events not grouped by id"):
        seqecreate(id=[1, 2, 1], timestamp=[0, 0, 1], event=["A", "B", "C"])
    with pytest.raises(ValueError, match="Missing values not supported"):
        seqecreate(id=[1, 1], timestamp=[0, np.nan], event=["A", "B"])
    with pytest.raises(ValueError, match="same size"):
        seqecreate(id=[1, 1], timestamp=[0], event=["A", "B"])
    with pytest.raises(ValueError, match="Could not find an id argument"):
        seqecreate(pd.DataFrame({"person": [1], "when": [0], "what": ["A"]}))
    with pytest.warns(UserWarning, match="events contain"):
        seqecreate(id=[1, 1], timestamp=[0, 1], event=["A(B", "C,D"])
    assert seqecreate(
        id=[1, 1], timestamp=[0, 1], event=["A", "B"], end_event="NOPE"
    ).as_character() == ["(A)-1-(B)"]
