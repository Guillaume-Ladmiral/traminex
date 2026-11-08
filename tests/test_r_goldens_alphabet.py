"""R-golden coverage for the planned ``TraMineR::alphabet`` getter contract.

Generator: ``scripts/r_goldens/make_goldens_alphabet.R`` against TraMineR
2.2-14.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import load_golden_lines

from traminex import STSList, alphabet, seqdef, seqecreate


def _state_sequences() -> tuple[STSList, STSList]:
    no_missing = seqdef([["A", "B"], ["B", "A"]])
    missing_and_void = seqdef([["A", "*", "B"], ["B", "A", "%"]])
    return no_missing, missing_and_void


def _event_sequences():
    return seqecreate(
        pd.DataFrame(
            {
                "id": [1, 1, 2, 2],
                "timestamp": [0, 2, 0, 1],
                "event": ["B", "A", "C", "B"],
            }
        )
    )


def test_alphabet_state_missing_contract_matches_r_without_mutation():
    no_missing, missing_and_void = _state_sequences()
    original = list(no_missing.alphabet)
    original_with_markers = list(missing_and_void.alphabet)
    assert missing_and_void.nr in original_with_markers
    assert missing_and_void.void in original_with_markers

    assert alphabet(no_missing) == load_golden_lines("alphabet_A_wm0")
    with_missing = alphabet(no_missing, with_missing=True)
    assert with_missing == load_golden_lines("alphabet_A_wm1")
    assert with_missing.count(no_missing.nr) == 1
    assert no_missing.void not in with_missing
    assert no_missing.alphabet == original

    ragged_result = alphabet(missing_and_void, with_missing=True)
    assert ragged_result == load_golden_lines("alphabet_B_wm1")
    assert ragged_result.count(missing_and_void.nr) == 1
    assert missing_and_void.void not in ragged_result
    assert missing_and_void.alphabet == original_with_markers


def test_alphabet_returns_independent_state_list():
    no_missing, _ = _state_sequences()
    result = alphabet(no_missing, with_missing=True)
    result.append("changed")
    assert no_missing.alphabet == ["A", "B"]


def test_alphabet_uses_marker_metadata_without_hardcoded_symbols():
    custom = seqdef(
        [["A", "MISSING", "A"], ["A", "A", "VOID"]],
        nr="MISSING",
        void="VOID",
    )
    assert alphabet(custom) == ["A"]
    assert alphabet(custom, with_missing=True) == ["A", "MISSING"]

    marker_text_as_state = STSList(
        np.array([[0, 1]], dtype=int),
        alphabet=["A", "*"],
        nr="*",
    )
    assert alphabet(marker_text_as_state) == ["A", "*"]


def test_alphabet_event_dictionary_matches_r_and_ignores_with_missing():
    events = _event_sequences()
    expected = load_golden_lines("alphabet_E_events")

    assert alphabet(events) == expected
    result = alphabet(events, with_missing=True)
    assert result == expected
    result.append("changed")
    assert events.dictionary == expected


def test_alphabet_rejects_unsupported_input():
    with pytest.raises(TypeError, match="state sequence object or an event sequence object"):
        alphabet(["A", "B"])
