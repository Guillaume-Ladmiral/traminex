"""``seqmodst`` has no modal state where no cell is valid.

``seqmodst.R:28-32`` takes ``which(freq[, i] == max(freq[, i]))[1]``.  Where a
position has no valid cell -- every cell missing (with ``with.missing=FALSE``)
or void, or a zero total weight -- ``seqstatd``'s frequency column is all
``NaN``, the index is ``NA``, and the modal cell and its frequency are ``NA``;
such a modal sequence matches no input sequence.  traminex's ``argmax`` over
the all-NaN column returned the first state, so it reported ``A-A-B-B`` where
R has ``A-NA-B-B`` and ``A-A-A-A`` where R has ``NA-NA-NA-NA`` (2026-09-22
audit S3-07).  An ``STSList`` cannot hold an ``NA`` cell, so traminex marks the
position with the missing state (``nr``), with frequency NaN, and counts no
occurrence, as R does.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-07), live TraMineR 2.2-14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.sequences import STSList

GOLDENS = Path(__file__).parent / "r_goldens" / "audit20260922_s3a.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDENS.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, *_ = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _r_payload(case: str) -> list[str]:
    kind, payload = R[case]
    assert kind == "value", R[case]
    return payload.split(" ")


def _sts(rows, **kwargs) -> STSList:
    return tx.seqdef(np.array(rows, dtype=object), **kwargs)


def _mo() -> STSList:
    rows = [["A", None, "B", "B"], [None, None, None, None], ["B", None, "C", "C"]]
    return _sts(rows, alphabet=["A", "B", "C"], right=None)


def _mo_external() -> STSList:
    x = _mo()
    assert x.alphabet == ["A", "B", "C", "*"] and x.nr_code == 3 and x.void_code is None
    return STSList(x.data.copy(), alphabet=["A", "B", "C"], nr="*", void="%", nr_code=3)


def _mz() -> STSList:
    rows = [["A", "B", "A", "B"], ["B", "B", "A", "A"]]
    return _sts(rows, alphabet=["A", "B"], weights=[0, 0])


def _mr() -> STSList:
    return _sts([["A", "B", None], ["A", None, None]], alphabet=["A", "B"])


CASES = [
    ("missing", _mo, {}),
    ("missing", _mo_external, {}),
    ("missing_wm", _mo, {"include_missing": True}),
    ("missing_wm", _mo_external, {"include_missing": True}),
    ("zero_weights", _mz, {}),
    ("zero_weights_unweighted", _mz, {"weighted": False}),
    ("ragged", _mr, {}),
]


@pytest.mark.parametrize("case,make,kwargs", CASES)
def test_modal_states_match_r_with_missing_for_na(case, make, kwargs):
    x = make()
    modal = tx.seqmodst(x, **kwargs)
    expected = [
        x.nr if state == "NA" else state for state in _r_payload(f"s3_07_modst_{case}_states")
    ]
    assert [modal._state_for_code(int(code)) for code in modal.data[0]] == expected


@pytest.mark.parametrize("case,make,kwargs", CASES)
def test_modal_frequencies_and_occurrences_match_r(case, make, kwargs):
    modal = tx.seqmodst(make(), **kwargs)
    expected = [np.nan if v == "NA" else float(v) for v in _r_payload(f"s3_07_modst_{case}_freq")]
    np.testing.assert_allclose(modal.frequencies_, expected, rtol=0, atol=1e-15, equal_nan=True)
    assert modal.occurrences_ == int(_r_payload(f"s3_07_modst_{case}_occurrences")[0])


def test_a_sequence_missing_at_the_undefined_position_is_not_an_occurrence():
    """Row 1 (``A-*-B-B``) equals the displayed modal sequence cell by cell, but
    R's NA cell matches nothing, so R counts 0 occurrences."""
    modal = tx.seqmodst(_mo())
    assert modal.to_strings() == ["A-*-B-B"]
    assert modal.occurrences_ == 0
