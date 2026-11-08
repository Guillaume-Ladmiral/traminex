"""``seqstatf`` keeps R's missing-state row and its NaN percentages.

``seqstatf.R:10-22`` sums ``seqistatd(with.missing)``, whose ``nr`` column
exists whenever ``with.missing=TRUE`` (``seqistatd.R:12-14``), and computes
``Percent`` as ``Freq / sum(Freq) * 100``.  traminex added the ``nr`` row only
when the object carried a missing cell, so on missing-free data
``include_missing=True`` returned 3 rows where R returns 4 (with ``*`` at 0),
and with a zero total weight it returned ``Percent`` 0 where R's ``0/0``
gives NaN (2026-09-22 audit S3-09).  ``seqstatf`` follows R's state list
exactly (it is not marker-padded; see the seqmeant/seqstatd DIVERGENCE rows
in ``docs/API_PARITY_DECISIONS.csv``).

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-09), live TraMineR 2.2-14.
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


def _tf() -> STSList:
    rows = [["A", "A", "B"], ["B", "C", "C"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"])


def _tz() -> STSList:
    rows = [["A", "B"], ["B", "B"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B"], weights=[0, 0])


@pytest.mark.parametrize(
    "case,make,kwargs",
    [
        ("wm", _tf, {"include_missing": True}),
        ("zero_weight", _tz, {}),
        ("zero_weight_wm", _tz, {"include_missing": True}),
        ("zero_weight_unweighted", _tz, {"weighted": False}),
    ],
)
def test_seqstatf_matches_r_rows_and_values(case, make, kwargs):
    table = tx.seqstatf(make(), **kwargs)
    assert list(table.index) == _r_payload(f"s3_09_statf_{case}_rownames")
    assert list(table.columns) == ["Freq", "Percent"]
    expected = np.array(_r_payload(f"s3_09_statf_{case}"), dtype=float)
    np.testing.assert_allclose(
        table.to_numpy().ravel(order="F"), expected, rtol=0, atol=1e-12, equal_nan=True
    )


def test_both_marker_representations_give_the_same_table():
    plain = _tf()
    declared_outside = STSList(
        plain.data.copy(), alphabet=["A", "B", "C"], nr="*", nr_code=3, cpal=list(plain.cpal)
    )
    for include_missing in (False, True):
        expected = tx.seqstatf(plain, include_missing=include_missing)
        result = tx.seqstatf(declared_outside, include_missing=include_missing)
        assert list(result.index) == list(expected.index)
        np.testing.assert_allclose(result.to_numpy(), expected.to_numpy(), rtol=0, atol=0)
