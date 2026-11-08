"""``seqtransn(pweight=True, with_missing=True)`` is NaN across a missing state.

``seqtransn.R:15-27`` weighs each DSS transition by ``1 - tr[a, b]`` with
``tr <- seqtrate(seqdata)``, which has no row or column for the missing state,
and indexes it with ``seqasnum(dss) + 1``, which is ``NA`` for a missing cell.
A transition to or from the missing state therefore makes that sequence's
weighted count ``NA``.  traminex skipped such transitions and returned a finite
count (``0 0 0`` where R gives ``NA 0 NA``), silently treating them as
weightless (2026-09-22 audit S3-06).  traminex now returns NaN for those
sequences, as R does.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-06), live TraMineR 2.2-14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx
from traminex.sequences import STSList

GOLDENS = Path(__file__).parent / "r_goldens" / "audit20260922_s3a.txt"
ROWS = [
    ["A", None, "B", "B"],
    [None, None, None, None],
    ["B", None, "C", "C"],
    ["A", "B", "B", "C"],
    ["C", "C", None, None],
]


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDENS.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, *_ = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _r_values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", R[case]
    return np.array([np.nan if v == "NA" else float(v) for v in payload.split(" ")])


def _tm() -> STSList:
    return tx.seqdef(np.array(ROWS, dtype=object), alphabet=["A", "B", "C"], right=None)


def _external_twin(x: STSList) -> STSList:
    """The same object with the nr marker declared outside the alphabet."""
    assert x.alphabet == ["A", "B", "C", "*"] and x.nr_code == 3 and x.void_code is None
    return STSList(x.data.copy(), alphabet=["A", "B", "C"], nr="*", void="%", nr_code=3)


@pytest.mark.parametrize("representation", ["in_alphabet", "external"])
@pytest.mark.parametrize(
    "case,kwargs",
    [
        ("s3_06_transn_pweight_wm", {"with_missing": True, "pweight": True}),
        ("s3_06_transn_pweight_wm_norm", {"with_missing": True, "pweight": True, "norm": True}),
        ("s3_06_transn_pweight", {"pweight": True}),
        ("s3_06_transn_wm", {"with_missing": True}),
    ],
)
def test_seqtransn_matches_r_including_na(case, kwargs, representation):
    x = _tm() if representation == "in_alphabet" else _external_twin(_tm())
    np.testing.assert_allclose(
        tx.seqtransn(x, **kwargs), _r_values(case), rtol=0, atol=1e-15, equal_nan=True
    )
