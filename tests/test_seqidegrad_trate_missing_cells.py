"""``method="TRATE"`` does not bridge missing cells in the degradation family.

``seqidegrad.R:233-234`` weighs transitions with ``seqtrate(seqdata)``, which
pairs column ``t`` with column ``t + 1`` only (``seqtrate.R:88-123``); with
``with.missing=FALSE`` a pair touching a missing cell is not counted.  traminex
dropped the missing cells first and then paired the survivors, so ``A * A B``
contributed an ``A -> A`` and an ``A -> B`` transition and the rate matrix
differed from ``seqtrate``'s (2026-09-22 audit S3-03).  The public ``seqtrate``
was already right; the private transition count of ``seqidegrad``,
``seqprecarity`` and ``seqinsecurity`` was not.

Expected values: ``scripts/r_goldens/make_goldens_audit20260922_s3a.R``
(section S3-03), live TraMineR 2.2-14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDENS = Path(__file__).parent / "r_goldens" / "audit20260922_s3a.txt"
ROWS = [["A", None, "A", "B", "C"], ["C", "B", "B", "A", "A"], ["B", "C", "C", "A", "B"]]


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
    return np.array(payload.split(" "), dtype=float)


def _tb(weights=None):
    return tx.seqdef(np.array(ROWS, dtype=object), alphabet=["A", "B", "C"], weights=weights)


@pytest.mark.parametrize(
    "case,fn,kwargs",
    [
        ("s3_03_degrad_trate", tx.seqidegrad, {"method": "TRATE"}),
        ("s3_03_degrad_trate_nointegr", tx.seqidegrad, {"method": "TRATE", "spell_integr": False}),
        ("s3_03_degrad_trateplus", tx.seqidegrad, {"method": "TRATE+"}),
        ("s3_03_degrad_trate_inv", tx.seqidegrad, {"method": "TRATE", "weight_type": "INV"}),
        ("s3_03_degrad_trate_wm", tx.seqidegrad, {"method": "TRATE", "with_missing": True}),
        ("s3_03_prec_trate", tx.seqprecarity, {"method": "TRATE"}),
        ("s3_03_insec_trate", tx.seqinsecurity, {"method": "TRATE"}),
    ],
)
def test_trate_weighting_pairs_adjacent_cells_only(case, fn, kwargs):
    np.testing.assert_allclose(fn(_tb(), **kwargs), _r_values(case), rtol=1e-12, atol=1e-15)


def test_weighted_trate_weighting_pairs_adjacent_cells_only():
    result = tx.seqidegrad(_tb(weights=[1, 2, 0.5]), method="TRATE")
    np.testing.assert_allclose(
        result, _r_values("s3_03_degrad_trate_weighted"), rtol=1e-12, atol=1e-15
    )
