"""The ``weighted`` flag of ``seqmodst`` and ``seqmeant`` follows TraMineR.

TraMineR's ``seqmodst.R`` (lines 57-66) and ``seqmeant.R`` (lines 12-20, 49)
record ``weighted = FALSE`` when the object has no weights, when
``weighted = FALSE`` is requested, or when every weight equals 1 (the form old
TraMineR versions stored for unweighted objects). traminex's ``seqstatd``
already did so, but ``seqmodst().weighted_`` echoed the argument (True without
weights) and ``seqmeant().attrs["weighted"]`` was True for all-one weights
(2026-09-22 audit S3-16). Expected values:
``scripts/r_goldens/make_goldens_audit20260922_s3b.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s3b.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, _warnings = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _flag(case: str) -> bool:
    kind, payload = R[case]
    assert kind == "text" and payload in {"TRUE", "FALSE"}, (case, kind, payload)
    return payload == "TRUE"


def _number(case: str) -> float:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    return float(payload)


ROWS = [["A", "A", "B"], ["B", "B", "A"], ["A", "B", "B"]]
WEIGHTS = {"no_weights": None, "unit_weights": [1, 1, 1], "other_weights": [1, 2, 0.5]}


def _sts(name: str) -> tx.STSList:
    return tx.seqdef(np.array(ROWS, dtype=object), alphabet=["A", "B"], weights=WEIGHTS[name])


@pytest.mark.parametrize("name", sorted(WEIGHTS))
@pytest.mark.parametrize("weighted", [True, False])
def test_seqmodst_weighted_flag_and_nbseq_match_r(name, weighted):
    modal = tx.seqmodst(_sts(name), weighted=weighted)
    prefix = f"s316_{name}_weighted{int(weighted)}"
    assert modal.weighted_ is _flag(f"{prefix}_seqmodst_weighted")
    assert modal.nbseq_ == _number(f"{prefix}_seqmodst_nbseq")


@pytest.mark.parametrize("name", sorted(WEIGHTS))
@pytest.mark.parametrize("weighted", [True, False])
def test_seqmeant_weighted_flag_and_nbseq_match_r(name, weighted):
    mean_times = tx.seqmeant(_sts(name), weighted=weighted)
    prefix = f"s316_{name}_weighted{int(weighted)}"
    assert mean_times.attrs["weighted"] is _flag(f"{prefix}_seqmeant_weighted")
    assert mean_times.attrs["nbseq"] == _number(f"{prefix}_seqmeant_nbseq")
