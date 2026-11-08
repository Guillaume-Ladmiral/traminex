"""A one-element ``idxs`` vector is R's scalar in ``seqtab`` (audit S3-21).

R has no scalar type: ``0`` and ``c(0)`` are the same length-one vector.
``seqtab.R`` line 12 refuses ``NULL``, any negative value, and a vector
*longer than one* that contains a value below 1; line 50 then maps a first
element of 0, or any rank beyond the number of distinct sequences, to every
rank. So ``idxs = c(0)`` means all ranks. traminex accepted the Python scalar
``0`` but refused ``[0]`` (and ``np.array([0])``) as a vector containing a
non-positive rank. Expected values:
``scripts/r_goldens/make_goldens_audit20260922_s3c.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import traminex as tx

GOLDEN = Path(__file__).parent / "r_goldens" / "audit20260922_s3c.txt"


def _r() -> dict[str, tuple[str, str]]:
    out = {}
    for line in GOLDEN.read_text().splitlines():
        if line.startswith("#"):
            continue
        case, kind, payload, _warnings = line.split("\t")
        out[case] = (kind, payload)
    return out


R = _r()


def _ranks(case: str) -> list[str]:
    kind, payload = R[case]
    assert kind == "text", (case, kind, payload)
    return payload.split("|")


def _sts() -> tx.STSList:
    rows = [["A", "A", "B", "C"], ["B", "B", "C", "C"], ["A", "B", "A", "B"], ["A", "A", "B", "C"]]
    return tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B", "C"])


@pytest.mark.parametrize(
    ("case", "idxs"),
    [
        ("s321_idxs_scalar_0", 0),
        ("s321_idxs_vector_0", [0]),
        ("s321_idxs_vector_0", (0,)),
        ("s321_idxs_vector_0", np.array([0])),
        ("s321_idxs_vector_0", np.array(0)),
        ("s321_idxs_vector_2", [2]),
        ("s321_idxs_vector_7", [7]),
        ("s321_idxs_vector_2_1", [2, 1]),
        ("s321_idxs_vector_2_9", [2, 9]),
    ],
)
def test_seqtab_rank_selection_matches_r(case, idxs):
    assert tx.seqtab(_sts(), idxs=idxs)["sequence"].tolist() == _ranks(case)


@pytest.mark.parametrize(
    ("case", "idxs"),
    [("s321_idxs_vector_minus_1", [-1]), ("s321_idxs_vector_0_1", [0, 1])],
)
def test_seqtab_still_refuses_what_r_refuses(case, idxs):
    kind, message = R[case]
    assert kind == "error" and "idxs should be" in message
    with pytest.raises(ValueError, match="idxs"):
        tx.seqtab(_sts(), idxs=idxs)
