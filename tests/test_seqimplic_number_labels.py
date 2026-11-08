"""``seqimplic`` labels numeric and logical groups as R's ``factor()`` does.

TraMineRextras ``seqimplic.R`` takes ``levels(factor(group))``, which are
``as.character`` of the sorted unique values: doubles with 15 significant
digits (``2``, ``10``, ``1e+05``; ``0.1 + 0.2`` and ``0.3`` both become the
level ``0.3``) and logicals as ``FALSE``/``TRUE``. traminex used Python's
``str``, so its levels read ``2.0``, ``10.0``, ``100000.0``,
``0.30000000000000004`` or ``False``/``True`` (2026-09-22 audit S3-20). The
package's ``r_as_character`` already renders numbers as R does. Expected
values: ``scripts/r_goldens/make_goldens_audit20260922_s3c.R``.
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


def _text(case: str) -> list[str]:
    kind, payload = R[case]
    assert kind == "text", (case, kind, payload)
    return payload.split("|")


def _values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    return np.array([np.nan if v == "NA" else float(v) for v in payload.split(" ")])


GROUPS = {
    "fractional": [1.5, 2.0, 10.0, 2.0],
    "whole_doubles": [2.0, 10.0, 2.0, 10.0],
    "large": [1e5, 2.0, 1e5, 2.0],
    "float_neighbours": [0.1 + 0.2, 0.3, 1.0, 1.0],
    "logical": [True, False, True, False],
}


@pytest.mark.parametrize("name", sorted(GROUPS))
def test_seqimplic_levels_and_indices_match_r(name):
    rows = [["A", "B", "A"], ["B", "A", "A"], ["A", "A", "B"], ["B", "B", "B"]]
    x = tx.seqdef(np.array(rows, dtype=object), alphabet=["A", "B"])
    result = tx.seqimplic(x, GROUPS[name])
    assert result.levels == _text(f"s320_{name}_levels")
    np.testing.assert_allclose(
        result.indices.flatten(order="F"),
        _values(f"s320_{name}_indices"),
        rtol=0,
        atol=1e-10,
        equal_nan=True,
    )
