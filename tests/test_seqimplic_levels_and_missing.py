"""``seqimplic`` orders group levels and adds the missing state as R does.

TraMineRextras ``seqimplic.R`` (lines 69-94) takes ``levels(factor(group))``:
numbers in numeric order, text in the locale collation, a factor in its declared
order. traminex string-sorted the labels, so ``1, 2, 10`` came out as
``1, 10, 2`` and ``a, A, b`` as ``A, a, b``, which permuted the first axis of
``indices``. With ``with.missing = TRUE`` R always appends the missing state to
the alphabet (label ``"Missing"``, the missing colour, an all-``NA`` slice)
even when no cell is missing; traminex added it only when the data had one
(2026-09-22 audit S3-15). Found while verifying it: with ``na_rm=False`` an NA
group value is still not a level in R, and R's weighted sums over that row are
NA at every position it enters; traminex made ``"None"`` a level and returned
finite values. Expected values:
``scripts/r_goldens/make_goldens_audit20260922_s3b.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
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


def _text(case: str) -> list[str]:
    kind, payload = R[case]
    assert kind == "text", (case, kind, payload)
    return payload.split("|")


def _values(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    return np.array([np.nan if v == "NA" else float(v) for v in payload.split(" ")])


ROWS = [["A", "B", "A"], ["B", "A", "A"], ["A", "A", "B"], ["B", "B", "B"]]
N = None
CASES = {
    "numeric_groups": (ROWS, [1, 2, 10, 2], False),
    "character_case": (ROWS, ["b", "A", "a", "b"], False),
    "declared_factor": (ROWS, pd.Categorical(["x", "y", "x", "y"], categories=["y", "x"]), False),
    "missing_state_absent": (ROWS, [1, 2, 1, 2], True),
    "missing_state_present": (
        [["A", N, "A"], ["B", "A", "A"], ["A", "A", N], ["B", "B", "B"]],
        [1, 2, 1, 2],
        True,
    ),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_seqimplic_matches_r_levels_alphabet_and_indices(name):
    rows, group, with_missing = CASES[name]
    x = tx.seqdef(
        np.array(rows, dtype=object),
        alphabet=["A", "B"],
        right=None,
        cpal=["#7FC97F", "#BEAED4"],
        missing_color="darkgrey",
    )
    result = tx.seqimplic(x, group, with_missing=with_missing)
    prefix = f"s315_{name}"
    assert result.levels == _text(f"{prefix}_levels")
    assert result.alphabet == _text(f"{prefix}_alphabet")
    assert result.labels == _text(f"{prefix}_labels")
    assert result.cpal == _text(f"{prefix}_cpal")
    np.testing.assert_allclose(
        result.indices.flatten(order="F"),
        _values(f"{prefix}_indices"),
        rtol=0,
        atol=1e-10,
        equal_nan=True,
    )


def test_seqimplic_keeps_r_na_semantics_for_an_na_group_when_na_rm_is_false():
    # Found while verifying S3-15: NA is never a factor level, and with
    # na_rm=False R's weighted sums over the NA row are NA at every position
    # that row enters. traminex made "None" a level and returned finite values.
    x = tx.seqdef(np.array(ROWS, dtype=object), alphabet=["A", "B"])
    result = tx.seqimplic(x, [1, None, 2, 2], na_rm=False)
    assert result.levels == _text("s315_na_group_kept_levels")
    np.testing.assert_allclose(
        result.indices.flatten(order="F"),
        _values("s315_na_group_kept_indices"),
        rtol=0,
        atol=1e-10,
        equal_nan=True,
    )
