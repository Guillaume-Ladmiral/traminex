"""``seqsurv`` keeps only survival rows with a positive weighted event count.

TraMineRextras ``seqsurv`` returns a ``survival::survfit`` object whose
``summary`` (with its default ``censored = FALSE``) keeps only the times where
``n.event > 0``; with case weights, ``n.event`` is the weighted count, so a time
at which only zero-weight spells end is not a row. traminex emitted a row for
every time at which any spell ended, so zero-weight cases added rows with no
event (7 rows where R has 5), and a time reached only by zero-weight spells,
after every weighted spell had ended, gave a NaN survival row (0/0) (2026-09-22
audit S3-14). Source of truth: ``summary.survfit`` (``index <- fit$n.event >
0``) in survival >= 3.8.9; expected values from
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


def _rows(case: str) -> np.ndarray:
    kind, payload = R[case]
    assert kind == "value", (case, kind, payload)
    if payload == "none":
        return np.empty((0, 7), dtype=float)
    flat = np.array([np.nan if v == "NA" else float(v) for v in payload.split(" ")])
    return flat.reshape(-1, 7)


ROWS = [
    ["A", "A", "B", "B", "C", "C"],
    ["A", "B", "B", "C", "C", "C"],
    ["B", "B", "B", "A", "A", "C"],
    ["C", "C", "A", "A", "B", "B"],
    ["A", "A", "A", "A", "A", "A"],
    ["B", "C", "B", "C", "B", "C"],
]


def _sts(weights, rows=ROWS, alphabet=("A", "B", "C", "D")):
    return tx.seqdef(np.array(rows, dtype=object), alphabet=list(alphabet), weights=weights)


def _check(result, name: str, keys) -> None:
    assert list(result.curves) == list(keys)
    for key in keys:
        np.testing.assert_allclose(
            result.curves[key],
            _rows(f"s314_{name}_{key}"),
            rtol=0,
            atol=1e-10,
            equal_nan=True,
            err_msg=f"curve {key}",
        )


@pytest.mark.parametrize(
    ("name", "weights"),
    [("integer", [0, 0, 1, 1, 0, 2]), ("fractional", [0, 0, 1.5, 2, 0, 0.5])],
)
def test_state_curves_drop_times_without_a_weighted_event(name, weights):
    _check(tx.seqsurv(_sts(weights)), name, ["A", "B", "C"])


def test_group_curves_drop_times_without_a_weighted_event():
    result = tx.seqsurv(
        _sts([0, 0, 1, 1, 0, 2]),
        groups=["g1", "g1", "g2", "g2", "g1", "g2"],
        per_state=True,
        state=["A", "B"],
    )
    _check(result, "per_state", ["g1", "g2"])


def test_a_time_reached_only_by_zero_weight_spells_is_not_a_nan_row():
    x = _sts(
        [1, 0, 1],
        rows=[["A", "B", "B", "B"], ["A", "A", "A", "B"], ["B", "A", "B", "B"]],
        alphabet=("A", "B"),
    )
    result = tx.seqsurv(x)
    _check(result, "after_weighted_end", ["A", "B"])
    assert not np.isnan(result.curves["A"][:, 3]).any()
