"""``dissrep``'s mean distance for a representative with no weight, as in R.

Audit 2026-09-22, finding S4-17. When every case assigned to a representative
has weight 0, traminex reported its mean distance ``MD`` as 0. TraMineR 2.2-14
computes ``MD <- SD/na`` (``SOURCES/TraMineR/R/dissrep.R:151-153``), a plain
division, so that representative's ``MD`` is ``0/0 = NaN``: there is no weight
over which to average. Expected values come from live R through
``scripts/r_goldens/make_goldens_audit20260922_s4b.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import dissrep

FIXTURE = Path(__file__).with_name("r_goldens") / "dissrep_zero_weight_representative.txt"
COLUMNS = ["na", "na(%)", "nb", "nb(%)", "SD", "MD", "DC", "V", "Q"]


def _rows() -> dict[str, np.ndarray]:
    rows: dict[str, np.ndarray] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        key, value = line.split("\t")
        rows[key] = np.array([np.nan if v == "NA" else float(v) for v in value.split(" ")])
    return rows


def test_zero_weight_representative_matches_r():
    rows = _rows()
    positions = rows["positions"]
    D = np.abs(positions[:, None] - positions[None, :])

    result = dissrep(D, weights=rows["weights"], nrep=2)

    assert (result["representatives"] + 1).tolist() == rows["representatives"].astype(int).tolist()
    expected = rows["statistics"].reshape(len(COLUMNS), -1).T  # R: column-major
    assert np.isnan(expected[1, COLUMNS.index("MD")])
    np.testing.assert_allclose(
        result["statistics"][COLUMNS].to_numpy(), expected, rtol=1e-12, atol=R_GOLDEN_ATOL
    )
    np.testing.assert_allclose(result["quality"], rows["quality"][0], rtol=1e-12)
