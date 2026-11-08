"""``seqrep_grp`` leaves out cases whose group is missing, as TraMineRextras does.

Audit 2026-09-22, finding S4-10. A missing group value made ``seqrep_grp``
fail: ``None`` raised ``TypeError`` while sorting the levels, and ``NaN`` (never
equal to itself) became a one-case level and triggered "Each group must have 2
or more cases". TraMineRextras 0.6.9 ``seqrep.grp``
(``SOURCES/TraMineRextras/R/seqrep-grp.R:15-45``) keeps the case's ``NA``
through ``TraMineR:::group``, and ``xtabs(~grp)``, ``levels(factor(grp))`` and
``which(grp == level)`` all leave it out, so the case simply belongs to no
group. Expected statistics come from live R through
``scripts/r_goldens/make_goldens_audit20260922_s4.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL

from traminex import seqdef
from traminex.cluster import seqrep_grp

FIXTURE = Path(__file__).with_name("r_goldens") / "seqrep_grp_missing_group.txt"
STAT_COLUMNS = ["na", "na(%)", "nb", "nb(%)", "SD", "MD", "DC", "V", "Q"]


def _rows() -> dict[str, dict[str, list[str]]]:
    cases: dict[str, dict[str, list[str]]] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        case, key, value = line.split("\t")
        cases.setdefault(case, {}).setdefault(key, []).append(value)
    return cases


CASES = _rows()
SEQUENCES = seqdef(
    np.array([list(row) for row in CASES["input"]["rows"]], dtype=object),
    alphabet=["A", "B", "C"],
)
DISS = np.array([float(v) for v in CASES["input"]["diss"][0].split(" ")]).reshape(8, 8)

GROUPS = {
    "character-None": ("character", np.array(["x", "x", "x", "y", "y", None, "y", "x"], object)),
    "character-pd.NA": (
        "character",
        pd.array(["x", "x", "x", "y", "y", pd.NA, "y", "x"], dtype="string").to_numpy(),
    ),
    "numeric-NaN": ("numeric", np.array([1, 1, 1, 2, 2, np.nan, 2, 1], dtype=float)),
}


@pytest.mark.parametrize("label", sorted(GROUPS))
def test_missing_group_cases_are_left_out(label):
    case, group = GROUPS[label]
    expected = {
        key[len("stat_") :]: np.array([float(v) for v in values[0].split(" ")])
        for key, values in CASES[case].items()
        if key.startswith("stat_")
    }

    stats = seqrep_grp(SEQUENCES, group, diss=DISS, ret="stat")

    # Levels are compared in sorted order, as in test_r_goldens_seqrep_grp.py
    # (R names them " x", " 1", ...; a NaN-coded group is a float array here).
    assert len(stats) == len(expected)
    for got_key, want_key in zip(sorted(stats), sorted(expected)):
        want = expected[want_key].reshape(len(STAT_COLUMNS), -1).T  # R: column-major
        got = stats[got_key][STAT_COLUMNS].to_numpy(dtype=float)
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=R_GOLDEN_ATOL)
