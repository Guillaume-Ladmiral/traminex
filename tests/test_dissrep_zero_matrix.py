"""``dissrep`` on an all-zero dissimilarity matrix follows TraMineR.

Audit 2026-09-22, finding S4-08. When every dissimilarity is zero (all
sequences identical), traminex reported ``nb = 0`` and quality statistics of
0. TraMineR 2.2-14 sets ``nb`` to the total weight in that case (the ``diss0``
rule, ``SOURCES/TraMineR/R/dissrep.R:16`` and ``:178``) and divides plainly, so
``Q = (DC - SD) / DC * 100`` and the overall quality are ``0/0 = NaN``
(``dissrep.R:188``, ``:200``) -- there is no discrepancy to explain. Expected
values come from live R through
``scripts/r_goldens/make_goldens_audit20260922_s4.R``.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from traminex import dissrep

FIXTURE = Path(__file__).with_name("r_goldens") / "dissrep_zero_matrix.txt"
COLUMNS = ["na", "na(%)", "nb", "nb(%)", "SD", "MD", "DC", "V", "Q"]


def _cases() -> dict[str, dict[str, str]]:
    cases: dict[str, dict[str, str]] = {}
    for line in FIXTURE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        case, key, value = line.split("\t")
        cases.setdefault(case, {})[key] = value
    return cases


CASES = _cases()


def _optional(text: str):
    return None if text == "none" else text


@pytest.mark.parametrize("case", sorted(CASES))
def test_zero_matrix_statistics_match_r(case):
    row = CASES[case]
    n = int(row["n"])
    weights = _optional(row["weights"])
    nrep = _optional(row["nrep"])
    dmax = _optional(row["dmax"])

    result = dissrep(
        np.zeros((n, n)),
        weights=None if weights is None else [float(v) for v in weights.split(" ")],
        nrep=None if nrep is None else int(nrep),
        dmax=None if dmax is None else float(dmax),
    )

    expected_reps = [int(v) - 1 for v in row["representatives"].split(" ")]
    assert result["representatives"].tolist() == expected_reps
    values = np.array([float(v) for v in row["statistics"].split(" ")])
    expected = values.reshape(len(COLUMNS), -1).T  # R wrote the data frame column-major
    np.testing.assert_array_equal(result["statistics"][COLUMNS].to_numpy(), expected)
    assert math.isnan(float(row["quality"]))
    assert math.isnan(result["quality"])
