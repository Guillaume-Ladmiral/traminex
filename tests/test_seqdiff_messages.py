"""seqdiff lets no seqdist notice reach the caller, as R's seqdiff does (TX-AUD-20261019-059).

R runs seqdist in every window inside ``suppressMessages()`` (seqdiff.R:62), so seqdist's msg.warn notices
("with.missing set as FALSE", "Empty sequence(s)") never reach the caller. traminex surfaced the with.missing notice
once per window. ``scripts/r_goldens/make_goldens_seqdiff_messages.R`` records R's conditions (0 in both cases) and
results for data without and with missing cells (Review Focus 3: the results stay R's).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "seqdiff_messages.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
GROUP = np.array(["m", "f", "m", "f", "m", "f"])
DATA = {
    "missing_free": [
        list("AABBC"),
        list("ABBCC"),
        list("CCBAA"),
        list("BBCCA"),
        list("AAAAB"),
        list("CCCBB"),
    ],
    "with_missing_cells": [
        ["A", None, "B", "B", "C"],
        list("ABBCC"),
        ["C", "C", None, "A", "A"],
        list("BBCCA"),
        list("AAAAB"),
        list("CCCBB"),
    ],
}


def _values(text: str) -> np.ndarray:
    return np.array([float(v) for v in text.split()])


#: The generator's calls: data and seqdiff options. "chi2_breaks_step" adds seqdist's notice that step is ignored
#: when breaks is given, R's msg.warn.ign2 (seqdist.R:258; batch-9 review M-2).
CASES = {
    "missing_free": ("missing_free", dict(with_missing=True)),
    "with_missing_cells": ("with_missing_cells", dict(with_missing=True)),
    "chi2_breaks_step": (
        "missing_free",
        dict(seqdist_args=dict(method="CHI2", breaks=[(1, 2)], step=2)),
    ),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_no_notice_reaches_the_caller_and_results_are_r_s(case: str) -> None:
    assert int(GOLDEN.loc[case, "conditions"]) == 0
    data, options = CASES[case]
    x = tx.seqdef(np.array(DATA[data], dtype=object))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = tx.seqdiff(x, GROUP, cmprange=(0, 1), **options)
    assert [str(w.message) for w in caught] == []
    np.testing.assert_allclose(
        result["stat"].to_numpy().ravel(order="F"),
        _values(GOLDEN.loc[case, "stat"]),
        rtol=0,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        result["discrepancy"].to_numpy().ravel(order="F"),
        _values(GOLDEN.loc[case, "discrepancy"]),
        rtol=0,
        atol=1e-10,
    )
