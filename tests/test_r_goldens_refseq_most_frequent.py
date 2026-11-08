"""seqdist(refseq=0) picks the most frequent sequence as R does when frequencies tie.

R golden: ``scripts/r_goldens/make_goldens_refseq_most_frequent.R`` (``refseq_most_frequent.csv``).
R ranks tied sequences by the text of ``seqtab(seqnum(seqdata), idxs=1)``: zero-based state
codes, a missing cell written ``*`` whatever the nr spelling, void dropped (seqdist.R:514,
570-574; seqnum.R; seqtab.R:17, 35-45). Audit finding TX-AUD-20261019-001: traminex built the
key from one-based codes, so with ten or more states, or a missing cell, it chose another
reference sequence.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR

import traminex as tx

GOLDEN = pd.read_csv(R_GOLDENS_DIR / "refseq_most_frequent.csv", dtype={"sequence": str})
LETTERS = [chr(code) for code in range(ord("A"), ord("Z") + 1)]
CASES = {
    "eleven_states": {"alphabet": LETTERS[:11], "kwargs": {"method": "LCP"}},
    "missing_tie": {
        "alphabet": LETTERS[:9],
        "kwargs": {"method": "OM", "with_missing": True, "subm": "CONSTANT", "indel": 1.0},
    },
    "missing_tie_nr_letter": {
        "alphabet": LETTERS[:9],
        "nr": "Z",
        "kwargs": {"method": "OM", "with_missing": True, "subm": "CONSTANT", "indel": 1.0},
    },
    "missing_first": {
        "alphabet": LETTERS[:9],
        "kwargs": {"method": "OM", "with_missing": True, "subm": "CONSTANT", "indel": 1.0},
    },
    "weighted_tie": {"alphabet": LETTERS[:11], "kwargs": {"method": "HAM"}},
}
# seqtab.R:17 drops only rows made entirely of nr; a void cell counts as non-nr (batch-7
# review M-1). The golden writes these objects' cells, "*" for nr and "%" for void, which
# seqdef's STS route reads as the markers.
MARKER_ROWS = {"method": "OM", "with_missing": True, "subm": "CONSTANT", "indel": 1.0}
for _name in ("all_void_rows", "nr_then_state", "nr_state_void", "nr_void_row", "all_nr_row"):
    CASES[_name] = {"alphabet": ["A", "B"], "kwargs": MARKER_ROWS}


def _case(name: str):
    rows = GOLDEN[GOLDEN["case"] == name].sort_values("row")
    data = [[None if cell == "NA" else cell for cell in seq.split("-")] for seq in rows["sequence"]]
    return data, rows["weight"].to_numpy(dtype=float), rows["distance"].to_numpy(dtype=float)


@pytest.mark.parametrize("name", sorted(CASES))
@pytest.mark.parametrize("reverse", [False, True], ids=["input_order", "reversed_rows"])
def test_most_frequent_reference_matches_r(name: str, reverse: bool) -> None:
    data, weights, expected = _case(name)
    if reverse:  # R's choice does not depend on the input order of the tied sequences
        data, weights, expected = data[::-1], weights[::-1], expected[::-1]
    spec = CASES[name]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x = tx.seqdef(
            data,
            alphabet=spec["alphabet"],
            weights=None if np.all(weights == 1) else list(weights),
            nr=spec.get("nr", "*"),
        )
        actual = np.asarray(tx.seqdist(x, refseq=0, **spec["kwargs"]), dtype=float).ravel()
    np.testing.assert_allclose(actual, expected, rtol=0, atol=R_GOLDEN_ATOL)
