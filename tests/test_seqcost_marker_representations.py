"""``seqcost``/``seqsubm`` give the same costs for both marker representations.

An ``STSList`` whose ``nr_code``/``void_code`` lie outside the alphabet and its
in-alphabet twin hold the same sequences, and the in-alphabet side is what the
``seqcost`` R goldens pin. ``seqcost`` did not canonicalize, so on the external
side the non-time-varying ``INDELS``/``INDELSLOG`` costs failed with a broadcast
``ValueError``, time-varying ``INDELS`` with an ``IndexError``, and ``FUTURE``
with ``with_missing=True`` silently returned different costs when
``void_code == len(alphabet)`` -- the code the widening gives the new missing
state (2026-09-22 audit S1a-02: 90 of 432 configurations differed).
"""

from __future__ import annotations

import itertools
import warnings

import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex import STSList

ROWS = [
    ["A", "B", "B", "C", "C"],
    ["B", "B", "A", "A", None],
    ["C", "A", "A", "B", "B"],
    ["A", "A", "C", "C", "C"],
    ["B", "C", "C", None, None],
]
WEIGHTS = [1.0, 2.0, 0.5, 1.5, 1.0]


def _datasets() -> dict[str, STSList]:
    ragged = tx.seqdef(ROWS, alphabet=["A", "B", "C"], right="DEL", weights=WEIGHTS)
    missing_rows = [list(row) for row in ROWS]
    missing_rows[1][4] = "B"
    missing_rows[4][3:] = ["A", "A"]
    missing_rows[0][2] = None
    missing_rows[3][1] = None
    missing = tx.seqdef(missing_rows, alphabet=["A", "B", "C"], weights=WEIGHTS)
    both_rows = [list(row) for row in ROWS]
    both_rows[0][2] = None
    both = tx.seqdef(both_rows, alphabet=["A", "B", "C"], right="DEL", weights=WEIGHTS)
    return {"ragged": ragged, "missing": missing, "ragged+missing": both}


def _external(twin: STSList, nr_code: int, void_code: int) -> STSList:
    """The same sequences with the markers moved to external codes."""
    markers = {twin.nr_code, twin.void_code}
    real = [code for code in range(len(twin.alphabet)) if code not in markers]
    remap = {old: new for new, old in enumerate(real)}
    if twin.nr_code is not None:
        remap[twin.nr_code] = nr_code
    if twin.void_code is not None:
        remap[twin.void_code] = void_code
    external = STSList(
        data=np.vectorize(remap.get)(twin.data).astype(int),
        alphabet=[twin.alphabet[code] for code in real],
        nr=twin.nr,
        void=twin.void,
        nr_code=nr_code if twin.nr_code is not None else None,
        void_code=void_code if twin.void_code is not None else None,
        weights=None if twin.weights is None else np.asarray(twin.weights).copy(),
    )
    canonical = external._canonicalized()
    assert canonical.alphabet == twin.alphabet
    np.testing.assert_array_equal(canonical.data, twin.data)
    return external


def _costs(result: dict) -> dict[str, np.ndarray]:
    return {
        key: np.asarray(result[key], dtype=float)
        for key in ("indel", "sm", "prox")
        if key in result
    }


FEATURES = pd.DataFrame({"x": [1.0, 2.0, 4.0, 3.0]})


@pytest.mark.parametrize(
    "method", ["CONSTANT", "TRATE", "INDELS", "INDELSLOG", "FUTURE", "FEATURES"]
)
def test_both_marker_representations_give_the_same_costs(method):
    compared = 0
    for (name, twin), codes in itertools.product(_datasets().items(), [(7, 9), (3, 4), (4, 3)]):
        external = _external(twin, *codes)
        for with_missing, time_varying, weighted in itertools.product([False, True], repeat=3):
            if method in {"FUTURE", "FEATURES"} and time_varying:
                continue
            options = {
                "with_missing": with_missing,
                "time_varying": time_varying,
                "weighted": weighted,
            }
            if method == "FEATURES":
                features = FEATURES.iloc[: 3 + int(with_missing)]
                options["state_features"] = features
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                expected = _costs(tx.seqcost(twin, method, **options))
                got = _costs(tx.seqcost(external, method, **options))
            context = f"{name} codes={codes} {options}"
            assert got.keys() == expected.keys(), context
            for key in expected:
                assert got[key].shape == expected[key].shape, f"{context} {key}"
                np.testing.assert_array_equal(got[key], expected[key], err_msg=f"{context} {key}")
            compared += 1
    assert compared >= 36


def test_seqsubm_inherits_the_same_agreement():
    """``seqsubm`` returns ``seqcost(...)["sm"]`` for sequence objects."""
    twin = _datasets()["ragged+missing"]
    external = _external(twin, 4, 3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expected = np.asarray(tx.seqsubm(twin, "INDELS", with_missing=True), dtype=float)
        got = np.asarray(tx.seqsubm(external, "INDELS", with_missing=True), dtype=float)
    np.testing.assert_array_equal(got, expected)
