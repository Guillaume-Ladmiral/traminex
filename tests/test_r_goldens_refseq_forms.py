"""seqdist's refseq forms, as R reads them (TX-AUD-20261011-044).

R takes refseq as one index, a list of two index sets or a one-sequence object;
a vector of several indices, ``c(1, 2)``, is none of these and stops
(seqdist.R:67-99). traminex read ``refseq=[1, 2]`` as two singleton sets and
returned a 1x1 block, and refused an out-of-range or negative index with
IndexError where index sets got ValueError. A flat list, tuple or array of
numbers is now refused as R refuses the vector; two index sets are two
collections (R's ``list(1, 2)`` is ``([1], [2])``); every refusal is ValueError.
A NumPy integer was taken for a sequence object (it has a ``.data`` attribute),
and other types raised TypeError (batch-6 review M-4). Fixture
``scripts/r_goldens/make_goldens_refseq_forms.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "refseq_forms.csv",
    dtype=str,
    keep_default_na=False,
).set_index("case")
S3 = tx.seqdef(
    np.array([["A", "B", "A"], ["B", "B", "A"], ["A", "A", "B"]], dtype=object), alphabet=["A", "B"]
)


def _r_values(case: str) -> np.ndarray:
    assert GOLDEN.loc[case, "ok"] == "TRUE"
    return np.array([float(v) for v in GOLDEN.loc[case, "value"].split()])


@pytest.mark.parametrize("refseq", [[1, 2], (1, 2), np.array([1, 2]), [1]])
def test_a_vector_of_indices_is_refused_as_r_refuses_it(refseq) -> None:
    assert GOLDEN.loc["vector_1_2", "ok"] == "FALSE"
    with pytest.raises(ValueError, match="two collections"):
        tx.seqdist(S3, method="LCS", refseq=refseq)


@pytest.mark.parametrize(
    "case,refseq",
    [
        ("scalar_out_of_range", 9),
        ("scalar_negative", -1),
        ("scalar_fraction", 1.5),
        ("sets_out_of_range", ([1], [9])),
    ],
)
def test_every_refseq_refusal_is_a_value_error(case: str, refseq) -> None:
    assert GOLDEN.loc[case, "ok"] == "FALSE"
    with pytest.raises(ValueError):
        tx.seqdist(S3, method="LCS", refseq=refseq)


@pytest.mark.parametrize(
    "case,refseq",
    [("sets_singletons", ([1], [2])), ("sets_pair_and_one", ([1, 2], [3])), ("scalar_two", 2)],
)
def test_accepted_forms_equal_r(case: str, refseq) -> None:
    actual = np.asarray(tx.seqdist(S3, method="LCS", refseq=refseq), dtype=float)
    np.testing.assert_allclose(actual.ravel(order="F"), _r_values(case), atol=1e-12, rtol=0)


@pytest.mark.parametrize(
    "refseq", [np.int64(2), np.int32(2), np.float64(2.0)], ids=["int64", "int32", "float64"]
)
def test_a_numpy_number_is_one_index(refseq) -> None:
    actual = np.asarray(tx.seqdist(S3, method="LCS", refseq=refseq), dtype=float)
    np.testing.assert_allclose(actual.ravel(order="F"), _r_values("scalar_two"), atol=1e-12, rtol=0)


@pytest.mark.parametrize(
    "refseq",
    ["1", {1}, {"a": 1}, pd.Series([1, 2]), np.array(2)],
    ids=["string", "set", "dict", "series", "zero_dim_array"],
)
def test_other_types_are_refused_with_a_value_error(refseq) -> None:
    assert GOLDEN.loc["character_1", "ok"] == "FALSE"
    with pytest.raises(ValueError, match="refseq"):
        tx.seqdist(S3, method="LCS", refseq=refseq)
