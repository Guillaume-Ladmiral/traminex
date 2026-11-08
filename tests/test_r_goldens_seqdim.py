"""seqdim follows R's seqfcheck rule (TX-AUD-20261011-071).

traminex inspected only the first element: c("AAB", "A-B-C") gave [1, 2] where R
gives 2 3, a leading NA hid the compressed rows, a single "A-B-C" string raised,
and a 2 x 2 character matrix with a dash was decomposed into [4, 2] where R keeps
2 2. R runs seqfcheck on as.matrix(seqdata) (seqdim.R:5-15, seqfcheck.R:3-23):
numeric data is never decomposed, one-column character data is compressed when
any element holds "-" (else ":"), data with several columns is not decomposed;
then a vector is (1, length) and a matrix (nrow, ncol). Numbers mixed with
text are text, as R's c() makes them: [1, "A-B"] gives 2 2 as c(1, "A-B")
does (batch-5 review M-5). Fixture: ``scripts/r_goldens/make_goldens_seqdim.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import traminex as tx

GOLDEN = pd.read_csv(Path(__file__).resolve().parent / "r_goldens" / "seqdim_cases.csv").set_index(
    "case"
)
# The R inputs of the generator, as Python data (R matrices are column-major).
CASES = {
    "plain_then_compressed": ["AAB", "A-B-C"],
    "na_then_compressed": [None, "A-B-C"],
    "single_compressed": "A-B-C",
    "matrix_with_dash": np.array([["A-B", "D"], ["C", "E"]], dtype=object),
    "one_column_frame": pd.DataFrame({"x": ["A-B", "C-D-E"]}),
    "factor_column": pd.DataFrame({"x": pd.Categorical(["A-B", "C-D"])}),
    "numeric_matrix": np.array([[1, 3, 5], [2, 4, 6]]),
    "all_na_vector": [None, None, None],
    "colon_compressed": ["A:B", "C:D:E"],
    "plain_vector": ["A", "B", "C"],
    "numeric_vector_negative": [-1, 2, 3],
    "plain_scalar": "ABC",
    "mixed_number_text": [1, "A-B"],
    "mixed_negative_number_text": [-1, "A"],
}


def test_every_r_case_has_an_input() -> None:
    assert set(CASES) == set(GOLDEN.index)


@pytest.mark.parametrize("case", sorted(CASES))
def test_seqdim_matches_r(case: str) -> None:
    expected = [int(GOLDEN.loc[case, "rows"]), int(GOLDEN.loc[case, "columns"])]
    assert [int(value) for value in tx.seqdim(CASES[case])] == expected


def test_a_sequence_object_gives_its_shape() -> None:
    seqs = tx.seqdef(np.array([list("AAB"), list("ABB")], dtype=object), alphabet=["A", "B"])
    assert [int(value) for value in tx.seqdim(seqs)] == [2, 3]
