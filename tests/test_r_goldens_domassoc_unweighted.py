"""seqdomassoc(weighted=False) selects unweighted representatives (TX-AUD-20261011-011).

With weighted = FALSE, R's seqdomassoc means to pass unit weights to seqrep
(seqdomassoc.R:121-122), but seqrep has no weights argument and drops them
with diss given (seqrep.R:5-24), so its representatives keep the object
weights; for eq.group, seqrf keeps weighted = TRUE (seqdomassoc.R:86). traminex
passes weighted=False to both, a recorded correction. Its values are R's on the
same domains without weights, where R's selection is unweighted. Fixture:
``scripts/r_goldens/make_goldens_domassoc_unweighted.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_DIFFERENCE_ATOL, load_golden

import traminex as tx

HERE = Path(__file__).resolve().parent / "r_goldens"
GOLDEN = pd.read_csv(HERE / "domassoc_unweighted.csv", dtype=str, keep_default_na=False)
WEIGHTS = pd.read_csv(HERE / "domassoc_unweighted_weights.csv")["weight"].to_numpy(dtype=float)
ALPHABET = ["1", "2", "3"]


def _domains():
    domains = []
    for k in (1, 2, 3):
        codes = load_golden(f"md_input_s{k}").astype(int)
        cells = np.array([[ALPHABET[c] for c in row] for row in codes], dtype=object)
        domains.append(tx.seqdef(cells, alphabet=ALPHABET, weights=WEIGHTS))
    return domains


DISS = [load_golden(f"md_input_d{k}") for k in (1, 2, 3)]
OPTIONS = {
    "rep": {"rep_method": "rep", "seqrep_args": {"coverage": 0.8, "pradius": 0.1}},
    "eq.group": {"rep_method": "eq.group", "seqrf_args": {"k": 3}},
}


def _r(method: str, kind: str) -> pd.DataFrame:
    rows = GOLDEN[(GOLDEN["method"] == method) & (GOLDEN["kind"] == kind)]
    table = rows.pivot(index="pair", columns="name", values="value").astype(float)
    return table


@pytest.mark.parametrize("method", sorted(OPTIONS))
def test_weighted_false_gives_r_values_without_weights(method: str) -> None:
    ours = tx.seqdomassoc(_domains(), diss_dom=DISS, weighted=False, **OPTIONS[method])
    expected = _r(method, "unweighted")
    for column in ["df", "LRT", "v"]:
        np.testing.assert_allclose(
            ours.loc[expected.index, column].to_numpy(dtype=float),
            expected[column].to_numpy(),
            rtol=0,
            atol=1e-10,
        )


@pytest.mark.parametrize("method", sorted(OPTIONS))
def test_r_keeps_the_object_weights_with_weighted_false(method: str) -> None:
    assert not np.allclose(
        _r(method, "r_weighted_false")["LRT"].to_numpy(),
        _r(method, "unweighted")["LRT"].to_numpy(),
        rtol=0,
        atol=R_DIFFERENCE_ATOL,
    )
