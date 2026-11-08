"""A given seqrep_args/seqrf_args replaces seqdomassoc's defaults, as in R (TX-AUD-20261011-012).

R's seqdomassoc defaults ``seqrep.args = list(coverage = .8, pradius = .1)`` and
``seqrf.args = list(k = 20)`` (seqdomassoc.R:79-80) are replaced by the user's
list, so seqrep's and seqrf's own defaults apply to what it leaves out. traminex
merged the user's keys into its defaults (``pradius=0.2`` kept coverage 0.8).
Fixture: ``scripts/r_goldens/make_goldens_domassoc_args.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import load_golden

import traminex as tx

GOLDEN = pd.read_csv(
    Path(__file__).resolve().parent / "r_goldens" / "domassoc_args.csv",
    dtype=str,
    keep_default_na=False,
)
ALPHABET = ["1", "2", "3"]
DOMAINS = [
    tx.seqdef(
        np.array(
            [[ALPHABET[c] for c in row] for row in load_golden(f"md_input_s{k}").astype(int)],
            dtype=object,
        ),
        alphabet=ALPHABET,
    )
    for k in (1, 2, 3)
]
DISS = [load_golden(f"md_input_d{k}") for k in (1, 2, 3)]
CASES = {
    "rep_pradius": {"rep_method": "rep", "seqrep_args": {"pradius": 0.2}},
    "rep_default": {"rep_method": "rep"},
    "eq_sortv": {"rep_method": "eq.group", "seqrf_args": {"sortv": "mds", "k": 4}},
    "eq_k_only": {"rep_method": "eq.group", "seqrf_args": {"grp_meth": "first"}},
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_given_arguments_replace_the_defaults_as_in_r(case: str) -> None:
    ours = tx.seqdomassoc(DOMAINS, diss_dom=DISS, **CASES[case])
    rows = GOLDEN[GOLDEN["case"] == case]
    expected = rows.pivot(index="pair", columns="name", values="value").astype(float)
    for column in ["df", "LRT", "v"]:
        np.testing.assert_allclose(
            ours.loc[expected.index, column].to_numpy(dtype=float),
            expected[column].to_numpy(),
            rtol=0,
            atol=1e-10,
        )
