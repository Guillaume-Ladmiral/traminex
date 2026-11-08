"""seqpolyads gives a two-sided p-value for a negative U, where R exceeds 1 (TX-AUD-20261010-052).

R computes ``U.tp <- 2*pt(mean.U/se, T-1, lower.tail=F)`` (seqpolyads.R:184):
the two-sided p-value its factor 2 intends only when U is positive; a negative
U gives ``U.tp`` above 1, which no probability can be. traminex computes the
two-sided p-value ``2 * P(T > |U| / se)`` for either sign, a recorded
correction (ledger row seqpolyads). traminex draws its own random pairs, so its
U differs from R's; the test checks the definition on traminex's draws and
that R's recorded value is not a probability. This behaviour predates the
record, so the test characterizes it. Fixture:
``scripts/r_goldens/make_goldens_seqpolyads_negative_u.R``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

import traminex as tx

GOLDEN = pd.read_csv(Path(__file__).resolve().parent / "r_goldens" / "seqpolyads_negative_u.csv")
FIRST = tx.seqdef(np.array([["A"] * 3, ["B"] * 3] * 3, dtype=object), alphabet=["A", "B"])
SECOND = tx.seqdef(np.array([["B"] * 3, ["A"] * 3] * 3, dtype=object), alphabet=["A", "B"])


def test_r_gives_a_p_value_above_one_for_a_negative_u() -> None:
    assert (GOLDEN["U"] < 0).all()
    assert (GOLDEN["U.tp"] > 1).all()


def test_traminex_gives_the_two_sided_p_value() -> None:
    result = tx.seqpolyads([FIRST, SECOND], T=200, s=7)
    u = np.asarray(result["U"], dtype=float)
    assert (u < 0).all()
    se = np.std(result["random.dist"], ddof=1) / np.sqrt(200)
    expected = 2 * student_t.sf(np.abs(u) / se, 199)
    np.testing.assert_allclose(result["U.tp"], expected, rtol=0, atol=1e-15)
    assert (np.asarray(result["U.tp"]) <= 1).all()
