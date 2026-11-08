"""seqMD_optimize_weights never returns a zero-sum weight vector (TX-AUD-20261011-018).

A traminex extension (no R counterpart). The coordinate search could move every
channel weight to 0, where the combined distance is all zeros, report the result
as converged, and hand back weights that seqMD_distances then refuses ("weights
must have positive sum"); the weights were NumPy floats. A zero-sum vector is
no longer a candidate, and the weights are Python floats.
"""

from __future__ import annotations

import numpy as np

import traminex as tx
from traminex import SeqMD

GROUPS = [0, 0, 0, 1, 1, 1]
# Channel 1 is anti-aligned with the groups; channel 2 is the seed-0 draw of the
# reproduction that made both weights drop to zero.
CHANNEL_1 = tx.seqdef([["A"] * 4, ["B"] * 4, ["C"] * 4] * 2, alphabet=["A", "B", "C"])
CHANNEL_2 = tx.seqdef(
    np.random.default_rng(0).choice(["A", "B", "C"], size=(6, 4)).tolist(), alphabet=["A", "B", "C"]
)


def test_the_weights_have_a_positive_sum_and_are_usable() -> None:
    md = SeqMD([CHANNEL_1, CHANNEL_2])
    result = tx.seqMD_optimize_weights(md, GROUPS, method="HAM")
    assert sum(result["weights"]) > 0
    distances = tx.seqMD_distances(md, method="HAM", weights=result["weights"])
    assert np.asarray(distances).shape == (6, 6)


def test_the_weights_are_python_floats() -> None:
    result = tx.seqMD_optimize_weights(SeqMD([CHANNEL_1, CHANNEL_2]), GROUPS, method="HAM")
    assert all(type(weight) is float for weight in result["weights"])
    assert type(result["score"]) is float
