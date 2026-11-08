"""seqentrans' Avg.occ is the average number of occurrences by sequence (TX-AUD-20261019-016).

seqentrans.Rd: Avg.occ is the average number of occurrences by sequence. R divides Count by the
total weight of the sequences whatever weighting produced Count (seqentrans.R:11-14), so after
mining with ``weighted=False`` on sequences of weight 10 it reports a tenth of the average.
traminex divides by the total the counts were taken on (DIVERGENCE row of seqentrans); R's
values are recorded by ``scripts/r_goldens/make_goldens_seqentrans.R``
(``seqentrans_unweighted.csv``). R's seqentrans returns a modified copy and leaves its argument
unchanged; traminex does the same.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDENS_DIR

import traminex as tx

GOLDEN = pd.read_csv(R_GOLDENS_DIR / "seqentrans_unweighted.csv")


def _mined(weighted: bool):
    sequences = tx.seqecreate(
        id=[1, 1, 2, 2, 3], timestamp=[0, 1, 0, 2, 0], event=["A", "B", "A", "B", "A"]
    )
    sequences.weights = [10, 10, 10]
    return tx.seqefsub(sequences, pmin_support=0.1, weighted=weighted)


def _r(weighted: bool) -> pd.DataFrame:
    return GOLDEN[GOLDEN["weighted"] == weighted].reset_index(drop=True)


def test_unweighted_mining_averages_over_the_sequences() -> None:
    # (A) occurs in 3 of the 3 sequences, (A)-(B) and (B) in 2: averages 1, 2/3, 2/3.
    result = tx.seqentrans(_mined(False), avg_occ=True)
    r_values = _r(False)
    assert list(result.subseq) == r_values["subseq"].tolist()
    np.testing.assert_allclose(result.data["Count"], r_values["Count"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.data["Avg.occ"], [1.0, 2 / 3, 2 / 3], rtol=0, atol=1e-12)
    np.testing.assert_allclose(r_values["Avg.occ"], [0.1, 1 / 15, 1 / 15], rtol=0, atol=1e-12)


def test_weighted_mining_is_r_s() -> None:
    result = tx.seqentrans(_mined(True), avg_occ=True)
    r_values = _r(True)
    np.testing.assert_allclose(result.data["Count"], r_values["Count"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.data["Avg.occ"], r_values["Avg.occ"], rtol=0, atol=1e-12)


@pytest.mark.parametrize("weighted", [False, True])
def test_a_subset_keeps_the_total_its_counts_were_taken_on(weighted: bool) -> None:
    full = tx.seqentrans(_mined(weighted), avg_occ=True)
    subset = tx.seqentrans(_mined(weighted)[[0, 2]], avg_occ=True)
    np.testing.assert_allclose(
        subset.data["Avg.occ"], full.data["Avg.occ"].iloc[[0, 2]], rtol=0, atol=1e-12
    )


def test_the_argument_is_left_unchanged() -> None:
    mined = _mined(True)
    columns = list(mined.data.columns)
    result = tx.seqentrans(mined, avg_occ=True)
    assert result is not mined
    assert list(mined.data.columns) == columns
    assert {"ntrans", "nevent", "Avg.occ"} <= set(result.data.columns)


def test_dataframe_route_error_says_who_sets_the_total() -> None:
    """No traminex function sets ``seqe_weight_total``; the message said seqefsub did (batch-8 review M-2)."""
    frame = pd.DataFrame({"subseq": ["(A)-(B)"], "count": [2]})
    with pytest.raises(ValueError) as excinfo:
        tx.seqentrans(frame, avg_occ=True)
    message = str(excinfo.value)
    assert "seqefsub" not in message
    assert "the total the counts were taken on" in message
