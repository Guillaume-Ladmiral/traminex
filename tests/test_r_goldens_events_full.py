"""R-golden tests for the event framework -- ``seqentrans`` augmentation (T1).

``TraMineRextras::seqentrans`` does no mining; it augments an already-mined
``seqefsub`` table with ``ntrans`` / ``nevent`` / ``Avg.occ`` computed from the
subsequence *label strings*. We anchor traminex's augmentation to TraMineR's own
(structurally varied -- multi-event, multi-transition) labels, so the test is
faithful to the function under fix regardless of traminex's (separate) miner.

Generator: scripts/r_goldens/make_goldens_seqentrans.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden

from traminex import seqecreate, seqefsub, seqentrans


def test_seqentrans_ntrans_nevent_avgocc_match_r():
    """The augmentation reproduces R's ntrans/nevent/Avg.occ on R's own labels."""
    labels = (R_GOLDENS_DIR / "seqentrans_subseq.txt").read_text(encoding="utf-8").splitlines()
    labels = [s.strip() for s in labels if s.strip()]
    wtot = float((R_GOLDENS_DIR / "seqentrans_wtot.txt").read_text(encoding="utf-8").strip())
    golden = load_golden("seqentrans")  # columns: Count, ntrans, nevent, Avg.occ
    assert len(labels) == golden.shape[0]

    # Feed R's subsequence label strings + counts through traminex's augmentation
    # (R's seqentrans likewise operates purely on fsubseq$subseq).
    fsub = pd.DataFrame({"subseq": labels, "count": golden[:, 0]})
    fsub.attrs["seqe_weight_total"] = wtot

    ft = seqentrans(fsub, avg_occ=True)
    actual = np.column_stack(
        [
            ft["count"].to_numpy(),
            ft["ntrans"].to_numpy(),
            ft["nevent"].to_numpy(),
            ft["avg_occ"].to_numpy(),
        ]
    )
    np.testing.assert_allclose(actual, golden, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqentrans_native_pipeline_augments_table():
    """seqecreate -> seqefsub -> seqentrans augments the faithful miner's result.

    The faithful miner yields generic subsequence label strings ("(A)-(B)" etc.);
    seqentrans reads those labels (ntrans = number of "-"-transitions, nevent =
    ntrans - 1 + comma tokens) and Avg.occ = Count / total sequence weight.
    """
    tse = pd.DataFrame(
        {
            "id": [1, 1, 1, 2, 2, 3, 3, 3],
            "timestamp": [1, 2, 3, 1, 2, 1, 2, 3],
            "event": ["A", "B", "C", "A", "B", "A", "B", "C"],
        }
    )
    eseq = seqecreate(tse)
    fsub = seqefsub(eseq, min_support=1, max_k=3)
    # seqentrans returns an augmented copy and leaves its argument unchanged, as R's
    # copy-on-modify argument does (live R: the argument keeps Support and Count only;
    # TX-AUD-20261019-016).
    ft = seqentrans(fsub, avg_occ=True)

    assert ft is not fsub
    assert list(fsub.data.columns) == ["Support", "Count"]
    assert {"ntrans", "nevent", "Avg.occ"}.issubset(ft.data.columns)
    d = ft.data
    for i in range(len(d)):
        # Single-event transitions only here, so nevent == ntrans.
        assert d["nevent"].iloc[i] == d["ntrans"].iloc[i]
        assert d["Avg.occ"].iloc[i] == d["Count"].iloc[i] / 3.0  # 3 sequences, unweighted

    # "(A)-(B)" is present in all 3 sequences -> Count 3, Avg.occ 1.0.
    idx = ft.subseq.index("(A)-(B)")
    ab = ft.data.iloc[idx]
    assert ab["Count"] == 3.0
    assert ab["ntrans"] == 2 and ab["nevent"] == 2
    assert ab["Avg.occ"] == 1.0
