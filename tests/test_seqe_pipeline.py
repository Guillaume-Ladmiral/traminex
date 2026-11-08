import pandas as pd

import traminex as tx


def test_seqecreate_and_seqefsub_default():
    tse = pd.DataFrame(
        {
            "id": ["i1", "i1", "i1", "i2", "i2"],
            "timestamp": [1, 2, 4, 1, 2],
            "event": ["A", "B", "A", "A", "B"],
        }
    )
    e = tx.seqecreate(tse)
    fsub = tx.seqefsub(e, min_support=1, max_k=2)
    # The faithful miner returns a SubSeqElist with generic "(A)-(B)" labels.
    assert tx.is_subseqelist(fsub)
    assert "(A)" in fsub.subseq
    assert "(A)-(B)" in fsub.subseq
    assert {"Support", "Count"}.issubset(fsub.data.columns)
    assert len(fsub) == len(fsub.subseq) == fsub.data.shape[0]


def test_seqefsub_noncontig_maxgap_constraint():
    tse = pd.DataFrame(
        {
            "id": ["i1", "i1", "i1", "i2", "i2"],
            "timestamp": [1, 2, 4, 1, 3],
            "event": ["A", "B", "C", "A", "C"],
        }
    )
    e = tx.seqecreate(tse)
    ct = tx.events.seqeconstraint(max_gap=2)
    fsub = tx.seqefsub(e, min_support=1, max_k=3, constraint=ct)
    assert tx.is_subseqelist(fsub)
    # Support is a proportion in [0, 1].
    assert (fsub.data["Support"] >= 0).all() and (fsub.data["Support"] <= 1).all()


def test_seqefsub_pmin_support_proportional():
    tse = pd.DataFrame(
        {
            "id": ["i1", "i1", "i2", "i2", "i3"],
            "timestamp": [1, 2, 1, 2, 1],
            "event": ["A", "B", "A", "B", "A"],
        }
    )
    e = tx.seqecreate(tse)
    # (A) occurs in all 3 sequences -> Support 1.0, Count 3.
    fsub = tx.seqefsub(e, pmin_support=0.5, max_k=1)
    a = fsub.data.iloc[fsub.subseq.index("(A)")]
    assert a["Count"] == 3.0
    assert a["Support"] == 1.0
