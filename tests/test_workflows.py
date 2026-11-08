import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from traminex import (
    STS_to_TSE,
    TSE_to_STS,
    disstree2dot_explore,
    disstree_feature_importance,
    disstree_fit,
    disstree_leaf_summary,
    disstree_predict,
    events_to_state_sequences,
    seqdef,
    seqdiff,
    seqdist,
    seqecreate,
    seqefsub,
    seqerules,
    seqST,
    seqstatd,
    seqtree_fit,
    seqtree_leaf_summary,
    seqtree_predict,
)
from traminex.cluster import cluster_agglomerative, cluster_representatives, cluster_summary
from traminex.plotting import save_plot, sequence_index_plot_grouped, state_distribution_plot


def test_end_to_end_sequence_analysis_workflow():
    seqmat = np.array(
        [
            ["A", "A", "B", "C", "C"],
            ["A", "B", "B", "C", "C"],
            ["A", "A", "B", "B", "C"],
            ["C", "C", "B", "A", "A"],
            ["C", "B", "B", "A", "A"],
            ["C", "C", "A", "A", "A"],
        ],
        dtype=object,
    )
    groups = np.array([0, 0, 0, 1, 1, 1], dtype=int)

    x = seqdef(seqmat, start=2000)

    st = seqstatd(x)
    assert st["Frequencies"].shape == (x.seqlen, len(x.alphabet))
    assert np.allclose(st["Frequencies"].sum(axis=1), 1.0)

    turb = seqST(x, norm=True)
    assert turb.shape == (x.nseq,)
    assert np.all((turb >= 0.0) & (turb <= 1.0))

    D = seqdist(x, method="HAM")
    assert D.shape == (x.nseq, x.nseq)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)

    cl = cluster_agglomerative(D, n_clusters=2)
    assert cl.labels_.shape == (x.nseq,)
    assert np.unique(cl.labels_).size == 2

    summ = cluster_summary(x, cl.labels_, D)
    assert int(summ["size"].sum()) == x.nseq
    assert {"cluster", "size", "medoid_index"}.issubset(summ.columns)

    reps = cluster_representatives(x, cl.labels_, D)
    assert {"cluster", "rep_index", "rep_sequence"}.issubset(reps.columns)
    assert reps["rep_index"].between(0, x.nseq - 1).all()

    tree = seqtree_fit(x, groups, max_depth=2, min_samples_leaf=1, random_state=0)
    pred = seqtree_predict(tree, x)
    assert pred.shape == groups.shape
    leaf = seqtree_leaf_summary(tree)
    assert set(["leaf", "size"]).issubset(leaf.columns)
    assert int(leaf["size"].sum()) == x.nseq

    diff = seqdiff(x, groups)
    assert set(diff.keys()) == {"stat", "discrepancy"}
    assert not diff["stat"].empty
    assert "Total" in diff["discrepancy"].columns


def test_end_to_end_event_pipeline_workflow():
    tse = pd.DataFrame(
        {
            "id": [
                "p1",
                "p1",
                "p1",
                "p2",
                "p2",
                "p2",
                "p3",
                "p3",
                "p3",
            ],
            "timestamp": [2, 3, 4, 2, 3, 4, 2, 3, 4],
            "event": ["A", "B", "C", "A", "B", "D", "A", "B", "C"],
        }
    )

    events = seqecreate(tse)
    fsub = seqefsub(events, min_support=2, max_k=2)
    assert len(fsub) > 0
    assert {"Support", "Count"}.issubset(fsub.data.columns)
    assert "(A)-(B)" in fsub.subseq

    rules = seqerules(fsub)
    assert not rules.empty
    assert {"Rules", "Support", "Conf", "Lift"}.issubset(rules.columns)
    assert "(A) => (B)" in set(rules["Rules"])

    sts = events_to_state_sequences(events, fill="_", align_missing=True, time_grid=np.arange(1, 5))
    assert sts.nseq == 3
    assert sts.seqlen == 4
    assert sts.nr == "_"
    assert sts.nr_code is not None

    D = seqdist(sts, method="HAM", with_missing=True)
    assert D.shape == (sts.nseq, sts.nseq)
    assert np.allclose(D, D.T)


def test_end_to_end_plot_save_workflow(tmp_path):
    seqmat = np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
            ["C", "B", "A", "A"],
            ["C", "C", "A", "A"],
        ],
        dtype=object,
    )
    labels = np.array([0, 0, 1, 1], dtype=int)
    x = seqdef(seqmat, start=2010)

    fig1, ax1, order = sequence_index_plot_grouped(x, labels=labels)
    assert ax1 is not None
    assert order.shape == (x.nseq,)
    out1 = tmp_path / "index_grouped.png"
    save_plot((fig1, ax1), str(out1), dpi=120)
    assert out1.exists()
    assert out1.stat().st_size > 0
    plt.close(fig1)

    fig2, ax2 = state_distribution_plot(x, mode="dH")
    assert ax2 is not None
    out2 = tmp_path / "distribution.svg"
    save_plot(fig2, str(out2))
    assert out2.exists()
    assert out2.stat().st_size > 0
    plt.close(fig2)


def test_end_to_end_long_to_disstree_roundtrip_workflow(tmp_path):
    long_df = pd.DataFrame(
        {
            "id": [
                "p1",
                "p1",
                "p1",
                "p1",
                "p2",
                "p2",
                "p2",
                "p2",
                "p3",
                "p3",
                "p3",
                "p3",
                "p4",
                "p4",
                "p4",
                "p4",
            ],
            "time": [1, 2, 3, 4] * 4,
            "state": [
                "A",
                "A",
                "B",
                "C",
                "A",
                "B",
                "B",
                "C",
                "C",
                "C",
                "B",
                "A",
                "C",
                "B",
                "A",
                "A",
            ],
        }
    )
    y = np.array([0, 0, 1, 1], dtype=int)

    sts = TSE_to_STS(long_df, id_col="id", time_col="time", state_col="state")
    assert sts.nseq == 4
    assert sts.seqlen == 4

    D = seqdist(sts, method="HAM")
    assert D.shape == (4, 4)
    assert np.allclose(D, D.T)
    assert np.allclose(np.diag(D), 0.0)

    tree = disstree_fit(D, y, n_components=2, max_depth=2, min_samples_leaf=1, random_state=0)
    pred = disstree_predict(tree, D)
    assert pred.shape == y.shape

    leaf = disstree_leaf_summary(tree)
    assert {"leaf", "size"}.issubset(leaf.columns)
    assert int(leaf["size"].sum()) == sts.nseq

    fi = disstree_feature_importance(tree)
    assert {"feature", "importance"}.issubset(fi.columns)
    assert (fi["importance"] >= 0).all()

    dot_path = tmp_path / "disstree.dot"
    disstree2dot_explore(tree, str(dot_path))
    assert dot_path.exists()
    assert dot_path.stat().st_size > 0

    back = STS_to_TSE(sts, state_changes_only=False)
    assert {"id", "time", "state"}.issubset(back.columns)
    assert len(back) == sts.nseq * sts.seqlen
