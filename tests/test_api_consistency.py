import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from traminex import (
    STS_to_TSE,
    TSE_to_STS,
    seqdef,
    seqstatd,
    seqtree_fit,
    seqtree_predict,
)
from traminex.events import events_to_state_sequences, seqecreate
from traminex.plotting import sequence_index_plot, state_distribution_plot


def test_sts_to_tse_roundtrip_fullgrid():
    arr = np.array([["A", "A", "B", "B"], ["A", "B", "B", "B"]], dtype=object)
    x = seqdef(arr, start=1)
    tse = STS_to_TSE(x, state_changes_only=False, time_col="time")
    x2 = TSE_to_STS(tse, id_col="id", time_col="time", state_col="state")
    assert x2.data.shape == x.data.shape
    assert np.array_equal(x2.data, x.data)
    assert x2.alphabet == x.alphabet


def test_seqtree_predict_shape_and_values():
    # Small dataset with two labels
    arr = np.array([["A", "A", "B", "B"], ["B", "B", "B", "A"], ["A", "B", "A", "B"]], dtype=object)
    x = seqdef(arr)
    y = np.array([0, 1, 0])
    res = seqtree_fit(x, y)
    pred = seqtree_predict(res, x)
    assert pred.shape == y.shape
    # Predicting on training set should be equal to model's predictions
    from traminex.seqtree_explore import _features_from_sequences

    Xf, _ = _features_from_sequences(x)
    pred_ref = res.model.predict(Xf)
    assert np.array_equal(pred, pred_ref)


def test_events_to_state_sequences_align_missing():
    df = pd.DataFrame({"id": ["p1", "p1"], "timestamp": [2, 4], "event": ["X", "Y"]})
    events = seqecreate(df, id_col="id", time_col="timestamp", event_col="event")
    sts = events_to_state_sequences(events, fill="_", align_missing=True)
    # '_' must be treated as missing
    assert sts.nr == "_"
    assert sts.nr in sts.alphabet
    assert sts.nr_code == sts.alphabet.index("_")


def test_plotting_ax_passthrough():
    arr = np.array([["A", "B", "A"], ["B", "B", "A"]], dtype=object)
    x = seqdef(arr)
    fig1, ax1 = plt.subplots()
    fig_r, ax_r = sequence_index_plot(x, ax=ax1)
    assert ax_r is ax1
    assert fig_r is fig1
    fig2, ax2 = plt.subplots()
    fig_r2, ax_r2 = state_distribution_plot(x, ax=ax2)
    assert ax_r2 is ax2
    assert fig_r2 is fig2
    plt.close(fig1)
    plt.close(fig2)


def test_seqstatd_with_weights():
    arr = np.array([["A", "B"], ["B", "B"]], dtype=object)
    x = seqdef(arr)
    # weights: first sequence has weight 2, second has weight 1
    w = np.array([2.0, 1.0])
    st = seqstatd(x, include_missing=False, weights=w)
    freq = st["Frequencies"]
    # At time 0: states are [A,B] for rows -> weighted counts A=2, B=1 -> proportions [2/3, 1/3]
    a_idx = x.alphabet.index("A")
    b_idx = x.alphabet.index("B")
    assert np.isclose(freq[0, a_idx], 2.0 / 3.0)
    assert np.isclose(freq[0, b_idx], 1.0 / 3.0)
