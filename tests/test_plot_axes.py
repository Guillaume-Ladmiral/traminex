import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import traminex as tx


def test_sequence_index_plot_time_axis_labels():
    df = pd.DataFrame(
        {
            "t1": ["A", "B"],
            "t2": ["B", "B"],
            "t3": ["A", "B"],
            "t4": ["A", "A"],
            "t5": ["B", "A"],
        }
    )
    x = tx.seqdef(df, xtstep=2, start=2000)
    fig, ax = tx.sequence_index_plot(x)
    labs = [t.get_text() for t in ax.get_xticklabels()]
    # Expect 2000, 2002, 2004 labels (approx)
    assert labs[0] in ("2000", "2001")  # depending on matplotlib alignment
    assert any(l in labs for l in ["2002", "2003"]) and any(l in labs for l in ["2004", "2005"])
    plt.close(fig)
