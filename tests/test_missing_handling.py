import numpy as np
import pandas as pd

from traminex.stats import seqstatd


def test_nr_void_handling_in_stats():
    # include '*' as NR token in data
    wide = pd.DataFrame(
        {
            1: ["A", "*", "A"],
            2: ["B", "B", "*"],
        }
    )
    from traminex import seqdef

    x = seqdef(wide, nr="*", void="%")
    st_inc = seqstatd(x, include_missing=True)
    st_exc = seqstatd(x, include_missing=False)
    # Frequencies differ when excluding missing
    assert not np.allclose(st_inc["Frequencies"], st_exc["Frequencies"])


def test_nr_void_handling_in_seqistatd():
    wide = pd.DataFrame(
        {
            1: ["A", "*", "A"],
            2: ["B", "B", "*"],
        }
    )
    from traminex import seqdef
    from traminex.stats import seqistatd

    x = seqdef(wide, nr="*", void="%")
    d_inc = seqistatd(x, include_missing=True)
    d_exc = seqistatd(x, include_missing=False)

    nr_idx = x.nr_code
    assert nr_idx is not None
    assert d_inc[:, nr_idx].sum() > 0
    assert d_exc[:, nr_idx].sum() == 0
