import numpy as np

import traminex as tx
import traminex.stats as tx_stats
import traminex.stats_complexity as tx_stats_complexity


def test_stats_facade_exports_canonical_split_module_functions():
    # The split stats facade should expose the same canonical functions as top-level API.
    assert tx_stats.seqici is tx.seqici
    assert tx_stats.seqmeant is tx.seqmeant
    assert tx_stats.seqpm is tx.seqpm
    assert tx_stats_complexity.seqpm is tx.seqpm


def test_stats_facade_split_modules_work_together():
    arr = np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
            ["C", "C", "B", "A"],
        ],
        dtype=object,
    )
    x = tx.seqdef(arr, nr="*")

    # basic module path
    assert np.array_equal(tx_stats.seqlength(x), np.array([4, 4, 4]))
    # utils module path: seqnum recodes the alphabet to numeric strings
    assert list(tx_stats.seqnum(x).alphabet[:3]) == ["0", "1", "2"]
    # complexity path
    t = tx_stats.seqST(x, norm=True)
    assert t.shape == (x.nseq,)
