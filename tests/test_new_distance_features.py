import numpy as np
import pandas as pd
import pytest

import traminex as tx
from traminex.distances import seqsubm
from traminex.distances.costs import indel_from_subm


def make_simple_sts():
    df = pd.DataFrame(
        {
            "t1": ["A", "A"],
            "t2": ["B", "A"],
            "t3": ["B", "B"],
            "t4": ["C", "B"],
        }
    )
    return tx.seqdef(df)


def test_wham_matches_manual():
    x = make_simple_sts()
    # Manual weighted Hamming between two sequences
    w = np.array([1.0, 2.0, 1.0, 1.0])
    a = x.data[0]
    b = x.data[1]
    expected = float(np.dot(w, (a != b).astype(float)) / w.sum())
    D = tx.seqdist(x, method="WHAM", weights=w)
    assert np.isclose(D[0, 1], expected)


def test_dur_metrics_basic_properties():
    x = make_simple_sts()
    for metric in ("l1", "cos", "chi2"):
        D = tx.seqdist(x, method="DUR", metric=metric)
        assert D.shape == (x.nseq, x.nseq)
        assert np.allclose(np.diag(D), 0.0)
        assert np.all(D >= 0)


def test_jsd_trans_zero_for_identical():
    # Identical sequences yield 0 divergence
    arr = np.array([[0, 1, 0, 1], [0, 1, 0, 1]], dtype=int)
    D = tx.seqdist(arr, method="JSD_TRANS")
    assert np.isclose(D[0, 1], 0.0)


# R's diagnostics this test meets by design (TX-AUD-20261019-068).
@pytest.mark.filterwarnings(
    r"ignore:\s*\[!!\]\ At\ least\ one\ substitution\ cost\ greater\ than\ twice\ the\ indel\ cost\.:traminex.distances.python_kernels.CostDiagnosticMessage"
)
def test_seqsubm_constant_and_indelslog():
    x = make_simple_sts()
    # CONSTANT costs
    M = seqsubm(x.data, method="CONSTANT", cval=2.0)
    assert M.shape == (len(x.alphabet), len(x.alphabet))
    assert np.all(np.diag(M) == 0.0)
    assert np.all(M[np.triu_indices_from(M, 1)] == 2.0)
    # INDELSLOG matrix + derived indel
    M2 = seqsubm(x.data, method="INDELSLOG")
    indel = indel_from_subm(M2)
    D = tx.seqdist(x, method="OM", indel=indel, subm=M2)
    assert D.shape == (x.nseq, x.nseq)
