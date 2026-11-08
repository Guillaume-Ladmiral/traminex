import os
import warnings

import numpy as np
import pandas as pd
import pytest

from traminex import seqdef, seqdist
from traminex.distances.python_kernels import pairwise as py_pairwise
from traminex.stats import seqstatd

pytestmark = pytest.mark.skipif(
    not os.path.exists("tests/references"), reason="No references directory"
)


def test_parity_refs_if_present():
    base = "tests/references"
    input_path = os.path.join(base, "input_wide.csv")
    # Honesty guard: tests/references currently ships only a README, so without
    # this skip the test would PASS having asserted nothing (a misleading green).
    # Generate the reference CSVs (see tests/references/README.md) to activate it.
    if not os.path.exists(input_path):
        pytest.skip(
            "no reference CSVs in tests/references (see its README); "
            "this optional parity test is inert until they are generated"
        )
    freq_ref = os.path.join(base, "expected_seqstatd_freq.csv")
    dist_ref = os.path.join(base, "expected_dist_om.csv")
    trate_ref = os.path.join(base, "expected_trate_subm.csv")
    dhd_ref = os.path.join(base, "expected_dhd.csv")
    available = [path for path in (freq_ref, dist_ref, trate_ref, dhd_ref) if os.path.exists(path)]
    if not available:
        pytest.skip(
            "input_wide.csv has no complete expected comparison; numerical parity not verified"
        )
    missing = [
        path for path in (freq_ref, dist_ref, trate_ref, dhd_ref) if not os.path.exists(path)
    ]
    if missing:
        warnings.warn("partial reference coverage; missing " + ", ".join(missing), UserWarning)
    if os.path.exists(input_path) and os.path.exists(freq_ref):
        df = pd.read_csv(input_path)
        x = seqdef(df)
        st = seqstatd(x, include_missing=False)
        freq = st["Frequencies"]
        ref = np.loadtxt(freq_ref, delimiter=",")
        assert freq.shape == ref.shape
        assert np.allclose(freq, ref, atol=1e-6)
    if os.path.exists(input_path) and os.path.exists(dist_ref):
        df = pd.read_csv(input_path)
        x = seqdef(df)
        D = py_pairwise(x.data[:2], method="OM", indel=1.0, subm=2.0)
        refD = np.loadtxt(dist_ref, delimiter=",")
        assert D.shape == refD.shape
        assert np.allclose(D, refD, atol=1e-6)
    if os.path.exists(input_path) and os.path.exists(trate_ref):
        df = pd.read_csv(input_path)
        x = seqdef(df)
        from traminex.distances.costs import trate_substitution_matrix

        sm = trate_substitution_matrix(x.data, smoothing=1e-6, mode="symmetric", normalize="none")
        refsm = np.loadtxt(trate_ref, delimiter=",")
        assert sm.shape == refsm.shape
        assert np.allclose(sm, refsm, atol=1e-6)
    if os.path.exists(input_path) and os.path.exists(dhd_ref):
        df = pd.read_csv(input_path)
        x = seqdef(df)
        # R's default DHD uses transition costs, not the low-level frequency
        # extension selected by variant="mean". Keep the complete alphabet
        # when subsetting, just as R seq[1:2,] does.
        subset = seqdef(df.iloc[:2], alphabet=x.alphabet)
        Dd = np.asarray(seqdist(subset, method="DHD"))
        refd = np.loadtxt(dhd_ref, delimiter=",")
        assert Dd.shape == refd.shape
        assert np.allclose(Dd, refd, atol=2e-3)
