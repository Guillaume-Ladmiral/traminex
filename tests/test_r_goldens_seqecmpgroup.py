"""R-golden tests for seqecmpgroup -- the most discriminating subsequences between
groups (weighted group x presence chi-square). Pins the per-group Freq/Resid,
statistic, p.value (+ Bonferroni) and the by-statistic ordering against live
TraMineR 2.2-14. Groups are derived from seqeid (golden-pinned) so the test
reconstructs them identically; the 2-group case exercises Yates' correction.

Generator: scripts/r_goldens/make_goldens_seqecmpgroup.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden, load_golden_lines

from traminex import seqecmpgroup, seqecreate, seqefsub


def _fsub():
    df = pd.read_csv(R_GOLDENS_DIR / "actcal_tse.csv")
    levels = load_golden_lines("eseq_actcal_levels")
    df["event"] = pd.Categorical(df["event"].astype(str), categories=levels)
    ea = seqecreate(df, time_col="time")
    return seqefsub(ea, pmin_support=0.05, max_k=2)


def _group(fsub, k):
    # Mirrors the generator's seq_along(ea) %% k (1-based eseq position).
    return np.arange(1, len(fsub.eseq) + 1) % k


def _assert_cmp(cmp, tag):
    exp_subseq = load_golden_lines("seqecmpgroup_" + tag + "_subseq")
    exp_cols = load_golden_lines("seqecmpgroup_" + tag + "_cols")
    golden = load_golden("seqecmpgroup_" + tag)
    assert cmp.subseq == exp_subseq
    assert list(cmp.data.columns) == exp_cols
    actual = cmp.data[exp_cols].to_numpy(dtype=float)
    np.testing.assert_allclose(actual, golden, rtol=0, atol=R_GOLDEN_ATOL, equal_nan=True)


def test_seqecmpgroup_two_groups_yates():
    fsub = _fsub()
    group = _group(fsub, 2)
    _assert_cmp(seqecmpgroup(fsub, group), "g2")


def test_seqecmpgroup_three_groups():
    fsub = _fsub()
    group = _group(fsub, 3)
    _assert_cmp(seqecmpgroup(fsub, group), "g3")


def test_seqecmpgroup_bonferroni():
    fsub = _fsub()
    group = _group(fsub, 2)
    _assert_cmp(seqecmpgroup(fsub, group, method="bonferroni"), "g2_bonf")


def test_seqecmpgroup_labels_and_bonferroni_record():
    fsub = _fsub()
    group = _group(fsub, 2)
    cmp = seqecmpgroup(fsub, group, method="bonferroni")
    assert cmp.labels == ["0", "1"]
    assert cmp.bonferroni["used"] is True
    assert cmp.bonferroni["ntest"] == len(fsub.subseq)
