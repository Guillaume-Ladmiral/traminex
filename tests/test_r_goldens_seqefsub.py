"""R-golden tests for the faithful prefix-tree frequent-subsequence miner.

Pins the subsequence label set/order + Support/Count of ``seqefsub`` against live
TraMineR 2.2-14 on actcal.tse (several pmin.support/max.k) and a small synthetic
TSE under default / max.gap / CDIST_O constraints. Also pins the structure of the
faithful ``seqeconstraint``.

Generator: scripts/r_goldens/make_goldens_seqefsub.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_GOLDEN_ATOL, R_GOLDENS_DIR, load_golden, load_golden_lines

from traminex import seqecreate, seqefsub
from traminex.events import seqeconstraint


def _actcal_eseq():
    df = pd.read_csv(R_GOLDENS_DIR / "actcal_tse.csv")
    levels = load_golden_lines("eseq_actcal_levels")
    df["event"] = pd.Categorical(df["event"].astype(str), categories=levels)
    return seqecreate(df, time_col="time")


def _syn_eseq():
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"],
            "timestamp": [0, 1, 1, 3, 0, 2, 4, 0, 1, 2],
            "event": ["A", "B", "C", "A", "A", "B", "C", "A", "B", "C"],
        }
    )
    return seqecreate(tse)


def _assert_fsub_matches(fsub, tag):
    labels = load_golden_lines("seqefsub_" + tag + "_subseq")
    golden = load_golden("seqefsub_" + tag)  # columns Support, Count
    assert fsub.subseq == labels
    actual = np.column_stack([fsub.data["Support"].to_numpy(), fsub.data["Count"].to_numpy()])
    np.testing.assert_allclose(actual, golden, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqefsub_actcal_p05_k2():
    _assert_fsub_matches(seqefsub(_actcal_eseq(), pmin_support=0.05, max_k=2), "actcal_p05_k2")


def test_seqefsub_actcal_p02_k3():
    _assert_fsub_matches(seqefsub(_actcal_eseq(), pmin_support=0.02, max_k=3), "actcal_p02_k3")


def test_seqefsub_syn_default():
    _assert_fsub_matches(seqefsub(_syn_eseq(), min_support=1, max_k=3), "syn_k3")


def test_seqefsub_syn_maxgap2():
    fsub = seqefsub(_syn_eseq(), min_support=1, max_k=3, constraint=seqeconstraint(max_gap=2))
    _assert_fsub_matches(fsub, "syn_maxgap2")


def test_seqefsub_syn_cdisto():
    fsub = seqefsub(
        _syn_eseq(), min_support=1, max_k=3, constraint=seqeconstraint(count_method="CDIST_O")
    )
    _assert_fsub_matches(fsub, "syn_cdisto")


def test_seqefsub_str_subseq_syn():
    fsub = seqefsub(_syn_eseq(), str_subseq=["(A)", "(A)-(B)", "(B,C)", "(C)-(A)"])
    _assert_fsub_matches(fsub, "syn_strsub")


def test_seqefsub_str_subseq_actcal():
    fsub = seqefsub(
        _actcal_eseq(),
        str_subseq=["(FullTime)", "(NoActivity)-(Start)", "(Start,LowPartTime)"],
    )
    _assert_fsub_matches(fsub, "actcal_strsub")


def test_seqefsub_str_subseq_cdisto():
    # CDIST_O count.method: Count = occurrence count (distinct from the COBJ Support).
    fsub = seqefsub(
        _actcal_eseq(),
        str_subseq=["(NoActivity)", "(NoActivity)-(Start)"],
        constraint=seqeconstraint(count_method="CDIST_O"),
    )
    _assert_fsub_matches(fsub, "actcal_strsub_cdisto")


def test_seqefsub_str_subseq_unknown_event_errors():
    with pytest.raises(ValueError):
        seqefsub(_syn_eseq(), str_subseq=["(Z)"])


def test_seqeconstraint_structure_and_validation():
    c = seqeconstraint()
    assert (c.max_gap, c.window_size, c.age_min, c.age_max, c.age_max_end, c.count_method) == (
        -1,
        -1,
        -1,
        -1,
        -1,
        1,
    )
    # Name -> integer code (COBJ=1 .. CDIST=5).
    assert seqeconstraint(count_method="CDIST_O").count_method == 2
    assert seqeconstraint(count_method="CWIN").count_method == 3
    assert seqeconstraint(count_method="CMINWIN").count_method == 4
    assert seqeconstraint(count_method="CDIST").count_method == 5
    # age.max.end fills age.max when the latter is unset.
    assert seqeconstraint(age_max_end=5).age_max == 5
    # Validation (mirrors seqeconstraint.R).
    with pytest.raises(ValueError):
        seqeconstraint(max_gap=5, window_size=2)
    with pytest.raises(ValueError):
        seqeconstraint(age_min=5, age_max=2)
    with pytest.raises(ValueError):
        seqeconstraint(count_method="NOPE")
    with pytest.raises(ValueError):
        seqeconstraint(count_method=9)


def test_seqefsub_requires_a_support_threshold():
    with pytest.raises(ValueError):
        seqefsub(_syn_eseq())  # neither min_support nor pmin_support


def test_seqefsub_cwin_requires_window_size():
    with pytest.raises(ValueError):
        seqefsub(_syn_eseq(), min_support=1, constraint=seqeconstraint(count_method="CWIN"))
