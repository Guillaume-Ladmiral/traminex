"""R-golden tests for seqeapplysub -- the (n_seq x n_subseq) subsequence-count
matrix. This is where all SIX count methods (COBJ/CDIST_O/CWIN/CMINWIN/CDIST +
age) are pinned end to end against live TraMineR 2.2-14.

Generator: scripts/r_goldens/make_goldens_seqeapplysub.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from r_goldens_utils import R_GOLDEN_ATOL, load_golden, load_golden_lines

from traminex import seqecreate
from traminex.events import seqeapplysub, seqeconstraint, seqefsub


def _fsub():
    # Mirrors scripts/r_goldens/make_goldens_seqeapplysub.R.
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3", "s3"],
            "timestamp": [0, 1, 1, 3, 0, 2, 4, 0, 1, 2],
            "event": ["A", "B", "C", "A", "A", "B", "C", "A", "B", "C"],
        }
    )
    return seqefsub(seqecreate(tse), min_support=1, max_k=3)


def _assert_matrix(actual: pd.DataFrame, name: str):
    golden = load_golden("seqeapplysub_" + name)
    np.testing.assert_allclose(actual.to_numpy(), golden, rtol=0, atol=R_GOLDEN_ATOL)


def test_seqeapplysub_presence():
    _assert_matrix(seqeapplysub(_fsub(), method="presence"), "presence")


def test_seqeapplysub_count():
    _assert_matrix(seqeapplysub(_fsub(), method="count"), "count")


def test_seqeapplysub_age():
    _assert_matrix(seqeapplysub(_fsub(), method="age"), "age")


def test_seqeapplysub_cdist():
    _assert_matrix(seqeapplysub(_fsub(), method="CDIST"), "cdist")


def test_seqeapplysub_cwin_window3():
    m = seqeapplysub(_fsub(), constraint=seqeconstraint(window_size=3, count_method="CWIN"))
    _assert_matrix(m, "cwin_w3")


def test_seqeapplysub_cminwin_window3():
    m = seqeapplysub(_fsub(), constraint=seqeconstraint(window_size=3, count_method="CMINWIN"))
    _assert_matrix(m, "cminwin_w3")


def test_seqeapplysub_rules_presence():
    m = seqeapplysub(_fsub(), method="presence", rules=True)
    _assert_matrix(m, "rules_presence")


def test_seqeapplysub_age_fractional_truncates_like_r():
    # The age method accumulates the first-occurrence age into a C++ `int c`, so
    # fractional ages are truncated toward zero (1.9999 -> 1, 0.7 -> 0). The
    # integer-time fixture above cannot see this. Regression for that bug.
    tsef = pd.DataFrame(
        {
            "id": ["p1", "p1", "p2", "p2", "p3", "p3"],
            "timestamp": [0.7, 1.5, 2.9, 5.0, 1.9999, 2.0],
            "event": ["A", "B", "A", "B", "A", "B"],
        }
    )
    fsubf = seqefsub(seqecreate(tsef), min_support=1, max_k=2)
    assert fsubf.subseq == load_golden_lines("seqeapplysub_age_frac_cols")
    _assert_matrix(seqeapplysub(fsubf, method="age"), "age_frac")


def test_seqeapplysub_cwin_degenerate_window_returns_zero_not_crash():
    # RATIFY (documented divergence): a CWIN window wider than the allowed age span
    # makes the C++ engine call new double[negative] -> R crashes the whole process
    # (std::bad_array_new_length). traminex clamps the window length to 0 and returns
    # 0 (zero valid windows fit) -- the safe, sensible answer. Audit finding.
    m = seqeapplysub(
        _fsub(),
        constraint=seqeconstraint(window_size=10, count_method="CWIN", age_min=0, age_max_end=2),
    )
    assert np.all(np.asarray(m, dtype=float) == 0.0)


def test_seqeapplysub_labels_align():
    fsub = _fsub()
    m = seqeapplysub(fsub, method="presence")
    # columns = subsequence labels, rows = the event-sequence label strings.
    assert list(m.columns) == fsub.subseq
    assert list(m.index) == fsub.eseq.as_character()
