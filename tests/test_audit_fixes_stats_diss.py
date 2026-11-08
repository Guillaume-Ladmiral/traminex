"""Regressions for the 2026-09-02 audit's statistics/dissimilarity findings.

R anchors were taken from live TraMineR 2.2-14 with ``format(x, digits = 15)``.
"""

from __future__ import annotations

import inspect

import numpy as np
import pandas as pd
import pytest

import traminex
from traminex.sequences import STSList


def test_seqmodst_weighted_false_ignores_attached_weights():
    # Live R: seqmodst(weighted=FALSE) on weights c(1, 10) -> A-A with freq 0.5.
    data = np.array([[0, 0], [1, 1]])
    weighted = STSList(data.copy(), alphabet=["A", "B"], weights=np.array([1.0, 10.0]))
    modal = traminex.seqmodst(weighted, weighted=False)
    assert modal.to_strings() == ["A-A"]
    np.testing.assert_allclose(modal.frequencies_, [0.5, 0.5])
    modal_weighted = traminex.seqmodst(weighted, weighted=True)
    assert modal_weighted.to_strings() == ["B-B"]
    np.testing.assert_allclose(modal_weighted.frequencies_, [10.0 / 11.0, 10.0 / 11.0])


def test_seqcta_rejects_subseq_longer_than_sequences():
    # R fails on the same input with an unrelated names-length error; a clear
    # validation error replaces the pandas "Length mismatch" crash.
    x = STSList(np.array([[0, 1, 0, 1], [1, 1, 0, 0]]), alphabet=["A", "B"])
    with pytest.raises(ValueError, match="subseq"):
        traminex.seqcta(x, initial_state="A")  # default subseq=5 > seqlen=4
    result = traminex.seqcta(x, subseq=2, initial_state="A")
    assert isinstance(result, pd.DataFrame)


def test_fit_statistics_warns_on_degenerate_degrees_of_freedom():
    # R (dissrf.R:262): msg.warn("Negative  or zero df! ...") next to pvalue <- NA.
    from traminex.repseq_rf import _fit_statistics

    matrix = np.array([[0.0, 1.0], [1.0, 0.0]])
    with pytest.warns(UserWarning, match="df"):
        r_squared, f_statistic, p_value = _fit_statistics(
            matrix, 0, 0.5, exponent=1.0, k=2, weight_sum=2.0
        )
    assert np.isnan(p_value)


@pytest.mark.parametrize(
    ("k", "weight_sum", "text"),
    [(2, 2.0, "k = 2  wsum =  2"), (3, 2.5, "k = 3  wsum =  2.5")],
)
def test_the_degenerate_df_warning_uses_the_dissrf_text(k, weight_sum, text):
    """dissrf.R:262 calls msg.warn("Negative  or zero df! k =", format(k), " wsum = ",
    format(wsum)); msg.warn pastes its parts with spaces after " [!!]"
    (TraMineR-logging_helpers.R:20-25). traminex wrote "Negative or zero df! k = 2, wsum = 2.0"
    (TX-REP-20261023-008)."""
    from traminex.repseq_rf import _fit_statistics

    matrix = np.array([[0.0, 1.0], [1.0, 0.0]])
    with pytest.warns(UserWarning) as caught:
        _fit_statistics(matrix, 0, 0.5, exponent=1.0, k=k, weight_sum=weight_sum)

    assert [str(w.message) for w in caught] == [f" [!!] Negative  or zero df! {text}"]


def test_disstree_get_rules_preserves_r_level_order():
    # R (disstree.R) labels both kids when their parent is expanded and
    # factor(factor(...)) keeps that level order (verified in live R 4.6.1),
    # so a left leaf whose label sorts lexicographically AFTER the right one
    # must still come first.
    from types import SimpleNamespace

    from traminex.disstree_helpers import disstree_get_rules

    split = SimpleNamespace(
        varindex=0, breaks=None, labels=["zzz", "aaa"], index=[1, 2], na_group=None
    )
    left = SimpleNamespace(id=2, kids=None, split=None)
    right = SimpleNamespace(id=3, kids=None, split=None)
    root = SimpleNamespace(id=1, kids=[left, right], split=split)
    tree = SimpleNamespace(root=root, data=pd.DataFrame({"g": []}))
    assert disstree_get_rules(tree) == ["g %in% c('zzz')", "g %in% c('aaa')"]


def test_stars_and_quality_lines_match_r_formatting():
    from types import SimpleNamespace

    from traminex.disstree_viz import _quality_lines, _stars

    # R: as.character(cut(val, c(0, 0.001, 0.01, 0.05, 1), labels=...)).
    assert _stars(0.0005) == "***"
    assert _stars(0.05) == "*"
    assert _stars(0.5) == ""
    assert _stars(0.0) == "NA"
    assert _stars(2.0) == "NA"
    stat = pd.DataFrame(
        {"t0": [12.34, 0.567, 1.23], "p_value": [0.0005, 0.5, 0.04]},
        index=["Pseudo F", "Pseudo R2", "Levene"],
    )
    tree = SimpleNamespace(adjustment=SimpleNamespace(stat=stat))
    lines = _quality_lines(tree)
    # R's paste(n, ": ", value, " ", star, sep="") keeps the separator space
    # even when the star is empty. The values use R's
    # format(x, digits = getOption("digits") - 2), i.e. 5 significant digits
    # (tests/r_goldens/r_format_digits.csv: 12.34 -> "12.34"); this test
    # pinned 3 digits ("12.3") until the 2026-09-22 audit S6-03.
    assert lines == ["Pseudo F: 12.34 ***", "Pseudo R2: 0.567 ", "Levene: 1.23 *"]


def test_tree_dot_writers_default_to_r_use_title():
    from traminex.disstree_viz import disstree2dot, disstree2dotp, seqtree2dot

    for fn in (disstree2dot, disstree2dotp, seqtree2dot):
        assert inspect.signature(fn).parameters["use_title"].default is True, fn.__name__


def test_seqmtplot_whiskers_use_weighted_nbseq():
    # R plot.stslist.meant: df = nbseq - 1 with nbseq = sum of weights.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import t as student_t

    x = STSList(
        np.array([[0, 0], [1, 1], [0, 1]]),
        alphabet=["A", "B"],
        weights=np.array([2.0, 3.0, 5.0]),
    )
    fig, ax = traminex.seqmtplot(x, errbar=True)
    df_frame = traminex.seqmeant(x, serr=True)
    caps = [line for line in ax.lines if line.get_marker() == "_"]
    upper = max(np.max(line.get_ydata()) for line in caps)
    expected_q = float(student_t.ppf(0.975, 10.0 - 1.0))
    expected_upper = float(
        np.max(df_frame["Mean"].to_numpy() + expected_q * df_frame["SE"].to_numpy())
    )
    assert upper == pytest.approx(expected_upper, rel=1e-12)
    plt.close(fig)
