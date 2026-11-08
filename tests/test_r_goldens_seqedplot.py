"""R-golden tests for the ``seqedplot`` producing series (TraMineRextras
``seqedplot``). The plot itself is structure-tested elsewhere; here the two
numeric series it draws are pinned to live R exactly as ``seqedplot`` computes
them (TraMineR 2.2-14 / TraMineRextras 0.6.9 / survival 3.8.6).

- survival: per non-ignored event, the unweighted Kaplan-Meier curve of the time
  to the event's first occurrence, censored at ``seqelength`` (a weighted fixture
  proves the survival branch ignores weights);
- hazard: added in the hazard commit.

Generator: scripts/r_goldens/make_goldens_seqedplot.R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from r_goldens_utils import R_DIFFERENCE_ATOL, R_GOLDEN_ATOL, assert_matches_golden, load_golden

from traminex.events import seqecreate
from traminex.events._duration import (
    _seqedplot_ages_range,
    _seqedplot_cutpoints,
    _seqedplot_hazard,
    _seqedplot_survival,
    _subset_seqelist,
)

_OBS_LEN = [5, 5, 6, 4, 5]
_WTS = [1.5, 2, 0.5, 3, 1]


def _surv_tse():
    return pd.DataFrame(
        {
            "id": ["s1", "s1", "s2", "s3", "s3", "s4", "s4", "s5"],
            "timestamp": [1, 3, 2, 1, 4, 1, 3, 2],
            "event": ["A", "B", "A", "B", "A", "A", "A", "B"],
        }
    )


def _surv_eseq(weighted=False):
    e = seqecreate(_surv_tse())
    e.lengths = _OBS_LEN  # seqelength<-
    if weighted:
        e.weights = _WTS  # seqeweight<-
    return e


def test_survival_series_unweighted():
    series = dict(_seqedplot_survival(_surv_eseq()))
    assert_matches_golden(series["A"], "seqedplot_surv_A")
    assert_matches_golden(series["B"], "seqedplot_surv_B")


def test_survival_series_ignores_weights():
    # The survival branch calls survfit unweighted: weighted == unweighted golden.
    series = dict(_seqedplot_survival(_surv_eseq(weighted=True)))
    assert_matches_golden(series["A"], "seqedplot_surv_A_wtd")
    assert_matches_golden(series["B"], "seqedplot_surv_B_wtd")
    # and identical to the unweighted golden
    assert_matches_golden(series["A"], "seqedplot_surv_A")
    assert_matches_golden(series["B"], "seqedplot_surv_B")


def test_survival_series_order_and_shape():
    series = _seqedplot_survival(_surv_eseq())
    labels = [lab for lab, _ in series]
    assert labels == ["A", "B"]  # kept-event (dictionary) order
    for _, km in series:
        assert km.shape[1] == 5  # time, n.risk, n.event, surv, std.err


def test_survival_series_ignore_list():
    series = _seqedplot_survival(_surv_eseq(), ignore=["B"])
    labels = [lab for lab, _ in series]
    assert labels == ["A"]
    assert_matches_golden(dict(series)["A"], "seqedplot_surv_A")


def test_survival_requires_observation_time():
    # R: stop(...) when all seqelength == -1 (the default for TSE-created eseq).
    e = seqecreate(_surv_tse())  # lengths unset (-1)
    with pytest.raises(ValueError, match="observation time"):
        _seqedplot_survival(e)


# --------------------------------------------------------------------------- #
# Hazard series                                                               #
# --------------------------------------------------------------------------- #
_HWTS = [1.0, 2.0, 0.5, 1.5, 2.5, 1.0]
_HGROUP_IDX = [[0, 1, 2], [3, 4, 5]]  # group X = s1..s3, group Y = s4..s6


def _haz_tse():
    return pd.DataFrame(
        {
            "id": [
                "s1",
                "s1",
                "s1",
                "s2",
                "s2",
                "s3",
                "s3",
                "s4",
                "s4",
                "s4",
                "s5",
                "s5",
                "s6",
                "s6",
            ],
            "timestamp": [1, 3, 5, 2, 4, 1, 6, 1, 2, 4, 3, 5, 5, 2],
            "event": ["A", "A", "B", "A", "B", "B", "A", "A", "A", "A", "B", "B", "A", "B"],
        }
    )


def _haz_eseq():
    e = seqecreate(_haz_tse())
    e.weights = _HWTS  # seqeweight<-
    return e


def test_hazard_cutpoints_reproduce_r_cut():
    # The delicate part: R's cut(ages, breaks=N) interval boundaries. ages=(1,6).
    _, onx5 = _seqedplot_cutpoints((1.0, 6.0), 5)
    np.testing.assert_allclose(onx5, [2, 3, 4, 5, 6], rtol=0, atol=R_GOLDEN_ATOL)
    cp5, _ = _seqedplot_cutpoints((1.0, 6.0), 5)
    np.testing.assert_allclose(cp5, [0.995, 2, 3, 4, 5, 6], rtol=0, atol=R_GOLDEN_ATOL)
    _, onx3 = _seqedplot_cutpoints((1.0, 6.0), 3)
    np.testing.assert_allclose(onx3, [2.67, 4.33, 6], rtol=0, atol=R_GOLDEN_ATOL)


def test_hazard_global_ages_range():
    assert _seqedplot_ages_range(_haz_eseq()) == (1.0, 6.0)


def test_hazard_ungrouped_b5():
    onx, series = _seqedplot_hazard(_haz_eseq(), breaks=5)
    golden = load_golden("seqedplot_haz_all_b5")
    np.testing.assert_allclose(onx, golden[0], rtol=0, atol=R_GOLDEN_ATOL)
    d = dict(series)
    np.testing.assert_allclose(d["A"], golden[1], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(d["B"], golden[2], rtol=0, atol=R_GOLDEN_ATOL)


def _assert_grouped_hazard(breaks, tag):
    e = _haz_eseq()
    ages = _seqedplot_ages_range(e)  # global cutpoints shared across groups
    golden = load_golden(tag)
    for gi, idx in enumerate(_HGROUP_IDX):
        onx, series = _seqedplot_hazard(_subset_seqelist(e, idx), ages=ages, breaks=breaks)
        np.testing.assert_allclose(onx, golden[0], rtol=0, atol=R_GOLDEN_ATOL)
        d = dict(series)
        np.testing.assert_allclose(d["A"], golden[1 + gi * 2], rtol=0, atol=R_GOLDEN_ATOL)
        np.testing.assert_allclose(d["B"], golden[2 + gi * 2], rtol=0, atol=R_GOLDEN_ATOL)


def test_hazard_grouped_b5():
    _assert_grouped_hazard(5, "seqedplot_haz_b5")


def test_hazard_grouped_b3():
    _assert_grouped_hazard(3, "seqedplot_haz_b3")


def test_hazard_weights_are_load_bearing():
    # Unlike survival, the hazard branch uses seqeweight: dropping weights changes
    # the series (guards the Ph4 "ignore weights -> red" mutation at the source).
    e = _haz_eseq()
    _, weighted = _seqedplot_hazard(e, breaks=5)
    e_unit = seqecreate(_haz_tse())  # default unit weights
    _, unit = _seqedplot_hazard(e_unit, breaks=5)
    assert not np.allclose(dict(weighted)["A"], dict(unit)["A"], rtol=0, atol=R_DIFFERENCE_ATOL)


# --------------------------------------------------------------------------- #
# Corner cases (Ph6): ages clipping; an all-censored (never-occurring) event   #
# --------------------------------------------------------------------------- #
def test_survival_ages_clip():
    # ages=(1.5, 4.5): first-occurrence ages of 1 fall outside the window and are
    # censored (seqedplot.R:131-134).
    series = dict(_seqedplot_survival(_surv_eseq(), ages=(1.5, 4.5)))
    assert_matches_golden(series["A"], "seqedplot_surv_clip_A")
    assert_matches_golden(series["B"], "seqedplot_surv_clip_B")


def test_hazard_ages_clip():
    # ages=(1.5, 5.5): occurrences at ages 1 and 6 are dropped; the cut-points come
    # from the provided window, not the data range.
    onx, series = _seqedplot_hazard(_haz_eseq(), ages=(1.5, 5.5), breaks=5)
    golden = load_golden("seqedplot_haz_clip_b5")  # rows: onx, A, B
    np.testing.assert_allclose(onx, golden[0], rtol=0, atol=R_GOLDEN_ATOL)
    d = dict(series)
    np.testing.assert_allclose(d["A"], golden[1], rtol=0, atol=R_GOLDEN_ATOL)
    np.testing.assert_allclose(d["B"], golden[2], rtol=0, atol=R_GOLDEN_ATOL)


def test_survival_all_censored_event_is_flat():
    # An event that never occurs (a declared-but-absent level): every sequence is
    # censored, so summary(survfit) has no event rows -> an empty (0, 5) curve
    # (survival flat at 1). Verified against live R in the Ph6 corner audit.
    tse = pd.DataFrame(
        {
            "id": ["s1", "s1", "s2", "s3"],
            "timestamp": [1, 3, 2, 1],
            "event": pd.Categorical(["A", "B", "A", "B"], categories=["A", "B", "C"]),
        }
    )
    e = seqecreate(tse)
    e.lengths = [5, 5, 6]
    series = dict(_seqedplot_survival(e))
    assert series["C"].shape == (0, 5)  # never occurs -> no event rows
