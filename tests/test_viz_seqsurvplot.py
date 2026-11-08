"""RP7a: seqsurvplot (plot.stslist.surv) structure test.

One Kaplan-Meier survival step curve per state. The step's (time, surv) vertices
== the seqsurv series (golden seqsurv_A: stratum, time, n.risk, n.event, surv,
std.err); each curve starts at (0, 1).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from corestats_fixtures import fixture_a
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.stats_p2 import seqsurv
from traminex.visualization.plots.survival import seqsurvplot


def test_seqsurvplot_step_curves_match_seqsurv_golden():
    x = fixture_a()
    fig, ax = seqsurvplot(x)
    golden = load_golden("seqsurv_A")  # legacy six-column R summary fixture
    # group golden rows by stratum (1-indexed state code); strata 1,2,3 = A,B,C
    strata = np.unique(golden[:, 0]).astype(int)
    # Every state with spells has a curve, one whose spells are all censored a
    # flat one, and each runs on to its last follow-up time, as R's plot.survfit
    # draws them (tests/test_r_goldens_survival_extent.py, TX-AUD-20261010-042).
    follow_up = seqsurv(x).follow_up
    assert len(ax.lines) == len(follow_up) >= len(strata) == 3
    for line, s in zip(ax.lines, strata):
        rows = golden[golden[:, 0] == s]
        expected = np.vstack([[0.0, 1.0], np.column_stack([rows[:, 1], rows[:, 4]])])
        last = follow_up[line.get_label()]
        if last > expected[-1, 0]:
            expected = np.vstack([expected, [last, expected[-1, 1]]])
        np.testing.assert_allclose(line.get_xydata(), expected, rtol=0, atol=R_GOLDEN_ATOL)
    for line in ax.lines[len(strata) :]:
        assert np.all(line.get_ydata() == 1.0)
    # Confidence bounds now live in the seven-column analytical result; the
    # basic survival plot still draws only the documented step curves.
    assert not ax.collections


def test_seqsurvplot_survival_axis_is_zero_to_one():
    x = fixture_a()
    fig, ax = seqsurvplot(x)
    lo, hi = ax.get_ylim()
    assert lo == 0.0 and 1.0 <= hi <= 1.1
