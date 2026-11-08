"""R-faithful parity for ``seqsubm`` TRATE and CONSTANT substitution costs.

The oracle (traminex) must match TraMineR 2.2-14 for the cost matrices the family
relies on. The broad live-R golden suite (``test_r_goldens_distances.py``) covers
INDELSLOG / FUTURE / FEATURES but **not** TRATE or CONSTANT; the only R-anchored
TRATE check lived downstream in pyweightedcluster (rpy2, skipped without R). This
module pins TRATE/CONSTANT inside the oracle itself, with no live-R dependency, by
replicating R's exact formula from the vendored source.

R source (``SOURCES/TraMineR/R/seqcost.R`` + ``seqtrate.R``), non-time-varying:

    tr  = seqtrate(seqdata, weighted=..., lag=1)        # row-normalized rates
    cost[i, j] = cval - tr[i, j] - tr[j, i]   (cval = 2), diag(cost) = 0
    indel      = 0.5 * max(cost)

Crucially, ``seqtrate`` sets a *source state with no observed transitions*
(``PA == 0``) to a **zero** rate row (``if (PA == 0) tmat[x, ] <- 0``) — not a
uniform 1/K row. traminex previously injected ``smoothing = 1e-9`` by default,
which (a) perturbed every rate by ~1e-9 and (b) turned an empty source row into a
uniform ``1/K`` row, corrupting the cost of any state that never acts as a
transition source by up to ``1/K`` (~0.33 for K=3). Fixed: TRATE defaults to
``smoothing = 0`` and zeroes empty source rows exactly as R does.
"""

from __future__ import annotations

import numpy as np
import pytest
from corestats_fixtures import fixture_a as _cs_fixture_a
from corestats_fixtures import fixture_b as _cs_fixture_b
from corestats_fixtures import fixture_c as _cs_fixture_c
from r_goldens_utils import R_GOLDEN_ATOL, load_golden

from traminex.distances import seqsubm

# 6x8 ABC fixture shared with the other distance goldens (every state is a
# transition source) plus per-sequence weights (mirrors the pyweightedcluster
# A3 parity fixture so both sides exercise the same numbers).
_SEQS = np.array(
    [
        [0, 0, 0, 1, 1, 2, 2, 2],
        [0, 1, 0, 1, 2, 1, 2, 0],
        [2, 2, 2, 2, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 2, 2],
        [1, 1, 1, 1, 1, 1, 1, 1],
        [2, 0, 2, 0, 1, 2, 1, 0],
    ],
    dtype=int,
)
_WEIGHTS = np.array([1.0, 2.0, 0.5, 1.5, 3.0, 1.0], dtype=float)

# Asymmetric edge fixture: state 2 ("C") appears only in the final column, so it
# is never a transition *source* (R's PA == 0 for that row).
_SEQS_ABSENT_SOURCE = np.array(
    [
        [0, 0, 1, 1, 0, 2],
        [1, 0, 1, 0, 1, 2],
        [0, 1, 0, 1, 0, 2],
    ],
    dtype=int,
)


def _r_trate_cost(seqmat: np.ndarray, weights: np.ndarray | None, cval: float = 2.0) -> np.ndarray:
    """Replicate R seqcost(method='TRATE')$sm exactly (no smoothing, PA==0 -> 0 row)."""
    n, L = seqmat.shape
    K = int(seqmat.max()) + 1
    w = np.ones(n, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    PA = np.zeros(K, dtype=float)
    PAB = np.zeros((K, K), dtype=float)
    for i in range(n):
        for t in range(L - 1):
            x = int(seqmat[i, t])
            y = int(seqmat[i, t + 1])
            PA[x] += w[i]
            PAB[x, y] += w[i]
    tr = np.zeros((K, K), dtype=float)
    for x in range(K):
        if PA[x] > 0:
            tr[x] = PAB[x] / PA[x]
    cost = cval - (tr + tr.T)
    np.fill_diagonal(cost, 0.0)
    return cost


@pytest.mark.parametrize("weighted", [False, True])
def test_seqsubm_trate_matches_r_formula(weighted):
    """seqsubm(method='TRATE') == R seqcost(method='TRATE'); weighted path included."""
    w = _WEIGHTS if weighted else None
    sm = np.asarray(seqsubm(_SEQS, method="TRATE", weights=w), dtype=float)
    ref = _r_trate_cost(_SEQS, w)
    assert sm.shape == ref.shape
    np.testing.assert_allclose(sm, ref, rtol=0, atol=1e-12)
    # R: indel = 0.5 * max(cost)
    assert 0.5 * float(sm.max()) == pytest.approx(0.5 * float(ref.max()), abs=1e-12)


def test_seqsubm_trate_absent_source_state_is_r_faithful():
    """A state that is never a transition source gets R's zero rate row, not 1/K.

    This is the load-bearing asymmetric fixture: with the old default smoothing the
    empty source row became uniform 1/K, biasing every cost involving that state by
    ~1/K. The R formula zeroes the row, so cost[i, absent] = 2 - tr[i, absent] - 0.
    """
    sm = np.asarray(seqsubm(_SEQS_ABSENT_SOURCE, method="TRATE"), dtype=float)
    ref = _r_trate_cost(_SEQS_ABSENT_SOURCE, None)
    np.testing.assert_allclose(sm, ref, rtol=0, atol=1e-12)
    # Sanity: the absent-source state (2) really has no observed outgoing transition.
    assert not np.any(_SEQS_ABSENT_SOURCE[:, :-1] == 2)


def test_seqsubm_trate_positive_smoothing_keeps_empty_source_row_zero():
    """Even with explicit smoothing>0, an empty *source* row stays zero (R's PA==0).

    Probed in ``mode="row"`` where ``cost[i, j] = 1 - T[i, j]``: an absent source
    row therefore reads exactly 1.0 off-diagonal regardless of smoothing. The
    pre-fix uniform-1/K row would have read ``1 - 1/K = 0.667`` instead. This is
    the precise, smoothing-independent guarantee the fix provides; the symmetric
    cost still carries O(smoothing) from the *column* (j -> absent) rates of
    present rows, which is correct (R only smooths nothing, traminex regularizes
    only non-empty rows).
    """
    from traminex.distances.costs import trate_substitution_matrix

    K = int(_SEQS_ABSENT_SOURCE.max()) + 1
    absent = K - 1
    expected_row = np.ones(K, dtype=float)
    expected_row[absent] = 0.0  # diagonal
    for s in (0.0, 1e-3, 0.5):
        row_cost = trate_substitution_matrix(_SEQS_ABSENT_SOURCE, smoothing=s, mode="row")
        np.testing.assert_allclose(row_cost[absent], expected_row, rtol=0, atol=1e-12)


def test_seqsubm_constant_matches_r():
    """seqsubm(method='CONSTANT', cval=2) == R seqcost(method='CONSTANT'): cval off-diag, 0 diag."""
    K = int(_SEQS.max()) + 1
    for cval in (2.0, 3.5):
        sm = np.asarray(seqsubm(_SEQS, method="CONSTANT", cval=cval), dtype=float)
        ref = np.full((K, K), cval, dtype=float)
        np.fill_diagonal(ref, 0.0)
        np.testing.assert_allclose(sm, ref, rtol=0, atol=1e-12)


# ---------------------------------------------------------------------------
# T9: the public seqtrate() matrix vs a live-R golden (TraMineR 2.2-14).
#
# The replication checks above pin the TRATE *cost* path; this section pins the
# public transition-rate matrix itself across the weighted / with.missing /
# ragged-void paths, generated by scripts/r_goldens/make_goldens_seqtrate.R.
# Mind the seqtrate-void caveat: with.missing=TRUE adds only the nr code to the
# state set, never the structural void.
# ---------------------------------------------------------------------------


def _seqtrate_state_cols(x, include_missing: bool):
    return [
        i
        for i in range(len(x.alphabet))
        if i != x.void_code and (include_missing or i != x.nr_code)
    ]


@pytest.mark.parametrize(
    "builder, tag, wm",
    [
        (lambda: _cs_fixture_a(), "A", False),
        (lambda: _cs_fixture_a(weighted=True), "Aw", False),
        (lambda: _cs_fixture_b(), "B_wm0", False),
        (lambda: _cs_fixture_b(), "B_wm1", True),
        (lambda: _cs_fixture_c(), "C_wm0", False),
        (lambda: _cs_fixture_c(), "C_wm1", True),
    ],
)
def test_public_seqtrate_matches_live_r(builder, tag, wm):
    from traminex.stats import seqtrate

    x = builder()
    out = np.asarray(
        seqtrate(x, weighted=True, include_missing=wm, as_dataframe=False), dtype=float
    )
    np.testing.assert_allclose(out, load_golden(f"seqtrate_{tag}"), rtol=0, atol=R_GOLDEN_ATOL)
    # The state set excludes void: dimension == real alphabet (+ nr when wm).
    assert out.shape[0] == len(_seqtrate_state_cols(x, wm))
