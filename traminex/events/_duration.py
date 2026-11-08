"""Producing series for the event-duration plot ``seqedplot`` (BUILD-from-R).

``seqedplot`` (TraMineRextras ``seqedplot.R``) is a thin plotting layer over two
numeric series, one per non-ignored event:

- :func:`_seqedplot_survival` -- the (unweighted) Kaplan-Meier survival curve of
  the time to the event's **first** occurrence, with non-occurrence right-censored
  at the sequence's observation length (R ``survfit(Surv(time, status) ~ 1)``,
  ``seqedplot.R:201-235``);
- the hazard series (added in the hazard commit) -- the weighted mean number of
  occurrences per equal-width age period.

Both consume the event-age table :func:`~traminex.events._object.seqeage` and are
pinned to live R goldens (``test_r_goldens_seqedplot.py``); the renderer
(:func:`~traminex.visualization.plots.event_duration.seqedplot`) is
structure-tested against them. Pure Python, no compiled code.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..stats_p2 import _weighted_km
from ._object import SeqElist, seqeage, seqelength, seqeweight


def _kept_event_codes(eseq: SeqElist, ignore: Optional[Sequence[str]]) -> List[int]:
    """1-based codes of the events not in ``ignore`` (R ``which(!levels %in% ignore)``)."""
    drop = set() if ignore is None else {str(x) for x in ignore}
    return [i + 1 for i, lab in enumerate(eseq.dictionary) if lab not in drop]


def _subset_seqelist(eseq: SeqElist, idx: Sequence[int]) -> SeqElist:
    """A :class:`SeqElist` of the selected sequences (R ``seqe[gindex]``).

    Used to build per-group views for the renderer's ``group=`` faceting. The
    flat event views are rebuilt from the selected sequences so a counting
    helper reading them sees the group's events rather than an empty vector
    (cycle-8 audit F-16).
    """
    seqs = [eseq.sequences[i] for i in idx]
    flat_id: list[str] = []
    flat_time: list[float] = []
    flat_label: list[str] = []
    flat_groups: list[list[int]] = []
    for sequence in seqs:
        start = len(flat_id)
        for time, code in zip(sequence.times, sequence.codes):
            flat_id.append(str(sequence.id_code))
            flat_time.append(float(time))
            flat_label.append(eseq.dictionary[int(code) - 1])
        flat_groups.append(list(range(start, len(flat_id))))
    subset = SeqElist(
        seqs,
        eseq.dictionary,
        flat_id=flat_id,
        flat_time=np.asarray(flat_time, dtype=float),
        flat_label=flat_label,
    )
    # Each sequence's rows are recorded, so slicing the subset stays exact even
    # for a sequence without events (2026-09-22 audit S5-03).
    subset._flat_groups = flat_groups
    return subset


def _seqedplot_survival(
    eseq: SeqElist,
    ignore: Optional[Sequence[str]] = None,
    ages: Optional[Tuple[float, float]] = None,
) -> List[Tuple[str, np.ndarray]]:
    """Per-event survival series for ``seqedplot(type="survival")`` (BUILD-from-R).

    For each non-ignored event, take the age of its **first** occurrence (column 0
    of :func:`seqeage`, with -1 mapped to "never occurred" and -- when ``ages`` is
    given -- first-occurrence ages outside ``[ages[0], ages[1]]`` also treated as
    non-occurrences, mirroring ``seqedplot.R:131-134``). A sequence in which the
    event does not occur (within the window) is right-censored at its observation
    length ``seqelength`` (``seqedplot.R:205-208``). The curve is the **unweighted**
    Kaplan-Meier estimate ``survfit(Surv(time, status) ~ 1)`` -- ``seqeweight`` is
    deliberately not used here (``seqedplot.R:209``).

    Returns a list of ``(event_label, km)`` in kept-event (dictionary) order, where
    ``km`` is an ``(m, 5)`` array of ``time, n.risk, n.event, surv, std.err`` rows
    at the event times (R ``summary(survfit(...))``); an event that never occurs
    yields a ``(0, 5)`` array (a curve flat at survival 1).

    Raises ``ValueError`` when no observation time is set for any sequence (all
    ``seqelength == -1``), matching ``seqedplot.R:16-18``.
    """
    out: List[Tuple[str, np.ndarray]] = []
    for label, time, status in _seqedplot_times(eseq, ignore, ages):
        km = _weighted_km(time, status, np.ones_like(time, dtype=float))
        arr = np.array(km, dtype=float).reshape(-1, 5) if km else np.empty((0, 5), dtype=float)
        out.append((label, arr))
    return out


def _seqedplot_times(
    eseq: SeqElist,
    ignore: Optional[Sequence[str]] = None,
    ages: Optional[Tuple[float, float]] = None,
) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    """Per kept event, survfit's ``(time, status)``: the first occurrence or the censoring time."""
    endobs = seqelength(eseq)
    if np.all(endobs == -1):
        raise ValueError("You should set observation time to use type='survival' (see seqelength).")
    out: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for code in _kept_event_codes(eseq, ignore):
        am = seqeage(eseq, code).astype(float)
        am[am == -1] = np.nan
        if ages is not None:
            am[am < ages[0]] = np.nan
            am[am > ages[1]] = np.nan
        time = am[:, 0].copy()
        status = ~np.isnan(time)
        time[~status] = endobs[~status]
        out.append((eseq.dictionary[code - 1], time, status))
    return out


def _seqedplot_follow_up(
    eseq: SeqElist,
    ignore: Optional[Sequence[str]] = None,
    ages: Optional[Tuple[float, float]] = None,
) -> dict:
    """Each kept event's last survfit time, event or censoring, which plot.survfit draws to."""
    return {
        label: float(np.max(time)) if time.size else 0.0
        for label, time, _ in _seqedplot_times(eseq, ignore, ages)
    }


def _seqedplot_ages_range(
    eseq: SeqElist, ignore: Optional[Sequence[str]] = None
) -> Tuple[float, float]:
    """Global ``(minage, maxage)`` over all non-ignored events (R ``seqedplot.R:135-140``).

    The default ``ages`` window for both series: the smallest and largest occurrence
    age across every kept event and every sequence (``-1`` non-occurrences ignored).
    For grouped plots the renderer computes this once over the full set so all groups
    share the same cut-points.
    """
    vals = []
    for code in _kept_event_codes(eseq, ignore):
        am = seqeage(eseq, code).astype(float)
        am[am == -1] = np.nan
        vals.append(am.ravel())
    allv = np.concatenate(vals) if vals else np.array([np.nan])
    return float(np.nanmin(allv)), float(np.nanmax(allv))


def _seqedplot_cutpoints(ages: Tuple[float, float], breaks: int) -> Tuple[np.ndarray, np.ndarray]:
    """Bin break-points from ``cut(ages, breaks=N)`` (R ``seqedplot.R:154-158``).

    Faithful reproduction of ``base::cut.default`` for a single ``breaks`` count,
    followed by seqedplot's parse of the interval labels back to numbers:

    1. take ``N+1`` equally spaced points over ``[ages[0], ages[1]]`` and push the
       two ends out by ``dx/1000`` (``dx`` = the range; the degenerate ``dx == 0``
       case uses ``|ages[0]|`` or 1) -- so the data minimum/maximum fall strictly
       inside the first/last bin;
    2. round each break to the **fewest significant digits** (``dig.lab`` from 3 up
       to 12) that keep consecutive breaks distinct, exactly as ``cut`` formats its
       labels with ``formatC(..., format="g")`` -- the rounded values *are* the
       break-points seqedplot then re-parses.

    Returns ``(cutpoints, onx)``: the ``N+1`` break-points and the ``N`` upper edges
    (the abscissae seqedplot plots, ``cutpoints[1:]``).
    """
    lo, hi = float(ages[0]), float(ages[1])
    nb = int(breaks) + 1
    dx = hi - lo
    raw = np.linspace(lo, hi, nb)
    if dx == 0:
        d = abs(lo) if lo != 0 else 1.0
        raw = np.linspace(lo - d / 1000.0, hi + d / 1000.0, nb)
    else:
        raw[0] = lo - dx / 1000.0
        raw[-1] = hi + dx / 1000.0

    ch: List[str] = [_format_g(b, 3) for b in raw]
    for dig in range(3, 13):
        ch = [_format_g(b, dig) for b in raw]
        if all(ch[i] != ch[i + 1] for i in range(len(ch) - 1)):
            break
    cutpoints = np.array([float(c) for c in ch], dtype=float)
    onx = cutpoints[1:].copy()
    return cutpoints, onx


def _format_g(x: float, dig: int) -> str:
    """``formatC(x, digits=dig, width=1L)`` -- C ``printf`` ``%g`` with ``dig`` sig. digits."""
    return "%.*g" % (dig, float(x))


def _seqedplot_hazard(
    eseq: SeqElist,
    ignore: Optional[Sequence[str]] = None,
    ages: Optional[Tuple[float, float]] = None,
    breaks: int = 20,
) -> Tuple[np.ndarray, List[Tuple[str, np.ndarray]]]:
    """Per-event hazard series for ``seqedplot(type="hazard")`` (BUILD-from-R).

    For each non-ignored event, bin **every** occurrence age (all columns of
    :func:`seqeage`) into the equal-width age periods, sum the per-sequence weights
    ``seqeweight`` falling in each period and divide by the number of sequences --
    the weighted mean number of occurrences per period (``seqedplot.R:159-172``).
    Occurrence ages outside ``[ages[0], ages[1]]`` (and the ``-1`` non-occurrences)
    are dropped, matching R's ``cut`` of the age matrix.

    ``ages`` defaults to the global range over the given sequences
    (:func:`_seqedplot_ages_range`); the renderer passes the full-set range so every
    group shares the same cut-points. The bins are right-closed ``(lo, hi]`` exactly
    as R's ``cut(..., right=TRUE)``.

    Returns ``(onx, series)`` where ``onx`` is the ``N`` period upper edges and
    ``series`` is a list of ``(event_label, ony)`` in kept-event order, each ``ony``
    a length-``N`` array of per-period means.
    """
    if ages is None:
        ages = _seqedplot_ages_range(eseq, ignore)
    cutpoints, onx = _seqedplot_cutpoints(ages, breaks)
    n_periods = len(onx)
    weights = seqeweight(eseq)
    n_seq = len(eseq)

    series: List[Tuple[str, np.ndarray]] = []
    for code in _kept_event_codes(eseq, ignore):
        am = seqeage(eseq, code).astype(float)
        am[am == -1] = np.nan
        am[am < ages[0]] = np.nan
        am[am > ages[1]] = np.nan
        ony = np.zeros(n_periods, dtype=float)
        for i in range(am.shape[0]):
            wi = weights[i]
            for v in am[i]:
                if not np.isnan(v):
                    # right-closed bin (cutpoints[b], cutpoints[b+1]] -> index b
                    b = int(np.searchsorted(cutpoints, v, side="left")) - 1
                    if 0 <= b < n_periods:
                        ony[b] += wi
        ony /= n_seq
        series.append((eseq.dictionary[code - 1], ony))
    return onx, series
