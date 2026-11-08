"""RP7a plot-series adapter (reshape-only; no numeric recomputation).

TraMineR separates *series computation* (``seqstatd`` / ``seqmeant`` / ``seqmodst``
/ ``seqsurv`` / ...) from *drawing* (``plot.stslist.*``). traminex already mirrors
the series half (R-golden'd in RP3-RP6). This module is the thin shared layer the
matplotlib renderers draw over: it reshapes the already-pinned series into the
bars / areas / steps / lines each plot needs, and resolves the index-plot sort
order. It computes **no** new numerics, so it carries no golden of its own -- its
inputs are golden'd and its reshaping is covered by the per-plot structure tests.

Sort resolvers mirror ``SOURCES/TraMineR/R/plot.stslist.R`` (the ``sortv`` block)
and ``seqplot.R``.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Literal

import numpy as np
import pandas as pd

from .._categorical_order import r_factor, r_order_keys
from .._rformat import r_match_arg
from ..sequences import STSList
from ..stats import seqstatd

__all__ = ["IndexSelection", "resolve_sort_order", "select_index_rows"]


@dataclass(frozen=True)
class IndexSelection:
    """Immutable row-selection result used by the index renderer."""

    source_rows: np.ndarray
    display_rows: np.ndarray
    weights: np.ndarray
    labels: tuple[object, ...]


def resolve_sort_order(x, sortv=None, *, weighted: bool = True) -> np.ndarray:
    """Row permutation for an index plot, faithful to ``plot.stslist`` ``sortv``.

    Parameters
    ----------
    x : STSList
        State sequence object (uses ``x.data``, the integer state-code matrix).
    sortv : None | "from.start" | "from.end" | array-like, optional
        - ``None``: identity order (no sorting).
        - ``"from.start"``: lexicographic order with the *first* time position as
          the primary key (``order`` over columns ``1..L`` in R).
        - ``"from.end"``: lexicographic order with the *last* time position as the
          primary key (``order`` over columns ``L..1`` in R).
        - array-like (one value per sequence): stable ``argsort`` of that vector
          (R ``order(sortv)``; factors are taken in their integer coding).
    weighted : bool, default True
        Accepted for signature parity with R; the order itself does not depend on
        weights (weights drive bar widths in the renderer, not the ordering).

    Returns
    -------
    np.ndarray
        Integer indices giving the plotted (top-to-bottom) row order.
    """
    # Cells rank as R's factor levels c(alphabet, nr, void), whatever codes
    # the markers are stored under (plot.stslist.R:62-67; TX-AUD-20261011-006).
    data = x._canonical_code_ranks(void_as_missing=False)[np.asarray(x.data)]
    n, length = data.shape
    if sortv is None:
        return np.arange(n)
    if isinstance(sortv, str):
        if length == 0:
            return np.arange(n)
        if sortv == "from.start":
            # primary key = first column -> it must be the LAST lexsort key
            keys = [data[:, j] for j in range(length - 1, -1, -1)]
            return np.lexsort(keys)
        if sortv == "from.end":
            # primary key = last column -> it must be the LAST lexsort key
            keys = [data[:, j] for j in range(length)]
            return np.lexsort(keys)
        raise ValueError("sortv string must be 'from.start' or 'from.end'")
    codes = getattr(sortv, "codes", None)
    if codes is None and isinstance(getattr(sortv, "dtype", None), pd.CategoricalDtype):
        # A pandas Series of category dtype keeps its codes under .cat: R orders a factor
        # by its level codes, so the Series is read like a Categorical (TX-AUD-20261019-019).
        codes = np.asarray(sortv.cat.codes)
    arr = np.asarray(codes if codes is not None else sortv)
    if arr.ndim != 1 or arr.shape[0] != n:
        raise ValueError(
            "sortv must contain one value for each sequence or be " "'from.start' / 'from.end'"
        )
    if codes is None and (arr.dtype.kind in "OUS" or not isinstance(sortv, np.ndarray)):
        order_keys = r_order_keys(list(sortv))
        if order_keys is not None:
            return np.argsort(order_keys, kind="stable")
    if codes is not None:
        # pandas categoricals encode missing factor values as -1, whereas R's
        # order(as.integer(factor)) keeps NA values at the end by default.
        arr = arr.astype(np.int64, copy=True)
        arr[arr < 0] = np.iinfo(np.int64).max
    return np.argsort(arr, kind="stable")


def _readonly(values, *, dtype=None) -> np.ndarray:
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def select_index_rows(
    x: STSList,
    *,
    sortv: object = None,
    sampling: int | None = None,
    sample_meth: Literal["systematic", "random"] = "systematic",
    weighted: bool = True,
    random_state: int | np.random.Generator | None = None,
) -> IndexSelection:
    """Sort and optionally subsample rows using ``plot.stslist`` mechanics.

    Sorting is applied to the full candidate set before selection. Systematic
    positions are the zero-based equivalent of R's truncating
    ``as.integer(seq(1, n, length.out=sampling))``. Random positions are drawn
    without replacement from an isolated NumPy generator, then sorted for
    display. ``sampling=None`` deliberately retains every row, following the R
    documentation rather than the R implementation's historical first-ten
    ``idxs`` default. ``sample_meth`` may be a unique prefix (``"sys"``), as R's
    ``match.arg`` takes it.
    """
    # plot.stslist.R:16 reads match.arg(sample.meth, c("systematic", "random")): an
    # exact name or a unique prefix (TX-REP-20261023-007).
    matched = r_match_arg(sample_meth, ("systematic", "random"))
    if matched is None:
        raise ValueError("sample_meth must be 'systematic' or 'random'")
    if sampling is not None:
        if isinstance(sampling, bool) or not isinstance(sampling, Integral):
            raise TypeError("sampling must be a positive integer or None")
        if sampling <= 0:
            raise ValueError("sampling must be a positive integer or None")
        sampling = int(sampling)

    candidate_rows = np.asarray(resolve_sort_order(x, sortv), dtype=int)
    candidate_count = candidate_rows.size

    raw_weights = getattr(x, "weights", None)
    if not weighted or raw_weights is None:
        all_weights = np.ones(x.nseq, dtype=float)
    else:
        all_weights = np.asarray(raw_weights, dtype=float)
        if all_weights.ndim != 1 or all_weights.shape[0] != x.nseq:
            raise ValueError("weights must contain one value for each sequence")
        if not np.all(np.isfinite(all_weights)):
            raise ValueError("weights must contain only finite values")
        if np.any(all_weights < 0):
            raise ValueError("weights must be non-negative")

    candidate_weights = all_weights[candidate_rows]
    if sampling is None or sampling >= candidate_count:
        positions = np.arange(candidate_count, dtype=int)
    elif matched == "systematic":
        positions = np.linspace(1, candidate_count, num=sampling).astype(int) - 1
    else:
        positive_count = int(np.count_nonzero(candidate_weights > 0))
        if candidate_weights.sum() <= 0:
            raise ValueError("weights must have a positive total for random sampling")
        if positive_count < sampling:
            raise ValueError("random sampling requires at least sampling positive weights")
        generator = np.random.default_rng(random_state)
        probabilities = candidate_weights / candidate_weights.sum()
        positions = np.sort(
            generator.choice(
                candidate_count,
                size=sampling,
                replace=False,
                p=probabilities,
            )
        )

    source_rows = candidate_rows[positions]
    ids = getattr(x, "ids", None)
    if ids is None:
        labels: tuple[object, ...] = tuple(int(row) + 1 for row in source_rows)
    else:
        labels = tuple(ids[int(row)] for row in source_rows)

    return IndexSelection(
        source_rows=_readonly(source_rows, dtype=int),
        display_rows=_readonly(np.asarray(x.data)[source_rows]),
        weights=_readonly(all_weights[source_rows], dtype=float),
        labels=labels,
    )


def _factor_panels(group) -> tuple[list[list[int]], list[str]]:
    """Rows and titles of R ``factor(group)``'s used levels.

    Shared by the grouped renderers (``seqedplot.R:53-66`` uses ``factor``).

    Numbers sort numerically and print as R prints them (``2``, not ``2.0``),
    text sorts in R's collation (ICU, batch-2 review M-1), a categorical keeps
    its declared level order, text mixed with numbers is refused, and missing
    values belong to no level (2026-09-22 audit S6-06, S6-07).
    """
    titles, codes = _r_factor_of(group)
    rows = [np.flatnonzero(codes == code).tolist() for code in range(len(titles))]
    return rows, titles


def _r_factor_of(group) -> tuple[list[str], np.ndarray]:
    """R ``factor(group)``'s used levels as R names them, and each case's level code (-1 missing).

    A categorical keeps its declared order; numbers equal to 15 significant
    digits are one level, as R's ``factor()`` matches the values' text
    (``r_factor``, TX-REP-20261022-006).
    """
    if isinstance(getattr(group, "dtype", None), pd.CategoricalDtype):
        categorical = pd.Categorical(group)
        labels, codes, _ = r_factor(
            np.asarray(categorical, dtype=object), list(categorical.categories)
        )
        return labels, codes
    labels, codes, _ = r_factor(list(group))
    return labels, codes


def transversal_entropy(x: STSList) -> np.ndarray:
    """Return the already-verified ``seqstatd`` entropy series as an array."""
    return np.asarray(seqstatd(x)["Entropy"], dtype=float)


def grouped_transversal_entropy(
    x: STSList, group
) -> tuple[tuple[object, ...], tuple[np.ndarray, ...]]:
    """Return non-missing factor groups and their ``seqstatd`` entropy curves.

    Groups are R's ``factor(group)`` levels (``seqplot-tentrop.R:10``): text in
    R's collation (ICU, batch-2 review M-1), numbers by value, a categorical in
    its declared order; text mixed with numbers is refused.
    """
    level_labels, codes = _r_factor_of(group)
    if len(codes) != x.nseq:
        raise ValueError("group must contain one value for each sequence")
    if not level_labels:
        raise ValueError("group must contain at least one non-missing group")

    labels: list[object] = []
    curves: list[np.ndarray] = []
    for code, label in enumerate(level_labels):
        rows = np.flatnonzero(codes == code)
        labels.append(label)
        subset = STSList(
            data=x.data[rows],
            alphabet=list(x.alphabet),
            labels=None if x.labels is None else list(x.labels),
            cpal=None if x.cpal is None else list(x.cpal),
            xtlab=None if x.xtlab is None else list(x.xtlab),
            xtstep=x.xtstep,
            tick_last=x.tick_last,
            nr=x.nr,
            void=x.void,
            nr_code=x.nr_code,
            void_code=x.void_code,
            weights=None if x.weights is None else np.asarray(x.weights)[rows],
            start=x.start,
            ids=None if x.ids is None else [x.ids[int(row)] for row in rows],
            missing_color=x.missing_color,
        )
        curves.append(transversal_entropy(subset))
    return tuple(labels), tuple(curves)
