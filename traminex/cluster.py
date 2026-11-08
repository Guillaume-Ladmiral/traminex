from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Literal, Optional

import numpy as np
import pandas as pd

from ._categorical_order import r_factor
from ._deprecated import MISSING, resolve_deprecated_alias
from .sequences import STSList, _require_stslist


@dataclass
class ClusterResult:
    labels_: np.ndarray
    medoid_index: Optional[int] = None


def medoid_index(D: np.ndarray) -> int:
    """Find the medoid (most central point) in a distance matrix.

    The medoid is the observation with the smallest sum of distances to all
    other observations, making it the most representative point.

    Parameters
    ----------
    D : np.ndarray
        Square distance matrix of shape (n, n)

    Returns
    -------
    int
        Index of the medoid observation (0-based)
    """
    s = D.sum(axis=1)
    return int(np.argmin(s))


def cluster_agglomerative(D: np.ndarray, n_clusters: int) -> ClusterResult:
    """Perform agglomerative hierarchical clustering on distance matrix.

    Uses average linkage with precomputed distances. Ward linkage is not
    applicable with precomputed distances.

    Parameters
    ----------
    D : np.ndarray
        Square distance matrix of shape (n, n)
    n_clusters : int
        Number of clusters to form

    Returns
    -------
    ClusterResult
        Object with attributes:
        - labels_: cluster assignment for each observation
        - medoid_index: zero-based row of the overall medoid

    Raises
    ------
    RuntimeError
        If scikit-learn is not installed
    """
    # Imported here: importing scikit-learn sets KMP_* variables in the process
    # environment (TX-AUD-20261011-034).
    try:
        from sklearn.cluster import AgglomerativeClustering
    except ImportError:
        raise RuntimeError("scikit-learn is required for clustering") from None
    model = AgglomerativeClustering(
        n_clusters=n_clusters,
        metric="precomputed",
        linkage="average",
    )
    labels = model.fit_predict(D)
    return ClusterResult(labels_=labels, medoid_index=medoid_index(D))


def _label_value(label: object) -> object:
    """A cluster label as given, as a Python scalar (an integer label stays an ``int``)."""
    return label.item() if isinstance(label, np.generic) else label


@_require_stslist
def cluster_summary(x, labels: np.ndarray, D: Optional[np.ndarray] = None) -> pd.DataFrame:
    """Summarize clusters: size, medoid (optional), and per-state mean durations.

    If D is provided, computes medoid per cluster using within-cluster distances.
    Returns a DataFrame with one row per cluster and columns:
      - cluster, size, medoid_index (zero-based row of the cluster's medoid in
        ``x``), and state_<state> columns with mean durations
    """
    from .stats import seqistatd

    # The per-state column set is keyed on ``len(x.alphabet)``, so a declared
    # *external* marker code (outside ``range(len(alphabet))``) produced fewer
    # columns than the identical in-alphabet twin. Canonicalizing at entry makes
    # both representations emit the same columns (cycle-10 N-14).
    if isinstance(x, STSList):
        x = x._canonicalized()

    labels = np.asarray(labels)
    K = len(x.alphabet)
    durations = seqistatd(x, include_missing=False)
    rows = []
    for c in np.unique(labels):
        idx = np.where(labels == c)[0]
        size = idx.size
        med = None
        if D is not None and size > 0:
            sub = D[np.ix_(idx, idx)]
            mlocal = medoid_index(sub)
            med = int(idx[mlocal])
        row: Dict[str, object] = {
            "cluster": _label_value(c),
            "size": int(size),
            "medoid_index": (None if med is None else med),
        }
        mean_dur = durations[idx].mean(axis=0) if size > 0 else np.zeros(K)
        for k in range(K):
            row[f"state_{x.alphabet[k]}"] = float(mean_dur[k])
        rows.append(row)
    return pd.DataFrame(rows).sort_values("cluster").reset_index(drop=True)


@_require_stslist
def cluster_representatives(
    x, labels: np.ndarray, D: np.ndarray, *, sep: str = "-"
) -> pd.DataFrame:
    """Return representative (medoid) sequence per cluster with its index and string.

    Columns: cluster, rep_index (zero-based row of the medoid in ``x``), rep_sequence
    """
    labels = np.asarray(labels)
    reps: List[Dict[str, object]] = []
    strings = x.to_strings(sep=sep)
    for c in np.unique(labels):
        idx = np.where(labels == c)[0]
        if idx.size == 0:
            continue
        sub = D[np.ix_(idx, idx)]
        mlocal = medoid_index(sub)
        gidx = int(idx[mlocal])
        reps.append({"cluster": _label_value(c), "rep_index": gidx, "rep_sequence": strings[gidx]})
    return pd.DataFrame(reps).sort_values("cluster").reset_index(drop=True)


@_require_stslist
def seqrep(
    x,
    *,
    diss: Optional[np.ndarray] = None,
    criterion: Literal["density", "freq", "dist", "prob", "mscore", "random"] = "density",
    score: Optional[np.ndarray] = None,
    decreasing: bool = True,
    coverage: float = 0.25,
    nrep: Optional[int] = None,
    pradius: float = 0.10,
    dmax: Optional[float] = None,
    weighted: bool = True,
    trep=MISSING,
    tsim=MISSING,
    dist_matrix=MISSING,
    random_state: int | np.random.Generator | None = None,
    **dist_kwargs,
) -> dict:
    """Representative sequences of a set of sequences — R-faithful.

    Port of TraMineR 2.2-14 ``seqrep`` (``SOURCES/TraMineR/R/seqrep.R``): build
    (or accept) the pairwise dissimilarity matrix, then delegate to
    :func:`dissrep` to extract the representative set, returning the
    representative sub-sequences together with the dissrep quality statistics.
    Anchored to live R by ``tests/test_r_goldens_dissrep.py``.

    Parameters
    ----------
    x : STSList
        State sequence object (``seqdef``).
    diss : np.ndarray, optional
        Pairwise dissimilarity matrix. If ``None`` it is computed from ``x`` via
        :func:`traminex.seqdist` (pass the method etc. as ``**dist_kwargs``).
    criterion : {"density", "freq", "dist", "prob", "mscore", "random"}, default "density"
        Representativeness criterion. ``prob`` ranks by sequence likelihood;
        ``mscore`` ranks by the summed cross-sectional frequency of each
        trajectory's observed states.
    score, decreasing, coverage, nrep, pradius, dmax
        Passed through to :func:`dissrep` (R defaults: ``coverage=0.25``,
        ``pradius=0.10``).
    weighted : bool, default True
        Use the sequence weights attached to ``x`` when present.
    random_state : int or numpy.random.Generator, optional
        Random-number source for ``criterion="random"``, passed to
        :func:`dissrep` (R draws from its global stream).
    **dist_kwargs
        Forwarded to :func:`traminex.seqdist` when ``diss`` is ``None``. With
        ``diss`` given they have nothing to act on and raise ``TypeError`` (R
        drops them silently; TX-AUD-20261010-015).

    Returns
    -------
    dict
        ``rep_indices`` (0-based indices of the representatives), the
        ``representatives`` code rows, and the dissrep ``statistics`` /
        ``quality`` / ``distances`` / ``rep_group``.
    """
    from .diss_p1 import dissrep

    coverage = resolve_deprecated_alias(
        coverage,
        trep,
        canonical_name="coverage",
        alias_name="trep",
        function_name="seqrep",
        default=0.25,
    )
    pradius = resolve_deprecated_alias(
        pradius,
        tsim,
        canonical_name="pradius",
        alias_name="tsim",
        function_name="seqrep",
        default=0.10,
    )
    diss = resolve_deprecated_alias(
        diss,
        dist_matrix,
        canonical_name="diss",
        alias_name="dist_matrix",
        function_name="seqrep",
        default=None,
    )

    seqmat = np.asarray(x.data)
    weights = None
    if weighted:
        w = getattr(x, "weights", None)
        if w is not None:
            weights = np.asarray(w, dtype=float)

    if diss is None:
        from .distances import seqdist

        diss = seqdist(x, **dist_kwargs)
    elif dist_kwargs:
        raise TypeError(
            "seqrep got keyword(s) it cannot use when diss is given: " f"{sorted(dist_kwargs)}"
        )
    diss = np.asarray(diss, dtype=float)

    resolved_score = None if score is None else np.asarray(score, dtype=float)
    if resolved_score is None and criterion == "prob":
        from .stats_complexity import seqlogp

        # R seqrep calls seqlogp(seqdata) without forwarding its own weighted
        # flag. The attached sequence weights therefore fit this score in both
        # the weighted and unweighted representative-selection routes.
        resolved_score = np.asarray(seqlogp(x), dtype=float)
        decreasing = False
    elif resolved_score is None and criterion == "mscore":
        from .stats_basic import seqstatd

        # Exact TraMineR.mscore: sum (not mean) the cross-sectional frequency
        # of each regular state at each rectangular position. Missing and void
        # cells do not belong to statelist and therefore add nothing.
        frequencies = np.asarray(seqstatd(x)["Frequencies"], dtype=float)
        marker_codes = {code for code in (x.nr_code, x.void_code) if code is not None}
        resolved_score = np.zeros(x.nseq, dtype=float)
        for i, row in enumerate(seqmat):
            for time, code in enumerate(row):
                code = int(code)
                if code not in marker_codes and code < frequencies.shape[1]:
                    resolved_score[i] += frequencies[time, code]
        decreasing = True

    # dissrep's own criterion chain computes (and would ignore any passed)
    # scores for the density/dist/freq routes, so seqrep no longer duplicates
    # that work here (cycle-8 audit F-45); it precomputes only the prob/mscore
    # scores dissrep cannot derive itself.
    rep = dissrep(
        diss,
        criterion=criterion,  # type: ignore[arg-type]
        score=resolved_score,
        decreasing=decreasing,
        coverage=coverage,
        nrep=nrep,
        pradius=pradius,
        dmax=dmax,
        weights=weights,
        random_state=random_state,
    )
    # dissrep owns the effective score for its density/frequency/distance and
    # random routes.  Preserve that public result instead of returning a scalar
    # NaN when seqrep did not precompute prob/mscore scores itself.
    if resolved_score is None:
        resolved_score = np.asarray(rep["scores"], dtype=float)
    idx = rep["representatives"]
    distance_matrix = diss[:, idx].astype(float, copy=True)
    if idx.size > 1:
        groups = np.asarray(rep["rep_group"], dtype=int)
        for row, group in enumerate(groups):
            distance_matrix[row, np.arange(idx.size) != group] = np.nan
    distance_to_rep = np.nanmin(distance_matrix, axis=1)
    idx_rep = [
        np.flatnonzero((np.asarray(rep["rep_group"], dtype=int) == group) & (distance_to_rep == 0))
        for group in range(idx.size)
    ]
    metadata = dict(rep["metadata"])
    metadata.update(
        {
            "criterion": criterion,
            "weighted": bool(
                weighted and weights is not None and not np.all(np.asarray(weights) == 1)
            ),
        }
    )
    return {
        "rep_indices": idx,
        "representatives": seqmat[idx],
        "scores": np.asarray(resolved_score, dtype=float),
        "statistics": rep["statistics"],
        "quality": rep["quality"],
        "distances": rep["distances"],
        "distance_matrix": distance_matrix,
        "rep_group": rep["rep_group"],
        "idx_rep": idx_rep,
        "metadata": metadata,
    }


def seqrep_grp(
    x,
    group=None,
    *,
    diss: Optional[np.ndarray] = None,
    ret: str = "stat",
    criterion: Literal["density", "freq", "dist", "prob", "mscore", "random"] = "density",
    coverage: float = 0.25,
    weighted: bool = True,
    with_missing: bool = False,
    mdis=MISSING,
    **dist_kwargs,
):
    """Group-wise representative sequences with a shared ``dmax`` (R ``seqrep.grp``).

    Faithful port of TraMineRextras 0.6.9 ``seqrep.grp`` (GPL): runs
    :func:`seqrep` separately on each group, using a single shared
    ``dmax = max(diss)`` across groups so the per-group quality statistics are
    comparable. Each group must contain at least two cases.

    Parameters
    ----------
    x : STSList
        State sequence object.
    group : array-like or None, default None
        None puts all cases in one group. Group membership (one entry per sequence). Levels are taken in sorted
        order, as R's ``factor`` does. A sequence whose group is missing
        (``None``, ``NaN``, ``pd.NA``) belongs to no group and is left out, as
        in R.
    diss : np.ndarray, optional
        Full pairwise dissimilarity matrix. If ``None`` it is computed from ``x``
        via :func:`traminex.seqdist` (pass ``method`` etc. as ``**dist_kwargs``).
    ret : {"stat", "rep", "both"}, default "stat"
        ``"stat"`` returns the per-group quality statistics, ``"rep"`` the
        per-group :func:`seqrep` results (their ``rep_indices`` are zero-based
        positions within the group's own sequences), ``"both"`` a ``(rep, stat)``
        tuple.
    criterion, coverage, weighted
        Forwarded to :func:`seqrep` (R defaults: ``coverage=0.25``).
    **dist_kwargs
        Distance options when computing ``diss``; ``nrep``, ``score``,
        ``decreasing``, ``pradius``, ``trep`` and ``tsim`` are forwarded to
        representative selection. ``dmax`` is shared and cannot be overridden.

    Returns
    -------
    dict | tuple
        Keyed by the (stringified) group level.
    """
    diss = resolve_deprecated_alias(
        diss,
        mdis,
        canonical_name="diss",
        alias_name="mdis",
        function_name="seqrep_grp",
        default=None,
    )
    if not isinstance(with_missing, (bool, np.bool_)):
        raise TypeError("with_missing must be a boolean")
    if ret not in ("stat", "rep", "both"):
        raise ValueError("ret should be one of 'stat', 'rep', 'both'")
    if not isinstance(x, STSList):
        raise TypeError("x must be a state sequence object")
    representative_keys = {"nrep", "score", "decreasing", "pradius", "trep", "tsim", "random_state"}
    representative_options = {
        key: dist_kwargs.pop(key) for key in list(dist_kwargs) if key in representative_keys
    }
    if "dmax" in dist_kwargs:
        raise TypeError("dmax is computed from the full distance matrix and cannot be overridden")
    # Check accepted distance keys even when a supplied matrix bypasses seqdist.
    from .distances._validation import validate_pairwise_kwargs

    distance_controls = {
        "method",
        "weighted",
        "full_matrix",
        "with_missing",
        "void",
        "deduplicate",
        "condensed",
        "refseq",
        "ref_sets",
    }
    validate_pairwise_kwargs(
        {"OMOPT": "OM", "LCSOPT": "LCS"}.get(
            str(dist_kwargs.get("method", "OM")).upper(), dist_kwargs.get("method", "OM")
        ),
        {key: value for key, value in dist_kwargs.items() if key not in distance_controls},
        dispatcher=True,
        layer="seqdist",
    )
    group = np.ones(x.nseq, dtype=int) if group is None else np.asarray(group)
    if group.ndim != 1 or group.size != x.nseq:
        raise ValueError("group must be one-dimensional with one entry per sequence")
    score = representative_options.get("score")
    if score is not None:
        score = np.asarray(score, dtype=float)
        if score.ndim != 1 or score.size != x.nseq:
            raise ValueError("score must have one value per sequence")
    # A case whose group is missing belongs to no group: R's xtabs(~grp),
    # levels(factor(grp)) and which(grp == level) all leave NA out
    # (TraMineRextras seqrep-grp.R:15-45).
    observed_cases = np.flatnonzero(~pd.isna(group))
    observed = group[observed_cases]
    # R: levels(factor(grp)) -- text in R's collation (TX-AUD-20261010-039), numbers
    # named as R prints them and grouped to 15 significant digits (TX-AUD-20261019-046).
    # R's names carry a leading space (TraMineR:::group pastes "" and the level,
    # seqrep-grp.R:16); the keys here are the level names themselves.
    labels, codes, _ = r_factor(observed)
    for code in range(len(labels)):
        if int(np.sum(codes == code)) < 2:
            raise ValueError("Each group must have 2 or more cases. At least one group has only 1.")

    if diss is None:
        from .distances import seqdist

        diss = seqdist(x, with_missing=bool(with_missing), **dist_kwargs)
    diss = np.asarray(diss, dtype=float)
    if diss.shape != (x.nseq, x.nseq):
        raise ValueError("diss must be a square matrix with one row per sequence")
    dmax = float(diss.max())  # shared across groups

    stats: Dict[str, object] = {}
    reps: Dict[str, object] = {}
    for code, key in enumerate(labels):
        ig = observed_cases[codes == code]
        sub = x.copy()
        sub.data = x.data[ig]
        if getattr(x, "weights", None) is not None:
            sub.weights = np.asarray(x.weights, dtype=float)[ig]
        sub_diss = diss[np.ix_(ig, ig)]
        options = dict(representative_options)
        if score is not None:
            options["score"] = score[ig]
        r = seqrep(
            sub,
            diss=sub_diss,
            criterion=criterion,
            coverage=coverage,
            dmax=dmax,
            weighted=weighted,
            **options,
        )
        stats[key] = r["statistics"]
        reps[key] = r

    if ret == "stat":
        return stats
    if ret == "rep":
        return reps
    return reps, stats
