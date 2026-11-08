"""Exploratory sequence-tree helpers: sequence features fed to a scikit-learn tree.

.. warning:: **These are NOT the TraMineR ``seqtree``.** ``seqtree_fit``
   computes per-sequence features (time in each state, transitions between
   states, entropy and number of transitions) and fits an unweighted
   scikit-learn ``DecisionTreeClassifier`` to a label. It is not Studer &
   Ritschard's weighted-inertia regression tree on a dissimilarity matrix.
   The R-faithful, weighted ``seqtree`` lives in :mod:`traminex.seqtree`
   (see DUPLICATES.md A8 / migration plan Phase 2.4).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from ._deprecated import _caller_stacklevel
from .disstree_explore import _import_pyplot
from .sequences import STSList, _real_state_codes, _require_stslist
from .stats import seqistatd


def _features_from_sequences(x: STSList, kind: str = "durations+transitions+summary"):
    # Feature dimensions, labels and entropy normalization must use the same
    # codebook as the statistics, which already fold declared external markers.
    # The features are over the alphabet's states: a missing or void cell has
    # no column, joins no transition and is outside the entropy
    # (TX-AUD-20261010-081).
    x = x._canonicalized()
    real = _real_state_codes(x)
    states = [x.alphabet[code] for code in real]
    K = len(real)
    # Each cell's position among the states, -1 for a marker cell.
    position = np.full(max(len(x.alphabet), int(x.data.max(initial=0)) + 1), -1, dtype=np.int64)
    position[real] = np.arange(K)
    cells = position[np.asarray(x.data, dtype=np.int64)]
    feats: List[np.ndarray] = []
    names: List[str] = []
    if "durations" in kind:
        D = seqistatd(x, include_missing=False).astype(float)[:, real]
        feats.append(D)
        names.extend([f"dur_{s}" for s in states])
    if "transitions" in kind:
        # Per-sequence transition counts between two states (flattened K×K)
        transition_rows: List[np.ndarray] = []
        for i in range(x.nseq):
            a = cells[i, :-1]
            b = cells[i, 1:]
            both = (a >= 0) & (b >= 0)
            idx = a[both] * K + b[both]
            counts = np.bincount(idx, minlength=K * K).astype(float).reshape(K, K)
            transition_rows.append(counts.ravel())
        transition_matrix = np.vstack(transition_rows)
        feats.append(transition_matrix)
        names.extend([f"tr_{i}_{j}" for i in states for j in states])
    if "summary" in kind:
        # Per-sequence summaries: normalized entropy and transition-count proxy
        ent_vec = []
        for i in range(x.nseq):
            row = cells[i]
            state_counts = np.bincount(row[row >= 0], minlength=K).astype(float)
            s = state_counts.sum()
            if s > 0:
                # A zero share contributes 0 (TX-AUD-20261010-095).
                p = state_counts[state_counts > 0] / s
                H = float(-np.sum(p * np.log(p)))
                Hn = H / (np.log(K) if K > 1 else 1.0)
            else:
                Hn = 0.0
            ent_vec.append(Hn)
        ent = np.array(ent_vec, dtype=float).reshape(-1, 1)
        # Changes between two adjacent state cells, over the width seqlen - 1; a
        # marker cell joins no transition (batch-12 review M-2). seqtransn_count
        # counts as R's seqtransn instead (TX-AUD-20261011-068).
        before, after = cells[:, :-1], cells[:, 1:]
        changes = ((before >= 0) & (after >= 0) & (before != after)).sum(axis=1)
        turb = (changes / max(x.seqlen - 1, 1)).astype(float).reshape(-1, 1)
        feats.append(ent)
        feats.append(turb)
        names.extend(["entropy_norm", "transition_count_norm"])
    X = np.hstack(feats) if feats else np.zeros((x.nseq, 0))
    return X, names


def _sklearn_tree() -> Any:
    """scikit-learn's tree module, imported when a helper needs it, or None.

    Importing scikit-learn sets KMP_* variables in the process environment, so
    the module imports it only where it is used (TX-AUD-20261011-034).
    """
    try:
        from sklearn import tree
    except ImportError:
        return None
    return tree


@dataclass
class SeqTreeResult:
    model: Any
    feature_names: List[str]
    labels: np.ndarray
    leaf_index: np.ndarray  # per-sequence leaf id


@_require_stslist
def seqtree_fit(
    x: STSList,
    y,
    *,
    features: str = "durations+transitions+summary",
    max_depth: Optional[int] = None,
    min_samples_leaf: int = 2,
    random_state: int = 0,
) -> SeqTreeResult:
    """Fit a scikit-learn decision tree to a label from per-sequence features (exploratory).

    Not the TraMineR ``seqtree`` (no weights, no inertia splits on a
    dissimilarity matrix); for the R-faithful weighted regression tree use
    :func:`traminex.seqtree.seqtree`.

    features: string containing tokens among durations, transitions, summary
    """
    sk_tree = _sklearn_tree()
    if sk_tree is None:
        raise RuntimeError("scikit-learn is required for seqtree_fit")
    X, names = _features_from_sequences(x, kind=features)
    clf = sk_tree.DecisionTreeClassifier(
        max_depth=max_depth, min_samples_leaf=min_samples_leaf, random_state=random_state
    )
    labels = np.asarray(y)
    clf.fit(X, labels)
    leaf = clf.apply(X)
    return SeqTreeResult(model=clf, feature_names=names, labels=labels, leaf_index=leaf)


def seqtree_leaf_summary(
    x: STSList | SeqTreeResult | None = None,
    result: SeqTreeResult | None = None,
) -> pd.DataFrame:
    """Return a DataFrame of leaf identifiers and sizes.

    Pass the fitted :class:`SeqTreeResult` as the sole argument. The historical
    ``seqtree_leaf_summary(x, result)`` form remains available with a
    deprecation warning; ``x`` was never used to calculate the summary.
    """
    if result is None:
        if not isinstance(x, SeqTreeResult):
            raise TypeError("seqtree_leaf_summary requires a SeqTreeResult")
        resolved = x
    elif x is None:
        resolved = result
    else:
        warnings.warn(
            "seqtree_leaf_summary(x, result) is deprecated; pass only result",
            DeprecationWarning,
            stacklevel=_caller_stacklevel(),
        )
        resolved = result

    uniques, counts = np.unique(resolved.leaf_index, return_counts=True)
    return (
        pd.DataFrame({"leaf": uniques, "size": counts}).sort_values("leaf").reset_index(drop=True)
    )


def seqtree_plot(result: SeqTreeResult, ax=None):
    """Plot the fitted decision tree (requires sklearn)."""
    if not isinstance(result, SeqTreeResult):
        raise TypeError("result must be a SeqTreeResult")
    sk_tree = _sklearn_tree()
    if sk_tree is None:
        raise RuntimeError("scikit-learn is required for seqtree_plot")
    if ax is None:
        plt = _import_pyplot("seqtree_plot")
        fig, ax = plt.subplots(figsize=(10, 6))
    sk_tree.plot_tree(
        result.model,
        feature_names=result.feature_names,
        class_names=[str(c) for c in np.unique(result.labels)],
        filled=True,
        ax=ax,
    )
    return ax


@_require_stslist
def seqtree_cv_fit(
    x: STSList,
    y,
    *,
    features: str = "durations+transitions+summary",
    param_grid: Optional[Dict[str, Any]] = None,
    cv: int = 5,
    scoring: Optional[str] = None,
    random_state: int = 0,
):
    """Cross-validated tree fit; returns best estimator and CV results."""
    sk_tree = _sklearn_tree()
    if sk_tree is None:
        raise RuntimeError("scikit-learn is required for seqtree_cv_fit")
    from sklearn.model_selection import GridSearchCV

    X, names = _features_from_sequences(x, kind=features)
    clf = sk_tree.DecisionTreeClassifier(random_state=random_state)
    if param_grid is None:
        param_grid = {"max_depth": [2, 3, 5, None], "min_samples_leaf": [1, 2, 4]}
    gs = GridSearchCV(clf, param_grid=param_grid, cv=cv, scoring=scoring)
    gs.fit(X, np.asarray(y))
    return gs.best_estimator_, names, gs.cv_results_


def seqtree_predict(
    result: SeqTreeResult, x: STSList, *, features: str = "durations+transitions+summary"
) -> np.ndarray:
    """Predict labels for new sequences using a fitted seqtree result.

    Recomputes features for x using the provided features string. If the
    resulting feature names differ in order from result.feature_names but contain
    the same set, columns are reordered to match. If a required feature is
    missing, raises ValueError.
    """
    if not isinstance(result, SeqTreeResult):
        raise TypeError("result must be a SeqTreeResult")
    if not isinstance(x, STSList):
        raise TypeError("x must be an STSList state sequence object")
    X_new, names_new = _features_from_sequences(x, kind=features)
    if names_new != result.feature_names:
        name_to_idx = {n: i for i, n in enumerate(names_new)}
        try:
            idx = [name_to_idx[n] for n in result.feature_names]
            X_new = X_new[:, idx]
        except KeyError as e:
            raise ValueError(
                "seqtree_predict: feature mismatch between fitted model and new data"
            ) from e
    return result.model.predict(X_new)
