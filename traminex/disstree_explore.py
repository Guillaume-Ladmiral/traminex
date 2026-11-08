"""Exploratory dissimilarity-tree helpers: an MDS embedding fed to a scikit-learn tree.

.. warning:: **These are NOT the TraMineR ``disstree``.** ``disstree_fit``
   embeds the distance matrix with classical MDS and fits an unweighted
   scikit-learn ``DecisionTreeClassifier`` to a label on the coordinates. It
   is not Studer & Ritschard's weighted-inertia regression tree. The
   R-faithful, weighted ``disstree`` lives in :mod:`traminex.disstree` (see
   DUPLICATES.md A8 / migration plan Phase 2.4).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
import pandas as pd


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


def _import_pyplot(feature: str) -> Any:
    """Import matplotlib.pyplot for ``feature`` or raise an ImportError naming the extra.

    The wording is the package facade's for its lazy plotting modules
    (TX-AUD-20261011-033); a bare ModuleNotFoundError named no extra.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        from . import _OPTIONAL_DEPENDENCY_EXTRAS

        extra = _OPTIONAL_DEPENDENCY_EXTRAS["matplotlib"]
        raise ImportError(
            f"matplotlib is required for {feature}; install the '{extra}' extra "
            f"(pip install 'traminex[{extra}]')"
        ) from exc
    return plt


def _cmdscale_fit(
    D: np.ndarray, n_components: int = 5
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Classical MDS fit with quantities needed for out-of-sample projection."""
    D = np.asarray(D, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("D must be a square distance matrix")
    n = D.shape[0]
    D2 = D**2
    row_mean_d2 = D2.mean(axis=1)
    grand_mean_d2 = float(D2.mean()) if D2.size else 0.0
    # Centering matrix
    J = np.eye(n) - np.ones((n, n)) / n
    # Double centered Gram from distances
    B = -0.5 * J.dot(D2).dot(J)
    # Eigen decomposition
    vals, vecs = np.linalg.eigh(B)
    idx = np.argsort(vals)[::-1]
    vals = vals[idx]
    vecs = vecs[:, idx]
    # Keep positive components
    pos = vals > 0
    eigvals = vals[pos][:n_components]
    eigvecs = vecs[:, pos][:, :n_components]
    if eigvals.size == 0:
        X = np.zeros((n, n_components))
        return X, eigvecs, eigvals, row_mean_d2, grand_mean_d2
    L = np.diag(np.sqrt(eigvals))
    X = eigvecs.dot(L)
    # pad if needed
    if X.shape[1] < n_components:
        pad = np.zeros((n, n_components - X.shape[1]))
        X = np.hstack([X, pad])
    return X, eigvecs, eigvals, row_mean_d2, grand_mean_d2


@dataclass
class DissTreeResult:
    model: Any
    coords: np.ndarray
    labels: np.ndarray
    # Store training-time MDS decomposition for stable out-of-sample projection
    eigvecs: np.ndarray
    eigvals: np.ndarray
    row_mean_d2: np.ndarray
    grand_mean_d2: float


def disstree_fit(
    D: np.ndarray,
    y,
    *,
    n_components: int = 5,
    max_depth: Optional[int] = None,
    min_samples_leaf: int = 2,
    random_state: int = 0,
) -> DissTreeResult:
    """Fit a scikit-learn decision tree to a label on a classical MDS embedding (exploratory).

    Not the TraMineR ``disstree`` (no weights, no inertia splits); for the
    R-faithful weighted regression tree use :func:`traminex.disstree.disstree`.
    """
    sk_tree = _sklearn_tree()
    if sk_tree is None:
        raise RuntimeError("scikit-learn is required for disstree_fit")
    X, eigvecs, eigvals, row_mean_d2, grand_mean_d2 = _cmdscale_fit(D, n_components=n_components)
    clf = sk_tree.DecisionTreeClassifier(
        max_depth=max_depth, min_samples_leaf=min_samples_leaf, random_state=random_state
    )
    y = np.asarray(y)
    clf.fit(X, y)
    return DissTreeResult(
        model=clf,
        coords=X,
        labels=y,
        eigvecs=eigvecs,
        eigvals=eigvals,
        row_mean_d2=row_mean_d2,
        grand_mean_d2=grand_mean_d2,
    )


def disstree_plot(result: DissTreeResult, ax=None):
    """Plot an exploratory sklearn discrepancy tree on a Matplotlib axes."""
    if not isinstance(result, DissTreeResult):
        raise TypeError("result must be a DissTreeResult")
    sk_tree = _sklearn_tree()
    if sk_tree is None:
        raise RuntimeError("scikit-learn is required for disstree_plot")
    if ax is None:
        plt = _import_pyplot("disstree_plot")
        fig, ax = plt.subplots(figsize=(10, 6))
    sk_tree.plot_tree(
        result.model, class_names=[str(c) for c in np.unique(result.labels)], filled=True, ax=ax
    )
    return ax


def disstree_leaf_summary(result: DissTreeResult) -> pd.DataFrame:
    """Return per-leaf counts similar to seqtree_leaf_summary."""
    if not isinstance(result, DissTreeResult):
        raise TypeError("result must be a DissTreeResult")
    leaf = result.model.apply(result.coords)
    uniq, cnt = np.unique(leaf, return_counts=True)
    return pd.DataFrame({"leaf": uniq, "size": cnt}).sort_values("leaf").reset_index(drop=True)


def disstree_feature_importance(result: DissTreeResult) -> pd.DataFrame:
    """Return feature importances with component names comp1..compK."""
    if not isinstance(result, DissTreeResult):
        raise TypeError("result must be a DissTreeResult")
    imp = getattr(result.model, "feature_importances_", None)
    if imp is None:
        return pd.DataFrame({"feature": [], "importance": []})
    names = [f"comp{i+1}" for i in range(result.coords.shape[1])]
    return (
        pd.DataFrame({"feature": names, "importance": imp})
        .sort_values("importance", ascending=False)
        .reset_index(drop=True)
    )


def disstree_predict(result: DissTreeResult, D_new: np.ndarray) -> np.ndarray:
    """Predict class labels for new samples given distances to training set.

    Uses classical-MDS out-of-sample projection with training eigen-structure.
    """
    if not isinstance(result, DissTreeResult):
        raise TypeError("result must be a DissTreeResult")
    n_train = result.coords.shape[0]
    D_new = np.asarray(D_new, dtype=float)
    if D_new.shape[1] != n_train:
        raise ValueError("D_new must be (n_new, n_train) distances to training set")
    n_components = result.coords.shape[1]
    n_new = D_new.shape[0]

    # If no positive eigenvalues were available at fit-time, all coordinates are zero.
    if result.eigvals.size == 0 or result.eigvecs.size == 0:
        X_new = np.zeros((n_new, n_components), dtype=float)
        return result.model.predict(X_new)

    # Classical MDS out-of-sample formula (Gower-style centering).
    d2 = D_new**2
    mean_d2_new = d2.mean(axis=1, keepdims=True)
    b = -0.5 * (d2 - result.row_mean_d2[np.newaxis, :] - mean_d2_new + float(result.grand_mean_d2))

    X_proj = b.dot(result.eigvecs) / np.sqrt(result.eigvals)[np.newaxis, :]
    if X_proj.shape[1] < n_components:
        pad = np.zeros((n_new, n_components - X_proj.shape[1]), dtype=float)
        X_new = np.hstack([X_proj, pad])
    else:
        X_new = X_proj[:, :n_components]
    return result.model.predict(X_new)


def disstree2dot_explore(result: DissTreeResult, out_file: str) -> None:
    """Export the *exploratory* (MDS + sklearn) disstree to a Graphviz DOT file.

    This is the sklearn ``export_graphviz`` over the exploratory
    :class:`DissTreeResult`; the R-faithful, tree-structure DOT export is
    :func:`traminex.disstree_viz.disstree2dot` (public ``disstree2dot``).
    Requires sklearn.tree.export_graphviz. out_file is a path to write the DOT.
    """
    if not isinstance(result, DissTreeResult):
        raise TypeError("result must be a DissTreeResult")
    sk_tree = _sklearn_tree()
    if sk_tree is None:
        raise RuntimeError("scikit-learn is required for disstree2dot_explore")
    feat_names = [f"comp{i+1}" for i in range(result.coords.shape[1])]
    class_names = [str(c) for c in np.unique(result.labels)]
    sk_tree.export_graphviz(
        result.model,
        out_file=out_file,
        feature_names=feat_names,
        class_names=class_names,
        filled=True,
    )
