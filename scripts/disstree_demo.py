#!/usr/bin/env python3
"""
Disstree demo CLI.

Usage:
  python scripts/disstree_demo.py --dist D.csv --labels y.csv [--n_components 5] [--max_depth 3] [--min_samples_leaf 2]
                                  [--plot tree.png] [--predict D_new.csv]

- D.csv: square distance matrix (n x n), numeric, CSV. Header optional.
- y.csv: labels vector (n,), CSV with either a headerless single column or a
  column headed 'label'.
- D_new.csv: distances from new points to training set (n_new x n), CSV (optional).

Outputs:
- Prints leaf summary and feature importances to stdout.
- Saves a tree plot if --plot is given.
- If --predict is provided, prints predictions for D_new.
"""

import argparse
import sys

import numpy as np
import pandas as pd

from traminex.disstree import (
    disstree_feature_importance,
    disstree_fit,
    disstree_leaf_summary,
    disstree_plot,
    disstree_predict,
)


def load_matrix(path: str) -> np.ndarray:
    df = pd.read_csv(path, header=None)
    try:
        return df.to_numpy(dtype=float)
    except ValueError:
        if df.shape[0] < 2:
            raise
        # A non-numeric first row is the documented optional column header.
        return df.iloc[1:].to_numpy(dtype=float)


def load_labels(path: str) -> np.ndarray:
    raw = pd.read_csv(path, header=None)
    if raw.shape[1] == 1:
        values = raw.iloc[:, 0]
        if len(values) and str(values.iloc[0]).strip().lower() == "label":
            values = values.iloc[1:]
        y = values.to_numpy()
    else:
        df = pd.read_csv(path)
        y = (df["label"] if "label" in df.columns else df.iloc[:, 0]).to_numpy()
    return np.asarray(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dist", required=True, help="Path to square distance CSV (n x n)")
    ap.add_argument("--labels", required=True, help="Path to labels CSV (n,)")
    ap.add_argument("--n_components", type=int, default=5)
    ap.add_argument("--max_depth", type=int, default=None)
    ap.add_argument("--min_samples_leaf", type=int, default=2)
    ap.add_argument("--random_state", type=int, default=0)
    ap.add_argument("--plot", type=str, default=None, help="Path to save tree plot (PNG)")
    ap.add_argument(
        "--predict",
        type=str,
        default=None,
        help="Path to new distances CSV (n_new x n) for prediction",
    )
    args = ap.parse_args()

    D = load_matrix(args.dist)
    y = load_labels(args.labels)

    if D.shape[0] != D.shape[1]:
        print("Error: --dist must be square (n x n).", file=sys.stderr)
        sys.exit(2)
    if y.shape[0] != D.shape[0]:
        print("Error: labels length must match n in distance matrix.", file=sys.stderr)
        sys.exit(2)

    res = disstree_fit(
        D,
        y,
        n_components=args.n_components,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
        random_state=args.random_state,
    )

    # Print summaries
    print("\n== Leaf summary ==")
    print(disstree_leaf_summary(res).to_string(index=False))

    print("\n== Feature importances ==")
    print(disstree_feature_importance(res).to_string(index=False))

    # Plot if requested
    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ax = disstree_plot(res)
        plt.tight_layout()
        plt.savefig(args.plot, dpi=150)
        plt.close()
        print(f"Saved tree plot to {args.plot}")

    # Predict if requested
    if args.predict:
        D_new = load_matrix(args.predict)
        if D_new.shape[1] != D.shape[0]:
            print(
                "Error: --predict matrix must have n columns equal to training size.",
                file=sys.stderr,
            )
            sys.exit(2)
        pred = disstree_predict(res, D_new)
        print("\n== Predictions for new samples ==")
        for i, p in enumerate(pred):
            print(f"new[{i}] -> {p}")


if __name__ == "__main__":
    main()
