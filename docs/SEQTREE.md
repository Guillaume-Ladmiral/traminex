# Sequence trees

`seqtree` is the R-faithful TraMineR sequence-tree entry point. It computes or
accepts sequence dissimilarities and delegates to the weighted, permutation-based
`disstree` implementation.

The separately named `seqtree_fit` is an exploratory sklearn model over
sequence-derived features. It is useful for experimentation but is not the
TraMineR algorithm.

## R-faithful API

- `seqtree(seqdata, predictors, *, diss=None, ...)` returns a `DissTree`.
  Where R's `seqtree(formula, data, ...)` names the sequence object and the
  covariates in a formula, traminex takes the state-sequence object and a
  predictor table as its first two arguments.
- `disstreeleaf`, `disstree_get_rules`, and `disstree_assign` inspect and apply
  fitted trees.
- `seqtree2dot` and `seqtreedisplay` export or render the tree.

## Exploratory sklearn API

- seqtree_fit(x, y, features="durations+transitions+summary", max_depth=None, min_samples_leaf=2, random_state=0)
  - x: STSList, y: labels (array-like)
  - features: combination of durations, transitions, summary
  - returns SeqTreeResult with model, feature_names, labels, leaf_index
- seqtree_leaf_summary(result)
  - Returns a small DataFrame with leaf sizes
  - The historical `seqtree_leaf_summary(x, result)` call remains temporarily
    available with a deprecation warning; `x` was never used by this helper

## Features

- durations: per-state durations (seqistatd)
- transitions: per-sequence transition counts (flattened K×K)
- summary: per-sequence normalized entropy of the state shares, and the number of
  changes between two adjacent state cells divided by the sequence width minus 1
  (a missing or void cell joins no transition; this is not `seqtransn_count`,
  which counts as R's `seqtransn`)

Do not use the exploratory API as a parity substitute for `seqtree`.

Tree fitting does not require the Graphviz executable. DOT rendering through
`seqtreedisplay`/`disstreedisplay` uses the `tree` extra and a system `dot` on
`PATH`. The maintained isolated tree profile renders a real SVG from an
installed sdist, so import-only Graphviz checks are not sufficient.
