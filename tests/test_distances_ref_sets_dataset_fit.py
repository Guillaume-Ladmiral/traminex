"""Cross-block distances must not re-estimate dataset-fitted parameters (F-001).

TraMineR fits substitution costs and state frequencies from the whole ``seqdata``
before selecting a reference set, so ``seqdist(s, sm="TRATE", refseq=list(1:2, 3:4))``
equals the corresponding block of the full-corpus matrix exactly.  The
``ref_sets`` route of :func:`traminex.distances.pairwise` must reproduce that.

The expected value here is derived from a corpus **strictly larger** than
``A | B``.  Deriving it from ``vstack([A, B])`` — as the previous cross-block test
did — cannot detect a re-fit by construction, because that is the very corpus a
re-fitting implementation estimates from.
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex.distances as distances

# Twenty sequences over a four-state alphabet.  Rows 0-3 (the reference sets) use a
# deliberately narrow slice of the alphabet, so a fit restricted to them differs
# sharply from a fit over all twenty rows.
_CORPUS = np.vstack(
    [
        np.array(
            [
                [0, 0, 1, 1, 0, 0, 1, 1],
                [0, 1, 1, 0, 0, 1, 1, 0],
                [1, 1, 0, 0, 1, 1, 0, 0],
                [0, 0, 1, 0, 0, 0, 1, 0],
            ],
            dtype=np.int32,
        ),
        np.random.default_rng(20260823).integers(0, 4, size=(16, 8)).astype(np.int32),
    ]
)
_SET_A = np.array([0, 1])
_SET_B = np.array([2, 3])

# The five configurations whose parameters are estimated from the corpus, with
# ``subm`` left at its default so the dataset-fitted route is the one exercised.
DATASET_FITTED = [
    ("OM", {"indel": 1.0, "subm": "TRATE"}),
    ("DHD", {}),
    ("CHI2", {}),
    ("EUCLID", {}),
    ("OMSTRAN", {"otto": 0.5}),
]

# Configurations whose distances depend only on the pair being compared.  These
# must keep the efficient rectangular path and stay byte-identical.
PAIR_LOCAL = [
    ("OM", {"indel": 1.0, "subm": 2.0}),
    ("DHD", {"subm": np.ones((4, 4)) - np.eye(4)}),
    ("LCS", {}),
    ("LCP", {}),
    ("RLCP", {}),
    ("HAM", {}),
    ("OMSPELL", {"indel": 1.0, "subm": 2.0, "expcost": 0.5, "tpow": 1.0}),
]


@pytest.mark.parametrize(("method", "options"), DATASET_FITTED, ids=[m for m, _ in DATASET_FITTED])
def test_ref_sets_cross_block_is_fitted_on_the_whole_corpus(method, options):
    full = distances.pairwise(_CORPUS.copy(), method=method, **options)
    expected = np.asarray(full)[np.ix_(_SET_A, _SET_B)]

    actual = distances.pairwise(_CORPUS.copy(), method=method, ref_sets=(_SET_A, _SET_B), **options)

    assert actual.shape == (len(_SET_A), len(_SET_B))
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize(("method", "options"), PAIR_LOCAL, ids=[m for m, _ in PAIR_LOCAL])
def test_ref_sets_cross_block_is_unchanged_for_pair_local_methods(method, options):
    full = distances.pairwise(_CORPUS.copy(), method=method, **options)
    expected = np.asarray(full)[np.ix_(_SET_A, _SET_B)]

    actual = distances.pairwise(_CORPUS.copy(), method=method, ref_sets=(_SET_A, _SET_B), **options)

    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize(("method", "options"), DATASET_FITTED, ids=[m for m, _ in DATASET_FITTED])
def test_pairwise_between_sets_refuses_a_dataset_fitted_configuration(method, options):
    """Given only A and B, the corpus to fit from is unknown, so refuse rather than guess."""
    with pytest.raises(ValueError, match="estimated from the whole corpus"):
        distances.pairwise_between_sets(
            _CORPUS[_SET_A].copy(), _CORPUS[_SET_B].copy(), method=method, **options
        )


@pytest.mark.parametrize(("method", "options"), PAIR_LOCAL, ids=[m for m, _ in PAIR_LOCAL])
def test_pairwise_between_sets_still_serves_pair_local_configurations(method, options):
    """The five shapes `weightedcluster._seqdist` delegates must keep working."""
    full = distances.pairwise(_CORPUS.copy(), method=method, **options)
    expected = np.asarray(full)[np.ix_(_SET_A, _SET_B)]

    actual = distances.pairwise_between_sets(
        _CORPUS[_SET_A].copy(), _CORPUS[_SET_B].copy(), method=method, **options
    )

    np.testing.assert_array_equal(actual, expected)


def test_ref_sets_refusal_names_the_supported_routes():
    with pytest.raises(ValueError) as excinfo:
        distances.pairwise_between_sets(
            _CORPUS[_SET_A].copy(), _CORPUS[_SET_B].copy(), method="CHI2"
        )
    message = str(excinfo.value)
    assert "pairwise(" in message
    assert "ref_sets" in message
    assert "seqdist" in message
