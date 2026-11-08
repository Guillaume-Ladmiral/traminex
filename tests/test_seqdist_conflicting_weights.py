"""Explicit sequence weights next to a weighted object are refused, not dropped.

``seqdist`` takes the sequence weights of CHI2/EUCLID-style and OMstran
distances (``weights=``) and of the TRATE cost fit (``trate_weights=``) from the
object when it carries them, and from the keyword otherwise. For a weighted
object the keyword was silently discarded, while the same keyword on a plain
array or an unweighted object was honoured (2026-09-22 audit S1a-06). R has no
such keyword -- weights live on the object and ``weighted=`` selects them -- so
the combination is now refused by name.
"""

from __future__ import annotations

import numpy as np
import pytest

import traminex as tx

X = np.array([[0, 1, 2, 1], [1, 1, 0, 2], [2, 0, 1, 1], [0, 2, 2, 0]])
CALLER = np.array([1.0, 3.0, 0.5, 2.0])


def _object(weights=None):
    return tx.seqdef(X, alphabet=[0, 1, 2], weights=weights)


@pytest.mark.parametrize(
    "method,extra,keyword",
    [
        ("CHI2", {}, "weights"),
        ("OMSTRAN", {"subm": "CONSTANT", "otto": 0.5}, "weights"),
        ("OM", {"subm": "TRATE"}, "trate_weights"),
    ],
)
def test_a_weighted_object_refuses_an_explicit_weight_keyword(method, extra, keyword):
    weighted = _object(weights=np.array([2.0, 1.0, 1.0, 1.0]))
    with pytest.raises(ValueError, match=keyword):
        tx.seqdist(weighted, method=method, **extra, **{keyword: CALLER})


@pytest.mark.parametrize(
    "method,extra,keyword",
    [("CHI2", {}, "weights"), ("OM", {"subm": "TRATE"}, "trate_weights")],
)
def test_an_unweighted_object_still_uses_the_keyword_like_a_plain_array(method, extra, keyword):
    from_object = tx.seqdist(_object(), method=method, **extra, **{keyword: CALLER})
    from_array = tx.seqdist(X, method=method, **extra, **{keyword: CALLER})
    np.testing.assert_array_equal(from_object, from_array)
