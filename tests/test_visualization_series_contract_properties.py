"""Property tests for ``docs/VISUALIZATION_SERIES_CONTRACT.md`` (RESIDUAL-5).

That contract is consumed by the sibling ``py_ggplotseq`` project, which reshapes
these outputs but must not recompute the numerics.  The cycle-4b audit reviewed it
at source and recorded that its non-mutation and metadata-completeness claims were
never independently property-tested.  This module tests them.

Two claims are exercised:

  * "Public analytical functions and series helpers do not mutate the caller"
    (line 41) -- every documented series helper is called and the whole input
    object is compared field-by-field against a pre-call deep snapshot.
  * The ``STSList`` field contract (lines 29-40) -- ``weights`` and ``ids``, when
    present, have length ``N``; codes are zero-based indexes into ``alphabet``;
    ``labels``/``cpal``, when present, follow the same state order; and derived
    one-row objects keep the applicable metadata.

The helper set is the one the contract names at line 19: ``seqdss``, ``seqtrate``,
``seqstatd``, ``seqmeant``, ``seqmodst``, ``seqtab``, ``seqrep``, ``seqrf``.
"""

from __future__ import annotations

import copy
import warnings

import numpy as np
import pytest

import traminex as tx

_RAW = np.array(
    [
        [None, "A", "B", "C", "A", None],
        ["A", "B", None, "C", "C", "B"],
        ["B", "B", "A", None, None, "C"],
        ["C", "A", "A", "B", "B", None],
        [None, None, "C", "C", "A", "A"],
        ["A", "C", "B", "B", "A", "C"],
    ],
    dtype=object,
)

_MUTABLE_FIELDS = (
    "data",
    "alphabet",
    "labels",
    "cpal",
    "xtlab",
    "weights",
    "ids",
    "issues",
)


def _sequences(weighted: bool = False):
    x = tx.seqdef(_RAW.copy(), alphabet=["A", "B", "C"])
    if weighted:
        x.weights = np.array([1.0, 2.0, 0.5, 1.5, 1.0, 3.0])
    x.ids = [f"p{i}" for i in range(x.nseq)]
    return x


def _snapshot(x) -> dict:
    return {name: copy.deepcopy(getattr(x, name, None)) for name in _MUTABLE_FIELDS}


def _assert_unchanged(before: dict, x, label: str) -> None:
    for name, original in before.items():
        current = getattr(x, name, None)
        if isinstance(original, np.ndarray):
            assert np.array_equal(original, np.asarray(current)), f"{label} mutated x.{name}"
        else:
            assert original == current, f"{label} mutated x.{name}"


def _quiet(call, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return call(*args, **kwargs)


def _distances(x):
    return np.asarray(_quiet(tx.seqdist, x, method="LCS", with_missing=True), dtype=float)


SERIES_HELPERS = [
    ("seqdss", lambda x: _quiet(tx.seqdss, x)),
    ("seqtrate", lambda x: _quiet(tx.seqtrate, x)),
    ("seqstatd", lambda x: _quiet(tx.seqstatd, x)),
    ("seqmeant", lambda x: _quiet(tx.seqmeant, x)),
    ("seqmodst", lambda x: _quiet(tx.seqmodst, x)),
    ("seqtab", lambda x: _quiet(tx.seqtab, x)),
    ("seqrep", lambda x: _quiet(tx.seqrep, x, diss=_distances(x))),
    ("seqrf", lambda x: _quiet(tx.seqrf, x, _distances(x), k=2)),
]


@pytest.mark.parametrize(("name", "call"), SERIES_HELPERS, ids=[n for n, _ in SERIES_HELPERS])
@pytest.mark.parametrize("weighted", [False, True], ids=["unweighted", "weighted"])
def test_series_helpers_do_not_mutate_the_caller(name, call, weighted):
    x = _sequences(weighted=weighted)
    before = _snapshot(x)

    call(x)

    _assert_unchanged(before, x, name)


@pytest.mark.parametrize(("name", "call"), SERIES_HELPERS, ids=[n for n, _ in SERIES_HELPERS])
def test_series_helpers_do_not_mutate_the_caller_when_called_twice(name, call):
    """Repeat invocation must be free of accumulated state as well."""
    x = _sequences(weighted=True)
    call(x)
    before = _snapshot(x)

    call(x)

    _assert_unchanged(before, x, f"{name} (second call)")


def test_stslist_field_contract_holds():
    x = _sequences(weighted=True)

    assert np.asarray(x.data).dtype.kind in {"i", "u"}
    assert np.asarray(x.data).shape == (x.nseq, x.seqlen)
    assert len(x.weights) == x.nseq
    assert len(x.ids) == x.nseq
    if x.labels is not None:
        assert len(x.labels) == len(x.alphabet)
    if x.cpal is not None:
        assert len(x.cpal) == len(x.alphabet)
    codes = np.asarray(x.data)
    assert int(codes.min()) >= 0
    assert int(codes.max()) < len(x.alphabet), "codes must index alphabet"


def test_seqmodst_returns_a_metadata_complete_one_row_object():
    """The contract calls this "a one-row STSList of length L" (line 81)."""
    x = _sequences()
    modal = _quiet(tx.seqmodst, x)

    assert modal.nseq == 1
    assert modal.seqlen == x.seqlen
    assert list(modal.alphabet) == list(x.alphabet)
    assert modal.nr_code == x.nr_code
    assert modal.void_code == x.void_code
    assert modal.start == x.start
    assert modal.xtstep == x.xtstep


def test_derived_objects_do_not_alias_the_source_metadata():
    """Line 100: derived arrays and mutable metadata are copies."""
    x = _sequences()
    modal = _quiet(tx.seqmodst, x)

    modal.alphabet.append("ZZZ")
    np.asarray(modal.data)[...] = 0

    assert "ZZZ" not in x.alphabet
    assert int(np.asarray(x.data).sum()) > 0


def test_malformed_input_is_refused_rather_than_repaired():
    """Line 43-45: malformed IDs and weights must raise, not be silently fixed."""
    x = _sequences()

    bad_weights = x.copy()
    with pytest.raises(ValueError, match="weights must have length"):
        bad_weights.weights = np.array([1.0, 2.0])

    with pytest.raises(TypeError):
        _quiet(tx.seqstatd, np.asarray(x.data))
