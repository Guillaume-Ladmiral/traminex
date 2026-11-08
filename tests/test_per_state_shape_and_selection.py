"""Per-state result shapes, and the documented STSList row-selection route.

F-003 -- ``STSList`` uses a single integer codebook in which ``x.alphabet`` appends
the ``nr`` and ``void`` markers after the real states, while ``alphabet(x)`` returns
the real states only.  Several per-state results therefore carry extra zero entries
that R never emits.  The existing golden tests filter those markers out before
comparing, so they are correct for values and structurally incapable of asserting
shape -- which is what let the ``seqmeant`` ledger row claim shape parity.  This
module asserts the *length* directly, against R's own row count.

TraMineR 2.2-14's shapes on the fixture below, recorded by
``scripts/r_goldens/make_goldens_audit20260922_s10.R`` in
``tests/r_goldens/inline_r_anchors.txt`` (the values are compared by the
per-function R goldens; only the shapes differ):

    seqmeant(with.missing=FALSE)   3 rows  A B C
    seqmeant(with.missing=TRUE)    4 rows  A B C *
    seqstatd(...)$Frequencies      3 x 6   (states x time)
    seqistatd(...)                 6 x 3
    cpal / stlab                   3
    seqstatf                       3 x 2   <- no padding
    seqtrate                       3 x 3   <- no padding

F-009 -- ``[.stslist`` and ``rbind.stslist`` claimed that "STSList indexing" and
"the public container protocol" replace R's S3 dispatch.  Neither exists:
``STSList`` defines no ``__len__``, ``__getitem__`` or ``__iter__``.  Row selection
is done by slicing ``.data`` and reconstructing, which is what those rows now say
and what the second half of this module exercises end-to-end into ``seqdist``.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from r_goldens_utils import r_anchor_shape

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
REAL_STATES = 3


def _sequences():
    return tx.seqdef(_RAW.copy(), alphabet=["A", "B", "C"])


def _quiet(call, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return call(*args, **kwargs)


# --------------------------------------------------------------------------
# F-003: shapes, asserted against R's row count rather than a filtered subset
# --------------------------------------------------------------------------


def test_the_internal_codebook_is_wider_than_the_real_alphabet():
    x = _sequences()
    assert list(tx.alphabet(x)) == ["A", "B", "C"]
    assert list(x.alphabet) == ["A", "B", "C", "*", "%"]


@pytest.mark.parametrize(
    ("include_missing", "r_case"),
    [(False, "seqmeant_without_missing"), (True, "seqmeant_with_missing")],
    ids=["missing-excluded", "missing-included"],
)
def test_seqmeant_is_marker_padded_relative_to_r(include_missing, r_case):
    """R returns statl (+ nr when with.missing); traminex returns the codebook."""
    r_rows = r_anchor_shape(r_case)[0]
    assert r_rows == REAL_STATES + int(include_missing)
    result = np.asarray(_quiet(tx.seqmeant, _sequences(), include_missing=include_missing))
    assert len(result) == len(_sequences().alphabet) == 5
    assert len(result) != r_rows, "shape parity would make this row EXACT again"
    # The extra entries are the markers, and they are zero when unused.
    assert float(result.ravel()[-1]) == 0.0


@pytest.mark.parametrize(
    ("name", "getter", "expected_axis"),
    [
        ("seqistatd", lambda x: np.asarray(_quiet(tx.seqistatd, x, include_missing=False)), 1),
        (
            "seqstatd",
            lambda x: np.asarray(_quiet(tx.seqstatd, x, include_missing=False)["Frequencies"]),
            1,
        ),
    ],
)
def test_per_state_tables_follow_the_internal_codebook(name, getter, expected_axis):
    values = getter(_sequences())
    assert values.shape[expected_axis] == 5, f"{name}: expected the padded codebook width"


@pytest.mark.parametrize("accessor", [tx.cpal, tx.stlab], ids=["cpal", "stlab"])
def test_metadata_accessors_follow_the_internal_codebook(accessor):
    assert len(accessor(_sequences())) == 5


@pytest.mark.parametrize(
    ("name", "getter"),
    [
        ("seqstatf", lambda x: np.asarray(_quiet(tx.seqstatf, x))),
        ("seqtrate", lambda x: np.asarray(_quiet(tx.seqtrate, x))),
    ],
)
def test_the_padding_is_not_uniform(name, getter):
    """seqstatf and seqtrate follow the real alphabet exactly as R does.

    This is the part that makes the padding a per-function decision rather than an
    unavoidable consequence of the codebook design, so it must be pinned.
    """
    assert getter(_sequences()).shape == r_anchor_shape(name)


# --------------------------------------------------------------------------
# F-009: the documented row-selection route
# --------------------------------------------------------------------------


def test_stslist_offers_no_container_protocol():
    """Pin the absence, so the ledger rationale cannot drift back to claiming one."""
    x = _sequences()
    with pytest.raises(TypeError):
        len(x)
    with pytest.raises(TypeError):
        x[0]
    with pytest.raises(TypeError):
        iter(x)


def test_row_selection_via_data_and_reconstruction_feeds_seqdist():
    """The route the ledger now documents, exercised end to end."""
    from traminex.sequences import STSList, _remap_issues

    x = _sequences()
    rows = np.array([0, 2])
    selected = STSList(
        data=np.asarray(x.data)[rows],
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
        ids=None if x.ids is None else [x.ids[i] for i in rows],
        issues=_remap_issues(x.issues, rows),
    )

    assert np.asarray(selected.data).shape == (2, 6)
    assert selected.nseq == 2
    assert list(tx.alphabet(selected)) == list(tx.alphabet(x))

    distances = np.asarray(_quiet(tx.seqdist, selected, method="LCS", with_missing=True))
    assert distances.shape == (2, 2)
    assert distances[0, 0] == 0.0
    assert distances[0, 1] == distances[1, 0]

    full = np.asarray(_quiet(tx.seqdist, x, method="LCS", with_missing=True))
    np.testing.assert_allclose(distances, full[np.ix_([0, 2], [0, 2])], rtol=0, atol=0)


def test_the_copy_and_slice_route_also_reconstructs_a_usable_object():
    """The lighter route the seqLLCS golden helper uses, made explicit."""
    x = _sequences()
    view = x.copy()
    view.data = np.asarray(x.data)[1:3]

    assert view.nseq == 2
    assert list(view.alphabet) == list(x.alphabet)
    distances = np.asarray(_quiet(tx.seqdist, view, method="LCS", with_missing=True))
    assert distances.shape == (2, 2)


# --------------------------------------------------------------------------
# Package-wide per-SEQUENCE orientation (found while completing RESIDUAL-4)
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name",
    ["seqient", "seqtransn", "seqlength", "seqici", "seqST"],
)
def test_per_sequence_indicators_are_one_dimensional(name):
    """R returns these as an n x 1 matrix; traminex returns a 1-D array of length n.

    TraMineR 2.2-14 on this fixture: each is 6x1 (recorded R shape). The
    convention is documented once in docs/API_CONSISTENCY.md rather than
    repeated in five decision rows.
    """
    assert r_anchor_shape(f"per_sequence_{name}") == (6, 1)
    values = np.asarray(_quiet(getattr(tx, name), _sequences()))
    assert values.ndim == 1
    assert values.shape == (6,)


@pytest.mark.parametrize("name", ["seqivolatility", "seqivardur"])
def test_the_indicators_r_already_returns_flat_are_unchanged(name):
    """R returns a plain vector for these two, so there is no orientation gap."""
    assert r_anchor_shape(f"per_sequence_{name}") == (6,)
    values = np.asarray(_quiet(getattr(tx, name), _sequences()))
    assert values.ndim == 1
    assert values.shape == (6,)
