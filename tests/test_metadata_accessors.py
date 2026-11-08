"""Focused public evidence for TraMineR metadata helper mappings."""

from __future__ import annotations

import numpy as np

from traminex import STSList, alphabet, cpal, seqdef, seqhasmiss, stlab


def test_cpal_and_stlab_return_independent_metadata_values() -> None:
    x = seqdef(
        [["A", "B"], ["B", "A"]],
        alphabet=["A", "B"],
        labels=["State A", "State B"],
        cpal=["#112233", "#445566"],
    )

    palette = cpal(x)
    labels = stlab(x)
    assert palette == ["#112233", "#445566"]
    assert labels == ["State A", "State B"]

    assert palette is not None
    palette.append("#778899")
    labels.append("State C")
    assert x.cpal == ["#112233", "#445566"]
    assert x.labels == ["State A", "State B"]


def test_python_metadata_mutation_reexpresses_r_setter_routes() -> None:
    x = seqdef(
        [["A", "B"]],
        alphabet=["A", "B"],
        labels=["State A", "State B"],
        cpal=["red", "blue"],
    )

    x.alphabet = ["X", "Y"]
    assert alphabet(x) == ["X", "Y"]
    # Unlike TraMineR's alphabet<- route, direct dataclass mutation does not
    # synchronize labels that previously equalled the alphabet.
    assert stlab(x) == ["State A", "State B"]

    x.cpal = ["black", "white"]
    x.labels = ["First", "Second"]
    assert cpal(x) == ["black", "white"]
    assert stlab(x) == ["First", "Second"]


def test_seqhasmiss_reexpresses_combined_missing_presence_without_mutation() -> None:
    x = STSList(
        data=np.array([[0, 2, 1], [1, 0, 3], [0, 1, 0]], dtype=int),
        alphabet=["A", "B", "*", "%"],
        labels=["A", "B", "Missing", "Void"],
        nr_code=2,
        void_code=3,
    )
    original = x.data.copy()

    assert seqhasmiss(x) is True
    per_sequence = seqhasmiss(x, per_sequence=True)
    assert isinstance(per_sequence, np.ndarray)
    assert per_sequence.dtype == np.bool_
    assert per_sequence.tolist() == [True, True, False]
    np.testing.assert_array_equal(x.data, original)


def test_seqhasmiss_without_declared_markers_is_false() -> None:
    x = seqdef([["A", "B"], ["B", "A"]])

    assert seqhasmiss(x) is False
    np.testing.assert_array_equal(
        seqhasmiss(x, per_sequence=True),
        np.array([False, False]),
    )
