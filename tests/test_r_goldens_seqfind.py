"""Public TraMineR 2.2-14 evidence for reference-major ``seqfind``."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from traminex import seqconc, seqdef, seqfind

GOLDEN_DIR = Path(__file__).resolve().parent / "r_goldens"


def _fixtures():
    reference = seqdef(
        np.array(
            [
                ["A", "B", "%", "%"],
                ["A", "B", "%", "%"],
                ["C", "*", "A", "%"],
                ["B", "A", "%", "%"],
                ["A", "C", "%", "%"],
            ],
            dtype=object,
        ),
        alphabet=["A", "B", "C"],
        nr="*",
    )
    target = seqdef(
        np.array(
            [
                ["A", "B", "%"],
                ["C", "*", "A"],
                ["A", "B", "%"],
                ["B", "A", "%"],
                ["C", "C", "%"],
            ],
            dtype=object,
        ),
        alphabet=["C", "B", "A"],
        nr="*",
    )
    return reference, target


def _real_alphabet(x) -> str:
    markers = {code for code in (x.nr_code, x.void_code) if code is not None}
    return ",".join(str(label) for code, label in enumerate(x.alphabet) if code not in markers)


def test_seqfind_fixture_matches_recorded_public_r_inputs():
    reference, target = _fixtures()
    golden = pd.read_csv(
        GOLDEN_DIR / "seqfind_public_inputs.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )

    for name, x in (("reference", reference), ("target", target)):
        rows = golden.loc[golden["object"] == name]
        assert rows["width"].astype(int).unique().tolist() == [x.seqlen]
        assert rows["conc"].tolist() == seqconc(x).tolist()
        assert rows["alphabet"].unique().tolist() == [_real_alphabet(x)]
        assert rows["void"].unique().tolist() == [x.void]

    assert reference.seqlen != target.seqlen
    assert _real_alphabet(reference) == "A,B,C"
    assert _real_alphabet(target) == "C,B,A"


def test_seqfind_matches_reference_major_public_r_indexes():
    reference, target = _fixtures()
    golden = pd.read_csv(
        GOLDEN_DIR / "seqfind_public_results.txt",
        sep="|",
        dtype=str,
        keep_default_na=False,
    )
    r_indexes = golden.loc[golden["case"] == "reference_major", "r_indexes"].item()
    expected = np.fromstring(r_indexes, dtype=int, sep=",") - 1

    actual = seqfind(reference, target)

    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(actual, np.array([0, 2, 0, 2, 1, 3]))


def test_seqfind_no_match_is_empty_integer_array():
    _, target = _fixtures()
    no_match = seqdef(
        np.array([["A", "C", "%", "%"]], dtype=object),
        alphabet=["A", "B", "C"],
        nr="*",
    )

    actual = seqfind(no_match, target)

    assert actual.dtype.kind == "i"
    assert actual.shape == (0,)


def test_seqfind_validates_both_sequence_objects_without_mutation():
    reference, target = _fixtures()
    reference_before = reference.data.copy()
    target_before = target.data.copy()
    with pytest.raises(TypeError, match="ref|reference|STSList"):
        seqfind(np.asarray(reference.data), target)
    with pytest.raises(TypeError, match="target|STSList"):
        seqfind(reference, np.asarray(target.data))

    seqfind(reference, target)
    np.testing.assert_array_equal(reference.data, reference_before)
    np.testing.assert_array_equal(target.data, target_before)
