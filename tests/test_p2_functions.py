"""Tests for Priority 2 parity functions."""

import numpy as np
from r_goldens_utils import r_anchor_values

from traminex import (
    dissrf,
    seqalign,
    seqdef,
    seqdist,
    seqgen,
    seqici,
    seqient,
    seqrecode,
    seqsubm,
)


def test_seqient_basic():
    """Test within-sequence entropy."""
    data = np.array(
        [
            ["A", "A", "A", "A"],  # Low entropy (all same)
            ["A", "B", "C", "D"],  # High entropy (all different)
            ["A", "A", "B", "B"],  # Medium entropy
        ]
    )
    sts = seqdef(data)
    ent = seqient(sts, norm=True)

    # Low < Medium < High
    assert ent[0] < ent[2] < ent[1]
    # All same should have entropy ~0
    assert ent[0] < 0.01


def test_seqient_unnormalized():
    """Test within-sequence entropy without normalization."""
    data = np.array(
        [
            ["A", "A", "A", "A"],
            ["A", "B", "A", "B"],
        ]
    )
    sts = seqdef(data)
    ent = seqient(sts, norm=False)

    assert ent[0] < ent[1]


def test_seqici_basic():
    """Test complexity index."""
    data = np.array(
        [
            ["A", "A", "A", "A"],  # Low complexity (no diversity, no change)
            ["A", "B", "C", "D"],  # High complexity (high diversity, many changes)
            ["A", "A", "A", "B"],  # Low complexity (low diversity, one change)
        ]
    )
    sts = seqdef(data)
    complexity = seqici(sts)

    # Low < Low < High
    assert complexity[0] < complexity[2] < complexity[1]
    # Sequence with no changes should have complexity 0
    assert complexity[0] == 0.0


def test_seqgen_basic():
    """Test random sequence generation."""
    data = seqgen(10, 5, ["A", "B", "C"], seed=42)

    assert data.shape == (10, 5)
    # All states should be from alphabet
    unique_states = np.unique(data)
    assert all(s in ["A", "B", "C"] for s in unique_states)


def test_seqgen_with_probabilities():
    """Test random generation with probabilities."""
    # Generate many sequences to test probability
    data = seqgen(1000, 10, ["A", "B", "C"], p=[0.7, 0.2, 0.1], seed=42)

    # Count occurrences
    unique, counts = np.unique(data, return_counts=True)
    props = dict(zip(unique, counts / counts.sum()))

    # Check rough proportions (with tolerance)
    assert abs(props.get("A", 0) - 0.7) < 0.05
    assert abs(props.get("B", 0) - 0.2) < 0.05
    assert abs(props.get("C", 0) - 0.1) < 0.05


def test_seqgen_reproducibility():
    """Test that seed produces reproducible results."""
    data1 = seqgen(10, 5, ["A", "B"], seed=42)
    data2 = seqgen(10, 5, ["A", "B"], seed=42)

    assert np.array_equal(data1, data2)


def test_seqrecode_basic():
    """Test alphabet recoding."""
    data = np.array(
        [
            ["FT", "FT", "PT", "UE"],
            ["PT", "PT", "FT", "FT"],
        ]
    )
    sts = seqdef(data)

    # Recode: FT/PT -> Work, UE -> NotWork
    recoded = seqrecode(sts, {"Work": ["FT", "PT"], "NotWork": ["UE"]})

    assert "Work" in recoded.alphabet
    assert "NotWork" in recoded.alphabet
    assert recoded.nseq == 2
    assert recoded.seqlen == 4


def test_seqrecode_with_otherwise():
    """Test recoding with 'otherwise' parameter."""
    data = np.array(
        [
            ["A", "B", "C", "D"],
            ["A", "A", "B", "B"],
        ]
    )
    sts = seqdef(data)

    # Recode A/B -> Group1, everything else -> Other
    recoded = seqrecode(
        sts,
        {
            "Group1": ["A", "B"],
        },
        otherwise="Other",
    )

    assert "Group1" in recoded.alphabet
    assert "Other" in recoded.alphabet
    # C and D should be recoded to 'Other'


def test_seqrecode_preserves_unmapped():
    """Test that unmapped states are preserved when otherwise=None."""
    data = np.array(
        [
            ["A", "B", "C"],
            ["A", "A", "B"],
        ]
    )
    sts = seqdef(data)

    # Recode only A/B
    recoded = seqrecode(
        sts,
        {
            "AB": ["A", "B"],
        },
        otherwise=None,
    )

    assert "AB" in recoded.alphabet
    assert "C" in recoded.alphabet


def test_seqalign_basic():
    """Test pairwise alignment."""
    data = np.array(
        [
            ["A", "B", "C"],
            ["A", "C", "C"],
        ]
    )
    sts = seqdef(data)
    sm = seqsubm(sts.data, method="CONSTANT", cval=2.0)

    alignment = seqalign(sts, (0, 1), indel=1.0, sm=sm)

    # Check basic properties
    assert alignment.total_cost > 0
    assert len(alignment.operation) > 0
    assert len(alignment.seq1) == len(alignment.seq2)
    assert len(alignment.seq1) == len(alignment.operation)
    assert len(alignment.seq1) == len(alignment.cost)


def test_seqalign_identical():
    """Test alignment of identical sequences."""
    data = np.array(
        [
            ["A", "B", "C"],
            ["A", "B", "C"],
        ]
    )
    sts = seqdef(data)
    sm = seqsubm(sts.data, method="CONSTANT", cval=2.0)

    alignment = seqalign(sts, (0, 1), indel=1.0, sm=sm)

    # Identical sequences should have cost 0
    assert alignment.total_cost == 0.0
    # All operations should be 'E' (equal)
    assert all(op == "E" for op in alignment.operation)


def test_seqalign_repr():
    """Test that alignment repr works."""
    data = np.array(
        [
            ["A", "B"],
            ["A", "A"],
        ]
    )
    sts = seqdef(data)
    sm = seqsubm(sts.data, method="CONSTANT", cval=2.0)

    alignment = seqalign(sts, (0, 1), indel=1.0, sm=sm)

    # Should not raise
    repr_str = repr(alignment)
    assert "Sequence Alignment" in repr_str
    assert "OM distance" in repr_str


def test_dissrf_basic():
    """RP4 T13b: faithful dissrf (grp_meth='prop') structure."""
    data = np.array(
        [
            ["A", "A", "B"],
            ["A", "B", "B"],
            ["A", "A", "B"],
            ["B", "B", "C"],
            ["B", "C", "C"],
        ]
    )
    sts = seqdef(data)
    D = seqdist(sts, method="OM", subm=2.0)

    result = dissrf(D, k=2)

    # Faithful return structure (R dissrf, grp.meth='prop').
    for key in ("medoids", "R2", "Fstat", "pvalue", "heights", "sizes"):
        assert key in result
    assert len(result["medoids"]) == 2
    assert 0 <= result["R2"] <= 1
    assert result["Fstat"] >= 0
    # prop heights are all 1/k.
    np.testing.assert_allclose(result["heights"], [0.5, 0.5])
    assert result["sizes"]["k"] == 2
    assert result["sizes"]["ncase"] == 5


def test_dissrf_auto_k():
    """Test dissrf with automatic k selection (k = min(floor(wsum/10), 100))."""
    data = seqgen(50, 5, ["A", "B", "C"], seed=42)
    sts = seqdef(data)
    D = seqdist(sts, method="HAM")

    result = dissrf(D)  # k=None -> floor(50/10)=5
    assert result["sizes"]["k"] == 5


def test_dissrf_with_weights():
    """Test dissrf with sequence weights (prop method)."""
    data = np.array(
        [
            ["A", "A"],
            ["B", "B"],
            ["C", "C"],
            ["A", "B"],
        ]
    )
    sts = seqdef(data, weights=[1.0, 2.0, 1.0, 1.5])
    D = seqdist(sts, method="HAM")

    result = dissrf(D, k=2, weights=sts.weights)
    # TraMineR 2.2-14 normalizes supplied dissrf weights to n cases (R's value
    # recorded by scripts/r_goldens/make_goldens_audit20260922_s10.R).
    assert result["sizes"]["wsum"] == r_anchor_values("dissrf_weighted_wsum")[0]
    assert len(result["medoids"]) == 2


def test_dissrf_first_grouping_is_supported():
    """R2 first-group remainder allocation remains available in the broad suite."""
    data = np.array([["A", "A"], ["B", "B"], ["C", "C"]])
    sts = seqdef(data)
    D = seqdist(sts, method="HAM")
    result = dissrf(D, k=2, grp_meth="first")
    np.testing.assert_array_equal(result["mdsk"], np.array([0, 0, 1]))
    assert result["grp_meth"] == "first"
