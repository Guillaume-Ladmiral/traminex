"""Tests for Priority 1 parity functions."""

import numpy as np

from traminex import (
    dissrep,
    seqdef,
    seqdist,
    seqfind,
    seqfpos,
    seqfposend,
    seqmeant,
    seqmodst,
    seqstatf,
    seqstatl,
    seqsubm,
)


def test_seqmodst_basic():
    """Test modal state sequence construction."""
    data = np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
            ["A", "A", "B", "B"],
        ]
    )
    sts = seqdef(data)
    modal = seqmodst(sts)

    # Modal sequence should be A-A-B-C (most frequent at each position)
    assert modal.to_strings()[0] == "A-A-B-C"
    assert modal.nseq == 1
    assert hasattr(modal, "frequencies_")
    assert hasattr(modal, "occurrences_")
    # One sequence exactly matches: the first
    assert modal.occurrences_ == 1


def test_seqmodst_weighted():
    """Test modal state with weights."""
    data = np.array(
        [
            ["A", "B"],
            ["A", "B"],
            ["B", "A"],
        ]
    )
    sts = seqdef(data, weights=[1.0, 1.0, 2.0])
    modal = seqmodst(sts, weighted=True)

    # With weights [1, 1, 2], position 0: A=2, B=2 (tie, first wins = A)
    # Position 1: B=2, A=2 (tie, first wins = B by alphabet order)
    # Actually A appears at pos 0 in rows 0,1 (weight 2 total), B in row 2 (weight 2)
    # So it's a tie - first by index wins
    assert modal.nseq == 1


def test_seqmeant_basic():
    """Test mean durations."""
    data = np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
        ]
    )
    sts = seqdef(data)
    result = seqmeant(sts)

    # State A: appears 3 times total (2+1), mean = 1.5
    # State B: appears 3 times total (1+2), mean = 1.5
    # State C: appears 2 times total (1+1), mean = 1.0
    assert result.loc["A", "Mean"] == 1.5
    assert result.loc["B", "Mean"] == 1.5
    assert result.loc["C", "Mean"] == 1.0


def test_seqmeant_with_serr():
    """Test mean durations with standard errors."""
    data = np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
        ]
    )
    sts = seqdef(data)
    result = seqmeant(sts, serr=True)

    assert "Var" in result.columns
    assert "Stdev" in result.columns
    assert "SE" in result.columns
    assert len(result) == 3  # A, B, C


def test_seqstatf_basic():
    """Overall (total) state frequency, R-faithful (not first-position)."""
    data = np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
            ["B", "A", "B", "B"],
        ]
    )
    sts = seqdef(data)
    result = seqstatf(sts)

    # Total time in each state across all positions: A=4, B=6, C=2 (sum=12).
    assert result.loc["A", "Freq"] == 4.0
    assert result.loc["B", "Freq"] == 6.0
    assert result.loc["C", "Freq"] == 2.0
    assert abs(result.loc["A", "Percent"] - 100.0 * 4 / 12) < 1e-9


def test_seqstatl_basic():
    """seqstatl lists the distinct states present (sorted), R-faithful."""
    data = np.array(
        [
            ["A", "A", "B", "C"],
            ["A", "B", "B", "C"],
            ["B", "A", "B", "B"],
        ]
    )
    sts = seqdef(data)
    assert seqstatl(sts) == ["A", "B", "C"]


def test_seqstatf_include_missing_flag():
    data = np.array(
        [
            ["*", "A", "B", "*"],
            ["A", "B", "*", "B"],
            ["*", "A", "B", "*"],
        ],
        dtype=object,
    )
    sts = seqdef(data, alphabet=["A", "B", "*"], nr="*")

    out_ex = seqstatf(sts, include_missing=False)
    out_in = seqstatf(sts, include_missing=True)

    # Overall counts: A=3, B=4 (real states); the 5 missing cells -> nr only when
    # include_missing. statl excludes nr/void unless with.missing.
    assert "*" not in out_ex.index
    assert out_ex.loc["A", "Freq"] == 3.0
    assert out_ex.loc["B", "Freq"] == 4.0
    assert np.isclose(float(out_ex["Percent"].sum()), 100.0)

    assert out_in.loc["*", "Freq"] == 5.0
    assert out_in.loc["A", "Freq"] == 3.0
    assert out_in.loc["B", "Freq"] == 4.0
    assert np.isclose(float(out_in["Percent"].sum()), 100.0)


def test_seqstatl_includes_missing_symbol():
    """seqstatl includes the missing symbol when present (R lists every symbol)."""
    data = np.array(
        [
            ["*", "A", "B", "*"],
            ["A", "B", "*", "B"],
            ["*", "A", "B", "*"],
        ],
        dtype=object,
    )
    sts = seqdef(data, alphabet=["A", "B", "*"], nr="*")
    # "*" sorts before the letters (ASCII 0x2A < 0x41).
    assert seqstatl(sts) == ["*", "A", "B"]


def test_seqfind_basic():
    """Test finding sequence occurrences."""
    data1 = np.array([["A", "B", "C"]])
    data2 = np.array(
        [
            ["A", "B", "C"],
            ["A", "B", "B"],
            ["A", "B", "C"],
            ["B", "B", "C"],
        ]
    )
    ref = seqdef(data1)
    target = seqdef(data2)
    result = seqfind(ref, target)

    # Sequences 0 and 2 match the reference
    assert len(result) == 2
    assert 0 in result
    assert 2 in result


def test_seqfind_no_match():
    """Test finding with no matches."""
    data1 = np.array([["A", "A", "A"]])
    data2 = np.array(
        [
            ["A", "B", "C"],
            ["B", "B", "C"],
        ]
    )
    ref = seqdef(data1)
    target = seqdef(data2)
    result = seqfind(ref, target)

    assert len(result) == 0


def test_seqfpos_basic():
    """Test first occurrence positions."""
    data = np.array(
        [
            ["A", "B", "C", "B"],
            ["B", "A", "C", "C"],
            ["C", "C", "A", "A"],
        ]
    )
    sts = seqdef(data)
    result = seqfpos(sts)

    # Sequence 0: A at pos 1, B at pos 2, C at pos 3
    assert result.loc[0, "A"] == 1
    assert result.loc[0, "B"] == 2
    assert result.loc[0, "C"] == 3

    # Sequence 1: B at pos 1, A at pos 2, C at pos 3
    assert result.loc[1, "B"] == 1
    assert result.loc[1, "A"] == 2
    assert result.loc[1, "C"] == 3

    # Sequence 2: C at pos 1, A at pos 3, B never occurs (0)
    assert result.loc[2, "C"] == 1
    assert result.loc[2, "A"] == 3
    assert result.loc[2, "B"] == 0


def test_seqfpos_include_missing_flag():
    data = np.array(
        [
            ["A", "*", "B", "*"],
            ["*", "A", "*", "B"],
        ],
        dtype=object,
    )
    sts = seqdef(data, alphabet=["A", "B", "*"], nr="*")

    out_ex = seqfpos(sts, include_missing=False)
    out_in = seqfpos(sts, include_missing=True)

    assert "*" not in out_ex.columns
    assert "*" in out_in.columns


def test_seqfposend_basic():
    """Test last occurrence positions."""
    data = np.array(
        [
            ["A", "B", "C", "B"],
            ["B", "A", "C", "C"],
            ["C", "C", "A", "A"],
        ]
    )
    sts = seqdef(data)
    result = seqfposend(sts)

    # Sequence 0: A last at pos 1, B last at pos 4, C last at pos 3
    assert result.loc[0, "A"] == 1
    assert result.loc[0, "B"] == 4
    assert result.loc[0, "C"] == 3

    # Sequence 2: C last at pos 2, A last at pos 4, B never occurs (0)
    assert result.loc[2, "C"] == 2
    assert result.loc[2, "A"] == 4
    assert result.loc[2, "B"] == 0


def test_seqfposend_include_missing_flag():
    data = np.array(
        [
            ["A", "*", "B", "*"],
            ["*", "A", "*", "B"],
        ],
        dtype=object,
    )
    sts = seqdef(data, alphabet=["A", "B", "*"], nr="*")

    out_ex = seqfposend(sts, include_missing=False)
    out_in = seqfposend(sts, include_missing=True)

    assert "*" not in out_ex.columns
    assert "*" in out_in.columns


def test_dissrep_density():
    """Test representative selection with density criterion."""
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
    sm = seqsubm(sts.data, method="CONSTANT", cval=2.0)
    D = seqdist(sts, method="OM", indel=1.0, subm=sm)
    result = dissrep(D, criterion="density", coverage=0.6, pradius=0.1)

    assert "representatives" in result
    assert "statistics" in result
    assert "quality" in result
    assert len(result["representatives"]) >= 1
    assert isinstance(result["statistics"], type(seqmeant(sts)))  # DataFrame


def test_dissrep_freq():
    """Test representative selection with frequency criterion."""
    data = np.array(
        [
            ["A", "B"],
            ["A", "B"],
            ["A", "B"],
            ["B", "C"],
        ]
    )
    sts = seqdef(data)
    D = seqdist(sts, method="HAM")
    result = dissrep(D, criterion="freq", nrep=2)

    # Should select the most frequent sequences
    assert len(result["representatives"]) == 2
    # First representative should be one of the A-B duplicates
    assert result["representatives"][0] in [0, 1, 2]


def test_dissrep_dist():
    """Test representative selection with distance criterion."""
    data = np.array(
        [
            ["A", "A", "A"],
            ["A", "B", "A"],
            ["B", "B", "B"],
        ]
    )
    sts = seqdef(data)
    D = seqdist(sts, method="HAM")
    result = dissrep(D, criterion="dist", nrep=1)

    # With dist criterion, select sequence with minimum sum of distances
    assert len(result["representatives"]) == 1
    # Middle sequence (1) has smallest sum of distances to others
    assert result["representatives"][0] == 1


def test_dissrep_with_weights():
    """Test representative selection with sequence weights."""
    data = np.array(
        [
            ["A", "B"],
            ["A", "B"],
            ["B", "C"],
        ]
    )
    sts = seqdef(data, weights=[1.0, 1.0, 3.0])
    D = seqdist(sts, method="HAM")
    result = dissrep(D, criterion="freq", weights=sts.weights, nrep=1)

    # With high weight on sequence 2, it might be selected
    # (depends on frequency calculation with weights)
    assert len(result["representatives"]) == 1
    assert "metadata" in result
    assert result["metadata"]["weights_sum"] == 5.0


def test_dissrep_md_uses_weighted_denominator():
    """Per-representative MD should be weighted by assigned sequence weights."""
    D = np.array(
        [
            [0.0, 1.0, 4.0, 4.0],
            [1.0, 0.0, 4.0, 4.0],
            [4.0, 4.0, 0.0, 1.0],
            [4.0, 4.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    weights = np.array([10.0, 1.0, 10.0, 1.0], dtype=float)
    result = dissrep(
        D,
        criterion="dist",
        nrep=2,
        pradius=0.0,
        dmax=4.0,
        weights=weights,
    )

    reps = result["representatives"]
    rep_group = result["rep_group"]
    stats = result["statistics"]

    assert len(reps) == 2

    for i, rep_idx in enumerate(reps):
        sel = rep_group == i
        d = D[sel][:, rep_idx]
        w = weights[sel]
        expected_md = float((d * w).sum() / w.sum())
        assert np.isclose(float(stats.loc[f"r{i+1}", "MD"]), expected_md)
        # V[i] = mean(disscenter(D[sel, sel], weights=w[sel])) -- R-faithful
        # within-group inertia, not the global-centre weighted mean.
        from traminex import disscenter

        sub = D[np.ix_(np.flatnonzero(sel), np.flatnonzero(sel))]
        expected_v = float(np.mean(disscenter(sub, weights=w)))
        assert np.isclose(float(stats.loc[f"r{i+1}", "V"]), expected_v)

    expected_total_md = float((result["distances"] * weights).sum() / weights.sum())
    assert np.isclose(float(stats.loc["Total", "MD"]), expected_total_md)
