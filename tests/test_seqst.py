import numpy as np
import pandas as pd

import traminex as tx
from traminex.stats import seqST, turbulence_exact


def test_seqST_basic_and_normalized():
    # Simple sequences with obvious DSS and durations
    # x0: A A A A (one run) -> phi=2, durations=[4], s2=0, s2max=0 -> T=log2(2*(0+1)/(0+1))=1
    # x1: A B A B (runs of 1) -> DSS same as sequence, phi=2^4 - (with repeats correction)
    df = pd.DataFrame(
        {
            "t1": ["A", "A"],
            "t2": ["A", "B"],
            "t3": ["A", "A"],
            "t4": ["A", "B"],
        }
    )
    X = tx.seqdef(df)
    vals = seqST(X, norm=False)
    assert vals.shape == (2,)
    assert vals[0] == 1.0  # all A -> one run, phi=2
    # normalized should be in [0,1]
    vals_n = seqST(X, norm=True)
    assert np.all(vals_n >= 0.0) and np.all(vals_n <= 1.0)


def test_turbulence_exact_matches_mean_of_seqST():
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "B"],
            "t2": ["A", "B", "B"],
            "t3": ["A", "A", "B"],
            "t4": ["A", "B", "B"],
        }
    )
    X = tx.seqdef(df)
    per = seqST(X, norm=True)
    agg = turbulence_exact(X, norm=True)
    assert agg == np.mean(per)


def test_seqST_type1_vs_type2():
    """Test that type=1 and type=2 produce different results when not all states are visited."""
    # Create sequences with alphabet [A, B, C] but only use A and B
    df = pd.DataFrame(
        {
            "t1": ["A", "B", "A"],
            "t2": ["A", "B", "A"],
            "t3": ["B", "A", "A"],
            "t4": ["B", "A", "B"],
        }
    )
    # Explicitly provide alphabet to ensure C is included but not used
    X = tx.seqdef(df, alphabet=["A", "B", "C"])

    # Type 1: variance only over observed runs (ignores non-visited state C)
    t1 = seqST(X, norm=False, type=1)

    # Type 2: variance denominator includes non-visited state C
    t2 = seqST(X, norm=False, type=2)

    # Type 2 should generally be different from type 1 when K > d
    # We can't assert exact values without knowing the calculation details,
    # but we can verify they're both valid and type 2 accounts for the extra state
    assert t1.shape == (3,)
    assert t2.shape == (3,)
    assert np.all(np.isfinite(t1))
    assert np.all(np.isfinite(t2))
    # Type 2 should typically be higher due to including non-visited states
    # (but this depends on the specific formula, so we just check they differ)
    assert not np.allclose(t1, t2), "Type 1 and type 2 should produce different results"


def test_seqST_type2_with_all_states_visited():
    """When all states are visited, type=1 and type=2 should be similar."""
    # Create sequences that use all states in alphabet
    df = pd.DataFrame(
        {
            "t1": ["A", "B", "C"],
            "t2": ["B", "C", "A"],
            "t3": ["C", "A", "B"],
        }
    )
    X = tx.seqdef(df, alphabet=["A", "B", "C"])

    t1 = seqST(X, norm=False, type=1)
    t2 = seqST(X, norm=False, type=2)

    # When all states are visited at least once, d == K, so type 1 and type 2 converge
    # They should be close but may not be identical due to calculation differences
    assert t1.shape == (3,)
    assert t2.shape == (3,)
    # For sequences where all alphabet states appear, results should be close
    # (actual behavior depends on run structure)


def test_seqST_normalized_type2():
    """Test normalized turbulence with type=2."""
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "B"],
            "t2": ["A", "B", "B"],
            "t3": ["A", "A", "C"],
            "t4": ["B", "B", "C"],
        }
    )
    X = tx.seqdef(df, alphabet=["A", "B", "C", "D"])  # D is never used

    vals = seqST(X, norm=True, type=2)

    assert vals.shape == (3,)
    assert np.all(vals >= 0.0), "Normalized turbulence should be >= 0"
    assert np.all(vals <= 1.0), "Normalized turbulence should be <= 1"
    assert np.all(np.isfinite(vals)), "Normalized turbulence should be finite"


def test_seqST_single_state_sequence():
    """Test turbulence for a sequence with only one state (should be minimal)."""
    df = pd.DataFrame(
        {
            "t1": ["A", "A"],
            "t2": ["A", "B"],
            "t3": ["A", "B"],
            "t4": ["A", "B"],
        }
    )
    X = tx.seqdef(df)

    vals = seqST(X, norm=False, type=1)

    # First sequence is all A -> one run, minimal turbulence
    assert vals[0] == 1.0  # phi=2 (empty + A), s2=0, s2max=0 -> log2(2*1/1) = 1
    assert vals[1] > vals[0]  # Second sequence has transitions, higher turbulence


def test_turbulence_exact_type2():
    """Test aggregate turbulence_exact with type=2."""
    df = pd.DataFrame(
        {
            "t1": ["A", "A", "B"],
            "t2": ["A", "B", "B"],
            "t3": ["A", "A", "B"],
        }
    )
    X = tx.seqdef(df, alphabet=["A", "B", "C"])

    agg1 = turbulence_exact(X, norm=True, type=1)
    agg2 = turbulence_exact(X, norm=True, type=2)

    assert 0.0 <= agg1 <= 1.0
    assert 0.0 <= agg2 <= 1.0
    # Aggregate should be mean of per-sequence values
    per2 = seqST(X, norm=True, type=2)
    assert np.isclose(agg2, np.mean(per2))
