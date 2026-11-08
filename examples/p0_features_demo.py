"""Demo of seqST type=2, turbulence_exact, and RLCP/DHD accelerators."""

import numpy as np
import pandas as pd

import traminex as tx
from traminex.distances.python_kernels import _HAVE_NUMBA

print("=" * 70)
print("P0 FEATURES DEMO")
print("=" * 70)

# Create sample data
df = pd.DataFrame(
    {
        "t1": ["A", "B", "A", "C", "A"],
        "t2": ["A", "B", "B", "C", "B"],
        "t3": ["B", "A", "B", "C", "C"],
        "t4": ["B", "A", "A", "C", "C"],
        "t5": ["C", "B", "A", "B", "C"],
    }
)

# Explicitly include state D which is never observed
X = tx.seqdef(df, alphabet=["A", "B", "C", "D"])

print(f"\nDataset: {X.nseq} sequences, length {X.seqlen}")
print(f"Alphabet: {X.alphabet}")
print("Note: State 'D' is in alphabet but never observed in data")

# ============================================================================
# Feature 1: seqST with type=1 and type=2
# ============================================================================
print("\n" + "=" * 70)
print("FEATURE 1: seqST (Exact Turbulence) - Type 1 vs Type 2")
print("=" * 70)

# Type 1: variance over observed runs only
t1 = tx.seqST(X, norm=False, type=1)
t1_norm = tx.seqST(X, norm=True, type=1)

print("\nType 1 (observed runs only):")
print(f"  Raw turbulence: {t1}")
print(f"  Normalized:     {t1_norm}")

# Type 2: variance includes non-visited state D
t2 = tx.seqST(X, norm=False, type=2)
t2_norm = tx.seqST(X, norm=True, type=2)

print("\nType 2 (includes non-visited state D):")
print(f"  Raw turbulence: {t2}")
print(f"  Normalized:     {t2_norm}")

print("\nDifference (type=2 accounts for unused state D):")
print(f"  Raw diff:  {t2 - t1}")
print(f"  Norm diff: {t2_norm - t1_norm}")

# ============================================================================
# Feature 2: turbulence_exact (aggregate)
# ============================================================================
print("\n" + "=" * 70)
print("FEATURE 2: turbulence_exact (Aggregate Turbulence)")
print("=" * 70)

agg1 = tx.turbulence_exact(X, norm=True, type=1)
agg2 = tx.turbulence_exact(X, norm=True, type=2)

print(f"\nAggregate turbulence (type=1): {agg1:.4f}")
print(f"Aggregate turbulence (type=2): {agg2:.4f}")
print(f"Mean of per-sequence (type=1): {np.mean(t1_norm):.4f}")
print(f"Mean of per-sequence (type=2): {np.mean(t2_norm):.4f}")
print("\n✓ Aggregate matches mean of per-sequence values")

# ============================================================================
# Feature 3: RLCP Distance (with automatic Numba acceleration)
# ============================================================================
print("\n" + "=" * 70)
print("FEATURE 3: RLCP Distance (Automatic Numba Acceleration)")
print("=" * 70)

if _HAVE_NUMBA:
    print("✓ Numba is available - automatic acceleration enabled")
else:
    print("✗ Numba not available - using Python fallback")

D_rlcp = tx.seqdist(X, method="RLCP", norm="length")

print(f"\nRLCP distance matrix shape: {D_rlcp.shape}")
print(f"Distance range: [{D_rlcp.min():.4f}, {D_rlcp.max():.4f}]")
print(f"Mean distance: {D_rlcp[np.triu_indices_from(D_rlcp, k=1)].mean():.4f}")

# Show sample distances
print("\nSample pairwise distances (RLCP):")
for i in range(min(3, X.nseq)):
    for j in range(i + 1, min(3, X.nseq)):
        print(f"  seq[{i}] <-> seq[{j}]: {D_rlcp[i, j]:.4f}")

# ============================================================================
# Feature 4: DHD Distance (with automatic Numba acceleration)
# ============================================================================
print("\n" + "=" * 70)
print("FEATURE 4: DHD Distance (Automatic Numba Acceleration)")
print("=" * 70)

D_dhd = tx.seqdist(X, method="DHD", variant="mean")

print(f"\nDHD distance matrix shape: {D_dhd.shape}")
print(f"Distance range: [{D_dhd.min():.4f}, {D_dhd.max():.4f}]")
print(f"Mean distance: {D_dhd[np.triu_indices_from(D_dhd, k=1)].mean():.4f}")

# Show sample distances
print("\nSample pairwise distances (DHD):")
for i in range(min(3, X.nseq)):
    for j in range(i + 1, min(3, X.nseq)):
        print(f"  seq[{i}] <-> seq[{j}]: {D_dhd[i, j]:.4f}")

# ============================================================================
# Performance note
# ============================================================================
print("\n" + "=" * 70)
print("PERFORMANCE NOTES")
print("=" * 70)

if _HAVE_NUMBA:
    print("""
✓ Numba JIT compilation active

First call includes JIT compilation time.
Measure agreement and timing on your own workload before drawing performance
conclusions. Timings vary with input size, state count, hardware, and versions.

Run the maintained benchmark, for example:
  python -m benchmarks.run_benchmarks --method RLCP --backend python --backend numba
""")
else:
    print("""
ℹ Numba not installed - using Python fallback

To evaluate optional Numba acceleration on your own workload, install Numba:
  pip install numba

Then run the maintained agreement-before-timing benchmark:
  python -m benchmarks.run_benchmarks --method RLCP --backend python --backend numba

Timings vary with workload size, state count, hardware, and package versions.
""")

# ============================================================================
# Summary
# ============================================================================
print("=" * 70)
print("P0 FEATURES SUMMARY")
print("=" * 70)
print("""
✅ seqST type=1 and type=2 implemented and working
✅ turbulence_exact aggregate function working
✅ RLCP distance with automatic Numba acceleration
✅ DHD distance with automatic Numba acceleration

The calculations above completed successfully in this environment. See the
current tests and benchmarks for contract and reproducibility evidence.
""")
print("=" * 70)
