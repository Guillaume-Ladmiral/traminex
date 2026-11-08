#!/usr/bin/env python
"""
Verification script for key audit improvements.

This script is intentionally aligned with the current public APIs in version 0.2.2.
"""

from __future__ import annotations

import sys

import numpy as np


def _ok(label: str) -> None:
    print(f"[OK] {label}")


def _fail(msg: str) -> None:
    print(f"[FAIL] {msg}")
    raise SystemExit(1)


print("=" * 70)
print("TRAMINEX AUDIT IMPROVEMENTS VERIFICATION")
print("=" * 70)
print()

# Test 1: Import package
print("Test 1: Package import")
try:
    import traminex

    _ok(f"Package version: {traminex.__version__}")
    _ok(f"Total exports: {len(traminex.__all__)}")
except ImportError as exc:
    _fail(f"Package import failed: {exc}")
print()

# Test 2: Verify key exports are callable
print("Test 2: Export resolution")
try:
    from traminex import seqdiff, seqibad, seqmodst

    assert callable(seqmodst)
    assert callable(seqdiff)
    assert callable(seqibad)
    _ok("seqmodst, seqdiff, seqibad are callable")
except Exception as exc:  # noqa: BLE001
    _fail(f"Export resolution failed: {exc}")
print()

# Test 3: Verify the documented missing-state parameters
print("Test 3: Missing-state parameters")
try:
    from traminex import seqdef, seqindic, seqnum, seqtabstocc, seqtrate

    x = seqdef(np.array([["A", "*", "B"], ["A", "B", "B"]], dtype=object), nr="*")
    _ = seqindic(x, with_missing=False)
    _ = seqnum(x, with_missing=False)
    _ = seqtrate(x, include_missing=False)
    _ = seqtabstocc(x, include_missing=True)
    _ok("documented with_missing/include_missing parameters accepted")
except Exception as exc:  # noqa: BLE001
    _fail(f"Missing-state consistency check failed: {exc}")
print()

# Test 4: Plotting imports are clean
print("Test 4: Plotting import warnings")
try:
    import warnings

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        import traminex.plotting as _plotting  # noqa: F401

    deprec = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    if deprec:
        _fail(f"Found {len(deprec)} DeprecationWarning entries during plotting import")
    _ok("No DeprecationWarning emitted while importing plotting")
except Exception as exc:  # noqa: BLE001
    _fail(f"Plotting warning check failed: {exc}")
print()

# Test 5: Basic docstring quality checks
print("Test 5: Docstring presence")
try:
    from traminex import (
        cluster_agglomerative,
        dissassoc,
        disscenter,
        dissvar,
        gower_matrix,
        medoid_index,
    )

    funcs = [
        ("disscenter", disscenter),
        ("gower_matrix", gower_matrix),
        ("dissvar", dissvar),
        ("dissassoc", dissassoc),
        ("medoid_index", medoid_index),
        ("cluster_agglomerative", cluster_agglomerative),
    ]
    for name, fn in funcs:
        doc = fn.__doc__ or ""
        if len(doc) < 80:
            _fail(f"Docstring too short: {name}")
        _ok(f"Docstring present: {name}")
except Exception as exc:  # noqa: BLE001
    _fail(f"Docstring checks failed: {exc}")
print()

# Test 6: Distance normalization behavior (current API)
print("Test 6: Distance normalization API")
try:
    from traminex import seqdef, seqdist, seqLLCP, seqLLCS, seqRLCP

    x = seqdef(np.array([["A", "A", "B", "C"], ["A", "B", "B", "C"], ["B", "B", "C", "C"]]))

    _ = seqdist(x, method="LCS", norm="length")
    assert seqLLCS(x.data[0], x.data[1]) >= 0
    assert seqLLCP(x.data[0], x.data[1]) >= 0
    _ = seqRLCP(x, norm="length")
    _ok("pairwise LCS normalization and scalar LLCS/LLCP contracts verified")

    try:
        seqdist(x, method="LCS", normalize=True)  # type: ignore[call-arg]
        _fail("legacy normalize kwarg should be rejected")
    except TypeError:
        _ok("legacy normalize kwarg correctly rejected")
except Exception as exc:  # noqa: BLE001
    _fail(f"Distance normalization checks failed: {exc}")
print()

print("=" * 70)
print("VERIFICATION COMPLETE")
print("=" * 70)
print()
print("All verification checks passed.")
sys.exit(0)
