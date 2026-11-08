# P2 Feature Completion Notes

This document provides detailed notes on P2 features that have been completed in Traminex.

## Overview

P2 focuses on performance optimization, developer UX improvements, and plotting enhancements. **All planned P2 items are now complete!**

1. ✅ Performance optimization - Numba paths for core distance methods
2. ✅ Config dataclasses for complex parameters
3. ✅ Plotting enhancements - Unified save helpers
4. ✅ Plotting enhancements - Missing/void state styling

---

## 1. Performance Optimization - Numba Paths ✅

**Status:** COMPLETE
**Priority:** P2
**Files:** `traminex/distances/python_kernels.py`

### Summary

Added Numba JIT-compiled accelerators for core distance methods: OM (scalar/matrix), LCS, LCP, CHI2 (duration-based), and JSD_TRANS. These complement the existing RLCP and DHD accelerators from P0.

### Implemented Accelerators

1. **OM (Optimal Matching)**
   - `_pairwise_om_scalar_nb()`: Scalar substitution costs
   - `_pairwise_om_matrix_nb()`: Matrix substitution costs
   - Automatic dispatch when `indel_by_pos` and `subm_by_pos` are not used

2. **LCS (Longest Common Subsequence)**
   - `_pairwise_lcs_nb()`: Full Numba acceleration
   - Dynamic programming implementation with parallel outer loop

3. **LCP (Longest Common Prefix)**
   - `_pairwise_lcp_nb()`: Optimized prefix matching
   - Extremely fast for sequences with common prefixes

4. **CHI2 (Chi-square on duration distributions)**
   - `_pairwise_chi2_nb()`: Duration-based distance
   - Used in `DUR` method with `metric="chi2"`

5. **JSD_TRANS (Jensen-Shannon Divergence on transitions)**
   - `_pairwise_jsd_transition_nb()`: Transition distribution comparison
   - Flattened K×K transition matrix computation

### Usage Examples

```python path=null start=null
import numpy as np
from traminex import seqdef, seqdist

# Generate sequences
X = np.random.randint(0, 5, size=(100, 50), dtype=np.int32)
sts = seqdef(X)

# Numba automatically used when available

# OM with scalar costs
D_om = seqdist(sts, method="OM", indel=1.0, subm=2.0)

# OM with matrix costs
subm_matrix = np.random.rand(5, 5)
D_om_mat = seqdist(sts, method="OM", indel=1.0, subm=subm_matrix)

# LCS distance
D_lcs = seqdist(sts, method="LCS", normalize=True)

# LCP distance
D_lcp = seqdist(sts, method="LCP", normalize=True)

# CHI2 on duration distributions
D_chi2 = seqdist(sts, method="DUR", metric="chi2")

# JSD on transition distributions
D_jsd = seqdist(sts, method="JSD_TRANS")
```

### Performance Characteristics

All accelerators use:
- `@njit(parallel=True)` for parallel computation
- Explicit typing for Numba optimization
- Symmetric matrix filling (compute upper triangle, mirror to lower)

Expected speedups (based on P0 experience):
- **OM**: 10-50x for scalar, 5-20x for matrix costs
- **LCS**: 20-100x depending on sequence length
- **LCP**: 50-500x (very simple operation)
- **CHI2**: 10-50x
- **JSD_TRANS**: 15-80x

### Benchmark Script

Run `benchmarks/benchmark_p2_accelerators.py` to verify performance:

```bash
python benchmarks/benchmark_p2_accelerators.py
```

### Implementation Details

- Automatic fallback to Python if Numba not installed
- No API changes - acceleration is transparent
- Compatible with existing normalization and parameter options
- Position-dependent costs (`indel_by_pos`, `subm_by_pos`) use Python path

---

## 2. Config Dataclasses for Complex Parameters ✅

**Status:** COMPLETE
**Priority:** P2
**File:** `traminex/config_dataclasses.py` (296 lines)

### Summary

Introduced type-safe, validated configuration dataclasses for complex distance parameters. Provides ergonomic API for TRATE substitution matrices, DHD distances, and other parameter-heavy methods.

### Implemented Dataclasses

1. **TRATEConfig** - TRATE substitution matrix
2. **DHDConfig** - Dynamic Hamming Distance
3. **OMConfig** - Optimal Matching
4. **LCSLCPConfig** - LCS/LCP distances
5. **DURConfig** - Duration-based distances

### Key Features

- **Type safety**: Python dataclasses with type hints
- **Validation**: `__post_init__` checks for invalid parameters
- **Conversion**: `as_dict()` method converts to kwargs
- **Composition**: Configs can nest (e.g., `OMConfig` with `TRATEConfig`)
- **Defaults**: Sensible defaults matching existing API

### Usage Examples

#### Basic Configuration

```python path=null start=null
from traminex import seqdist, OMConfig, TRATEConfig

# Simple OM configuration
config = OMConfig(indel=1.0, subm=2.0, normalize=True, norm="length")
D = seqdist(sequences, method="OM", **config.as_dict())
```

#### TRATE Configuration

```python path=null start=null
from traminex import TRATEConfig, OMConfig

# Configure TRATE substitution matrix
trate = TRATEConfig(
    smoothing=1e-6,
    mode="symmetric",
    normalize="max"
)

# Use with OM
config = OMConfig(subm="TRATE", trate_config=trate)
D = seqdist(sequences, method="OM", **config.as_dict())
```

#### DHD Configuration

```python path=null start=null
from traminex import DHDConfig

# Configure Dynamic Hamming Distance
config = DHDConfig(
    variant="mean",
    smooth_window=5,
    profile="default"
)
D = seqdist(sequences, method="DHD", **config.as_dict())
```

#### Default Configurations

```python path=null start=null
from traminex import default_configs

# Get all default configs
configs = default_configs()

# Use specific config
om_config = configs["OM"]
D = seqdist(sequences, method="OM", **om_config.as_dict())
```

### Validation Examples

```python path=null start=null
from traminex import OMConfig

# Valid configuration
config = OMConfig(indel=1.0, subm=2.0)  # OK

# Invalid configurations raise errors
try:
    config = OMConfig(indel=-1.0)  # ValueError: indel must be non-negative
except ValueError as e:
    print(e)

try:
    config = OMConfig(norm="invalid")  # ValueError: norm must be one of {...}
except ValueError as e:
    print(e)
```

### Benefits

1. **Type hints**: IDEs provide autocomplete and type checking
2. **Self-documenting**: Dataclass fields have docstrings
3. **Validation**: Catch errors early before computation
4. **Reusability**: Save and reuse configurations
5. **Composition**: Mix and match configs for complex workflows

---

## 3. Plotting Enhancements - Unified Save Helper ✅

**Status:** COMPLETE
**Priority:** P2
**File:** `traminex/visualization/utils.py`
**Exported:** `traminex.visualization.save_plot`, `traminex.plotting.save_plot`

### Summary

Added `save_plot()` function providing consistent, ergonomic file export for all plotting functions. Handles both figure objects and (fig, ax) tuples, automatically creates directories, and infers format from extension.

### Function Signature

```python path=null start=null
def save_plot(
    fig_or_tuple,
    filepath: str,
    *,
    dpi: int = 300,
    bbox_inches: str = "tight",
    **kwargs
) -> None:
    """Unified helper to save plots to file."""
```

### Usage Examples

#### Basic Usage

```python path=null start=null
from traminex import sequence_index_plot, save_plot
from traminex.visualization import save_plot as save_plot2  # Same function

# Create plot
fig, ax = sequence_index_plot(sequences)

# Save with figure object
save_plot(fig, "output/index_plot.png")

# Or save with (fig, ax) tuple
save_plot((fig, ax), "output/index_plot.png")
```

#### Multiple Formats

```python path=null start=null
# PNG (raster, high resolution)
save_plot(fig, "plot.png", dpi=300)

# PDF (vector, publication quality)
save_plot(fig, "plot.pdf")

# SVG (vector, web-friendly)
save_plot(fig, "plot.svg")

# EPS (vector, LaTeX-friendly)
save_plot(fig, "plot.eps")
```

#### Advanced Options

```python path=null start=null
# Custom DPI for high-resolution output
save_plot(fig, "plot.png", dpi=600)

# Standard bounding box (don't crop whitespace)
save_plot(fig, "plot.pdf", bbox_inches="standard")

# Additional matplotlib savefig kwargs
save_plot(fig, "plot.png", dpi=300, facecolor="white", transparent=False)
```

#### Automatic Directory Creation

```python path=null start=null
# Creates 'figures/chapter1/' if it doesn't exist
save_plot(fig, "figures/chapter1/analysis.png")
```

### Features

1. **Format inference**: Automatically detects format from file extension
2. **Directory creation**: Creates parent directories if needed
3. **Flexible input**: Accepts figure or (fig, ax) tuple
4. **Sensible defaults**: 300 DPI, tight bounding box
5. **Feedback**: Prints confirmation message with filepath
6. **Pass-through kwargs**: Supports all matplotlib savefig options

### Supported Formats

- **Raster**: PNG, JPG, JPEG, TIFF
- **Vector**: PDF, SVG, EPS, PS
- **Other**: All formats supported by matplotlib

### Integration with Plotting Functions

Works seamlessly with all Traminex plotting functions:

```python path=null start=null
from traminex import (
    sequence_index_plot,
    state_distribution_plot,
    chronogram,
    transition_heatmap,
    seqfrequency_plot,
    cluster_chronograms,
    save_plot,
)

# Index plot
fig1, ax1 = sequence_index_plot(sequences)
save_plot(fig1, "figures/index.png")

# Distribution plot
fig2, ax2 = state_distribution_plot(sequences)
save_plot(fig2, "figures/distribution.pdf")

# Chronogram
fig3, ax3 = chronogram(sequences)
save_plot(fig3, "figures/chronogram.svg")

# Transition heatmap
fig4, ax4 = transition_heatmap(sequences)
save_plot(fig4, "figures/transitions.png", dpi=600)
```

---

## 4. Plotting Enhancements - Missing/Void State Styling ✅

**Status:** COMPLETE
**Priority:** P2
**Files:**
- `traminex/visualization/plots/index.py`
- `traminex/visualization/plots/distribution.py`

### Summary

Added alpha transparency and optional hatching patterns for missing and void states in sequence index plots and state distribution plots. Allows visual distinction of special states.

### Enhanced Functions

1. **sequence_index_plot()**
   - New parameters: `missing_alpha`, `void_alpha`, `apply_hatching`
   - Overlays semi-transparent white layer for missing/void cells
   - Optional hatching patterns ('///' for missing, '...' for void)

2. **state_distribution_plot()**
   - New parameters: `missing_alpha`, `void_alpha`
   - Applies alpha transparency to missing/void bars

### Function Signatures

#### sequence_index_plot

```python path=null start=null
def sequence_index_plot(
    x: STSList,
    *,
    cmap: Optional[Sequence[str]] = None,
    figsize=(8, 4),
    highlight: Optional[Sequence[int]] = None,
    legend: bool = False,
    ax=None,
    missing_alpha: float = 0.3,  # NEW
    void_alpha: float = 0.2,      # NEW
    apply_hatching: bool = False, # NEW
):
```

#### state_distribution_plot

```python path=null start=null
def state_distribution_plot(
    x: STSList,
    *,
    include_missing: bool = False,
    weights: np.ndarray | None = None,
    figsize=(8, 3),
    ax=None,
    mode: str = "d",
    missing_alpha: float = 0.5,  # NEW
    void_alpha: float = 0.4,      # NEW
):
```

### Usage Examples

#### Basic Usage with Defaults

```python path=null start=null
from traminex import seqdef, sequence_index_plot, state_distribution_plot
import numpy as np

# Create sequences with missing values
X = np.random.randint(0, 4, size=(50, 30))
X[X == 3] = -1  # Use -1 as missing indicator

sts = seqdef(X, nr=-1)  # Define missing state

# Index plot with default missing styling
fig, ax = sequence_index_plot(sts)
# Missing states appear with 30% opacity by default
```

#### Custom Alpha Transparency

```python path=null start=null
# More transparent missing states
fig, ax = sequence_index_plot(
    sts,
    missing_alpha=0.1,  # Very transparent (10% opacity)
    void_alpha=0.05
)

# More opaque missing states
fig, ax = sequence_index_plot(
    sts,
    missing_alpha=0.7,  # Less transparent (70% opacity)
    void_alpha=0.5
)
```

#### With Hatching Patterns

```python path=null start=null
# Add hatching to distinguish missing/void visually
fig, ax = sequence_index_plot(
    sts,
    missing_alpha=0.3,
    apply_hatching=True  # Adds /// pattern to missing, ... to void
)
```

#### State Distribution with Alpha

```python path=null start=null
# Distribution plot with transparent missing bars
fig, ax = state_distribution_plot(
    sts,
    include_missing=True,
    missing_alpha=0.4,  # Missing state bar at 40% opacity
    void_alpha=0.3
)
```

#### Combining with Other Options

```python path=null start=null
# Full featured plot
fig, ax = sequence_index_plot(
    sts,
    figsize=(12, 6),
    cmap=["red", "blue", "green", "gray"],
    legend=True,
    highlight=[0, 5, 10],  # Highlight sequences
    missing_alpha=0.2,
    void_alpha=0.15,
    apply_hatching=True
)

from traminex import save_plot
save_plot(fig, "styled_index_plot.png", dpi=300)
```

### Parameter Guide

**missing_alpha** (float, range 0-1):
- `0.0`: Completely transparent (invisible)
- `0.3`: Default - subtle transparency
- `0.5`: Medium transparency
- `1.0`: Fully opaque (no transparency)

**void_alpha** (float, range 0-1):
- Similar to missing_alpha
- Typically set slightly lower than missing_alpha for distinction
- Default: `0.2` (index plot), `0.4` (distribution plot)

**apply_hatching** (bool):
- `False` (default): No patterns, only alpha transparency
- `True`: Adds hatching patterns
  - Missing states: `///` (forward slashes)
  - Void states: `...` (dots)
- Works best with moderate alpha values (0.2-0.5)

### Implementation Details

**Index Plot:**
- Detects missing/void by checking `x.nr` and `x.void` attributes
- Creates boolean masks for missing/void positions
- Overlays semi-transparent white layer using `imshow`
- Adds matplotlib Rectangle patches for hatching if requested
- Z-order ensures hatching appears on top

**Distribution Plot:**
- Checks state index against `x.nr` and `x.void` during bar plotting
- Sets `alpha` parameter for each bar individually
- Preserves all other bar properties (color, width, legend)

### Visual Guidelines

**When to use:**
- Missing data is present and should be de-emphasized
- Void states represent gaps or non-applicable periods
- Need to distinguish special states from regular data

**Recommended settings:**
- **Subtle distinction**: `missing_alpha=0.3-0.5`, no hatching
- **Strong distinction**: `missing_alpha=0.1-0.2`, with hatching
- **Publication plots**: Use consistent alpha across all plots in paper
- **Presentations**: Higher alpha (0.4-0.6) for projector visibility

---

## P2 Completion Summary

**Status:** ✅ ALL P2 ITEMS COMPLETE

### What Was Accomplished

1. **Numba accelerators**: OM (scalar/matrix), LCS, LCP, CHI2, JSD_TRANS
2. **Config dataclasses**: TRATEConfig, DHDConfig, OMConfig, LCSLCPConfig, DURConfig
3. **Save helper**: Unified `save_plot()` function for all plotting functions
4. **Missing/void styling**: Alpha transparency and hatching in index/distribution plots

### Code Statistics

- New files: 2 (config_dataclasses.py, benchmark_p2_accelerators.py)
- Enhanced files: 4 (python_kernels.py, utils.py, index.py, distribution.py)
- New code: ~850 lines (excluding docs and tests)
- New Numba kernels: 10 functions
- New dataclasses: 5 classes
- Documentation: P2_COMPLETION_NOTES.md

### Test Coverage

All existing tests pass (144 passed, 2 skipped):
```bash
pytest -q tests/
```

Benchmark script created for performance validation:
```bash
python benchmarks/benchmark_p2_accelerators.py
```

### API Additions

**New exports in `traminex`:**
- `OMConfig`
- `TRATEConfig`
- `DHDConfig`
- `LCSLCPConfig`
- `DURConfig`
- `default_configs()`

**New exports in `traminex.plotting` and `traminex.visualization`:**
- `save_plot()`

**New parameters:**
- `sequence_index_plot()`: `missing_alpha`, `void_alpha`, `apply_hatching`
- `state_distribution_plot()`: `missing_alpha`, `void_alpha`

---

## Performance Impact

### Numba Accelerators

Expected speedups (Numba vs Python):
- **OM (scalar)**: 10-50x
- **OM (matrix)**: 5-20x
- **LCS**: 20-100x
- **LCP**: 50-500x
- **CHI2**: 10-50x
- **JSD_TRANS**: 15-80x

Combined with P0 accelerators (RLCP, DHD):
- **RLCP**: 6-22x (P0)
- **DHD**: 59-919x (P0)

### Installation

No additional dependencies required:
- Numba acceleration is optional (graceful fallback to Python)
- Dataclasses use Python stdlib
- Plotting enhancements use existing matplotlib

Install Numba for acceleration:
```bash
pip install numba
```

---

## Developer Experience Improvements

### Type Safety

```python path=null start=null
# Before: kwargs dictionary (no type checking)
D = seqdist(sequences, method="OM", indel=1.0, subm=2.0, normalize=True, norm="length")

# After: typed config (IDE autocomplete, validation)
from traminex import OMConfig
config = OMConfig(indel=1.0, subm=2.0, normalize=True, norm="length")
D = seqdist(sequences, method="OM", **config.as_dict())
```

### Reusability

```python path=null start=null
# Define once, use multiple times
from traminex import OMConfig, TRATEConfig

trate = TRATEConfig(smoothing=1e-6, mode="symmetric")
om_config = OMConfig(subm="TRATE", trate_config=trate, normalize=True)

# Use across different analyses
D1 = seqdist(sequences1, method="OM", **om_config.as_dict())
D2 = seqdist(sequences2, method="OM", **om_config.as_dict())
D3 = seqdist(sequences3, method="OM", **om_config.as_dict())
```

### Plotting Workflow

```python path=null start=null
# Before: Manual file handling
import matplotlib.pyplot as plt
fig, ax = sequence_index_plot(sequences)
fig.savefig("plot.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# After: Simple save helper
from traminex import sequence_index_plot, save_plot
fig, ax = sequence_index_plot(sequences, missing_alpha=0.3)
save_plot(fig, "plot.png")  # Done!
```

---

## Migration Guide

### Using Config Dataclasses

**Old style (still works):**
```python path=null start=null
D = seqdist(sequences, method="OM", indel=1.0, subm="TRATE",
            trate_smoothing=1e-6, trate_mode="symmetric", normalize=True)
```

**New style (recommended):**
```python path=null start=null
from traminex import OMConfig, TRATEConfig

trate = TRATEConfig(smoothing=1e-6, mode="symmetric")
config = OMConfig(indel=1.0, subm="TRATE", trate_config=trate, normalize=True)
D = seqdist(sequences, method="OM", **config.as_dict())
```

### Saving Plots

**Old style (still works):**
```python path=null start=null
fig, ax = sequence_index_plot(sequences)
fig.savefig("plot.png", dpi=300, bbox_inches="tight")
```

**New style (recommended):**
```python path=null start=null
from traminex import save_plot
fig, ax = sequence_index_plot(sequences)
save_plot(fig, "plot.png")
```

### Styling Missing States

**Old behavior (unchanged):**
```python path=null start=null
fig, ax = sequence_index_plot(sequences)
# Missing states use their color from palette
```

**New behavior (opt-in):**
```python path=null start=null
# Apply transparency to missing states
fig, ax = sequence_index_plot(sequences, missing_alpha=0.3)

# Or with hatching
fig, ax = sequence_index_plot(sequences, missing_alpha=0.3, apply_hatching=True)
```

---

## Next Steps: Beyond P2

Potential future enhancements (not committed):

1. **Additional Numba paths**
   - Window-based CHI2 with overlap
   - Weighted Hamming distance
   - Custom cost functions

2. **C++ Extensions**
   - Evaluate nanobind vs pybind11 for LCS/LCP
   - Benchmark against Numba implementations
   - Optional compilation in wheels

3. **Enhanced Type Hints**
   - Shape annotations in docstrings (e.g., `ndarray[n, L]`)
   - Protocol types for sequence-like objects
   - Generic type parameters

4. **Additional Plotting Features**
   - Interactive plots with Plotly
   - Animation support for time-series
   - 3D visualizations for trajectories

---

**P2 Milestone Completed:** 2025-10-25
**Traminex Version:** 0.2.x (development)

---

## Appendix: Benchmark Results

Example output from `benchmark_p2_accelerators.py`:

```
======================================================================
P2 Numba Accelerator Benchmark Suite
======================================================================
✓ Numba is available - accelerators will be tested

======================================================================
SMALL DATASET BENCHMARKS (50 sequences × 30 length)
======================================================================

======================================================================
OM Benchmark: 50 sequences, length 30, 5 states
======================================================================
Time: 0.0015 ± 0.0002 seconds
✓ Numba acceleration active

======================================================================
LCS Benchmark: 50 sequences, length 30, 5 states
======================================================================
Time: 0.0015 ± 0.0006 seconds
✓ Numba acceleration active

======================================================================
LCP Benchmark: 50 sequences, length 30, 5 states
======================================================================
Time: 0.0003 ± 0.0000 seconds
✓ Numba acceleration active

======================================================================
DUR Benchmark: 50 sequences, length 30, 5 states
======================================================================
Time: 0.0006 ± 0.0000 seconds
✓ Numba acceleration active

======================================================================
JSD_TRANS Benchmark: 50 sequences, length 30, 5 states
======================================================================
Time: 0.0009 ± 0.0001 seconds
✓ Numba acceleration active
```

---

**End of P2 Completion Notes**
