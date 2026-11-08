# P1 (Priority 1) Partial Completion Summary

**Date**: 2025-10-25
**Status**: 🟡 PARTIALLY COMPLETE (2/6 major items)

## Overview

P1 focuses on feature breadth and robustness, expanding the project's capabilities beyond the core M1 baseline. Due to the extensive scope of P1, implementation has been prioritized based on impact and complexity.

## Completed Items ✅

### 1. Position-Dependent Substitution Matrix Support ✅

**Status**: COMPLETE

**What was done**:
- ✅ Verified existing `subm_by_pos` parameter support in `om_distance()` function
- ✅ Confirmed shape `[L,K,K]` array handling for time-specific substitution costs
- ✅ Python kernel already supports position-dependent substitution via `subm_by_pos`
- ✅ Implementation works correctly with DP algorithm (chooses time index by `min(i-1, j-1)`)

**Implementation**:
- `traminex/distances/python_kernels.py`: Lines 33-82
- Uses `sc_pos(i, j, t)` function for per-position substitution lookup

**Note**: Compiled backends (nanobind/pybind11) do not yet support `subm_by_pos` and will fall back to Python when this parameter is provided.

### 2. OMloc Wrapper Function ✅

**Status**: COMPLETE

**What was done**:
- ✨ Created comprehensive `omloc()` wrapper function
- ✨ Integrates `local_indel_by_pos()` for position-dependent indel costs
- ✨ Optionally uses `trate_substitution_matrix_timevarying()` for time-varying substitution
- ✨ Supports multiple indel methods: "change", "entropy"
- ✨ Supports multiple substitution modes: "TRATE", "CONSTANT", custom matrices, scalars
- ✨ Full parameter pass-through to underlying `pairwise()` function
- ✨ Three aliases: `omloc`, `OMloc`, `seqOMloc` (TraMineR compatibility)
- ✨ Comprehensive documentation with examples and usage notes
- ✨ 11 comprehensive tests, all passing

**Files Created**:
- `traminex/distances/omloc.py`: 197 lines with full implementation
- `tests/test_omloc.py`: 141 lines with 11 tests

**API**:
```python
omloc(seqmat,
      indel_method="change",           # "change" or "entropy"
      indel_base=1.0,
      indel_scale=1.0,
      indel_smooth=None,               # Optional smoothing window
      subm_timevarying=False,          # Position-specific substitution
      subm_base="TRATE",               # "TRATE", "CONSTANT", matrix, or scalar
      trate_smoothing=1e-9,
      trate_mode="symmetric",
      trate_normalize="none",
      norm="length",
      **kwargs)
```

**Examples**:
```python
import numpy as np
from traminex import omloc

seqmat = np.array([[0,0,1,1,2], [0,1,1,2,2], [2,2,1,0,0]])

# Position-dependent indel only (global substitution)
D1 = omloc(seqmat, subm_timevarying=False)

# Both position-dependent indel AND substitution
D2 = omloc(seqmat, subm_timevarying=True)

# Entropy-based indel localization
D3 = omloc(seqmat, indel_method="entropy")

# Custom substitution matrix
custom_subm = np.ones((3, 3)) * 1.5
np.fill_diagonal(custom_subm, 0.0)
D4 = omloc(seqmat, subm_base=custom_subm)
```

**Tests**:
- ✅ Basic functionality
- ✅ Time-varying vs. global substitution
- ✅ Different indel methods (change, entropy)
- ✅ Different substitution modes (TRATE, CONSTANT, custom)
- ✅ Scalar and matrix substitution
- ✅ Aliases (OMloc, seqOMloc)
- ✅ Normalization options
- ✅ Indel smoothing
- ✅ Larger datasets
- ✅ Condensed output format

## Remaining P1 Items 🔄

### 3. Extended seqefsub with Support/Confidence (TODO)

**Status**: NOT STARTED

**Scope**:
- Add `support_pct` calculation to frequent subsequence output
- Add `confidence` metric for rule generation
- Implement `lift` metric for association rules
- Simple rule generation from frequent subsequences

**Estimated Complexity**: Medium
**Priority**: Medium-High (useful for event sequence analysis)

### 4. Gap/Window Constraints for Events (TODO)

**Status**: NOT STARTED

**Scope**:
- Implement `min_gap` and `max_gap` parameters in `seqeconstraint`
- Extend `seqefsub` to respect gap constraints
- Add `id_weights` support for weighted subsequence mining
- Update event mining algorithms to handle temporal constraints

**Estimated Complexity**: Medium-High
**Priority**: Medium (enhances event analysis capabilities)

### 5. seqdiff with Bootstrapped CIs (TODO)

**Status**: NOT STARTED

**Scope**:
- Create `seqdiff()` function for group comparison
- Compute per-position state distribution differences
- Optional bootstrap resampling for confidence intervals
- Multiple comparison correction options
- Visualization helpers for difference plots

**Estimated Complexity**: Medium
**Priority**: High (important for group comparison analysis)

**Proposed API**:
```python
seqdiff(x, group_var,
        method="chi2",           # "chi2", "residuals", "prop_diff"
        bootstrap=False,
        n_bootstrap=1000,
        ci_level=0.95,
        correction="bonferroni")
```

### 6. Multi-Domain Weight Search (TODO)

**Status**: NOT STARTED

**Scope**:
- Weight grid search over multi-domain (multi-channel) distances
- Silhouette coefficient as optimization target
- Support for multiple distance metrics per domain
- Efficient evaluation over parameter grid
- Return optimal weights and corresponding silhouette score

**Estimated Complexity**: High
**Priority**: Medium-Low (specialized use case)

**Proposed API**:
```python
search_domain_weights(domains,           # List of STSList
                      labels,            # Cluster labels
                      weight_grid=None,  # Grid of weight combinations
                      method="OM",       # Distance method per domain
                      metric="silhouette")
```

## Technical Notes

### Position-Dependent Costs

The existing infrastructure already supports:
1. **Position-dependent indel**: via `indel_by_pos` array `[L]`
2. **Position-dependent substitution**: via `subm_by_pos` array `[L,K,K]`
3. **Time-varying TRATE**: via `trate_substitution_matrix_timevarying()`

The DP algorithm in `om_distance()` correctly handles both:
- Different indel costs at each position for gap penalties
- Different substitution costs at each position for state transitions

### OMloc Design Decisions

1. **Indel Method Selection**:
   - `"change"`: Recommended for life-course data with non-uniform transitions
   - `"entropy"`: Alternative based on state diversity at each position

2. **Substitution Options**:
   - `subm_timevarying=False`: Single global matrix (faster, simpler)
   - `subm_timevarying=True`: Position-specific matrices (slower, more accurate)

3. **Performance Considerations**:
   - Time-varying substitution requires storing L matrices of size K×K
   - Computation time increases linearly with sequence length
   - Memory usage: O(L×K²) for time-varying substitution

### Import Circular Dependency Fix

The `omloc` module initially had a circular import issue with `pairwise`. This was resolved by:
- Importing `pairwise` inside the `omloc()` function body
- Avoiding top-level imports from the same package's `__init__.py`

This pattern should be followed for future distance utility modules.

## Files Modified/Created

### Core Implementation
- ✨ `traminex/distances/omloc.py`: New OMloc wrapper module (197 lines)
- ✏️ `traminex/distances/__init__.py`: Added omloc exports
- ✏️ `traminex/__init__.py`: Added omloc, OMloc, seqOMloc to public API

### Tests
- ✨ `tests/test_omloc.py`: 11 comprehensive tests (141 lines)

### Documentation
- ✨ `P1_PARTIAL_COMPLETION.md`: This file

## Test Results

```
pytest tests/test_omloc.py -v
11 passed in 3.54s ✅

Full test suite:
pytest -q
81 passed, 2 skipped ✅
```

## Recommendations

### For Immediate Next Steps

1. **seqdiff Implementation** (Highest Priority):
   - Most requested feature for comparative analysis
   - Relatively self-contained implementation
   - High impact for researchers comparing groups

2. **Extended seqefsub** (Medium Priority):
   - Completes the event mining pipeline
   - Adds association rule mining capabilities
   - Moderate complexity, high value for event analysis

3. **Gap Constraints** (Medium Priority):
   - Enhances event pattern mining
   - Builds on existing seqefsub infrastructure
   - Useful for temporal pattern discovery

4. **Multi-Domain Weights** (Lower Priority):
   - Specialized use case
   - Can be implemented as separate utility module
   - Users can manually search weights in the meantime

### For Users

**Using OMloc**:
```python
import traminex as tx
import pandas as pd

# Load your sequence data
df = pd.read_csv("sequences.csv")
X = tx.seqdef(df)

# Basic OMloc with position-dependent indel
D = tx.omloc(X.data)

# Advanced: time-varying substitution
D_tv = tx.omloc(X.data, subm_timevarying=True)

# Customize localization
D_custom = tx.omloc(X.data,
                    indel_method="entropy",
                    indel_scale=0.8,
                    indel_smooth=5)
```

**Performance Tips**:
1. Use `subm_timevarying=False` for faster computation (default)
2. Enable time-varying substitution only when position-specific transitions matter
3. Smooth indel profiles (`indel_smooth=3-7`) for more stable costs
4. Normalize with `norm="length"` for interpretable distances

## Next Steps for P1 Completion

To complete P1, the following work remains:

1. **Implement seqdiff** (~200 lines + tests)
   - Per-position distribution comparison
   - Bootstrap confidence intervals
   - Multiple comparison correction
   - Visualization utilities

2. **Extend seqefsub** (~150 lines + tests)
   - Support percentage calculation
   - Confidence and lift metrics
   - Association rule generation

3. **Add gap constraints** (~100 lines + tests)
   - min_gap/max_gap in seqeconstraint
   - Update mining algorithms
   - id_weights support

4. **Multi-domain utilities** (~200 lines + tests)
   - Grid search implementation
   - Silhouette optimization
   - Results reporting

**Estimated Total**: ~650 lines of code + ~400 lines of tests + documentation

## Conclusion

**P1 is ~33% complete** (2/6 major items) with the most complex distance-related features implemented:

- ✅ Position-dependent substitution matrices (infrastructure verified)
- ✅ OMloc wrapper (comprehensive implementation with 11 tests)
- ⏳ Extended event mining features (pending)
- ⏳ Group comparison utilities (pending)
- ⏳ Multi-domain optimization (pending)

The foundation for advanced distance computation is now in place. The remaining items focus on statistical analysis, event mining, and multi-domain extensions.
