import numpy as np
import pytest

import traminex
from traminex import seqdef
from traminex.stats import seqtransn_count, turbulence_proxy


@pytest.mark.parametrize("length", [1023, 1024])
@pytest.mark.parametrize("dss", [False, True])
def test_subsequence_float_result_at_binary64_count_boundary(length, dss):
    codes = np.arange(length, dtype=np.int64).reshape(1, -1)
    codes.setflags(write=False)
    x = traminex.STSList(data=codes, alphabet=[f"S{i}" for i in range(length)])
    result = traminex.seqsubsn(x, dss=dss)
    expected = np.ldexp(1.0, 1023) if length == 1023 else np.inf
    # Every symbol is unique: each of the 2**length subsets has one subsequence.
    assert result.shape == (1,)
    assert result[0] == expected
    np.testing.assert_array_equal(x.data, codes)


@pytest.mark.parametrize("norm", [False, True])
def test_turbulence_ieee_result_for_overflowing_subsequence_count(norm):
    x = traminex.STSList(
        data=np.arange(1024).reshape(1, -1), alphabet=[f"S{i}" for i in range(1024)]
    )
    result = traminex.seqST(x, norm=norm)
    # All spell durations are 1: variance and maximum variance are zero.
    # R's raw log2(Inf) is Inf; normalized Inf/Inf is NaN.
    assert result.shape == (1,)
    if norm:
        assert np.isnan(result[0])
    else:
        assert np.isposinf(result[0])


def test_long_repeated_states_use_exact_recurrence_before_float_rounding():
    # For an alternating word of length n, the distinct count is F(n+3)-1.
    # This independent Fibonacci calculation establishes float overflow.
    previous, current = 0, 1
    for _ in range(1503):
        previous, current = current, previous + current
    assert previous - 1 > int(np.finfo(float).max)
    x = traminex.STSList(data=(np.arange(1500) % 2).reshape(1, -1), alphabet=["A", "B"])
    for dss in (False, True):
        assert np.isposinf(traminex.seqsubsn(x, dss=dss)[0])
    assert np.isposinf(traminex.seqST(x)[0])


def test_seqtransn_count_and_turbulence_proxy():
    mat = np.array(
        [
            [0, 0, 0, 0],  # no transitions
            [0, 1, 0, 1],  # transitions at every step
        ]
    )
    x = seqdef(mat)
    per = seqtransn_count(x, normalized=True)
    assert per.shape == (2,)
    assert per[0] == 0.0
    # 3 steps, changes at each step -> 1.0
    assert per[1] == 1.0
    # Global proxy: simple average
    t = turbulence_proxy(x, normalized=True)
    assert t == 0.5
