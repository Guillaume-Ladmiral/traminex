"""Direct calls to the private compiled kernels are validated and release the GIL (TX-AUD-20261011-045).

``traminex._pb_kernels`` is private, and the dispatcher validates what it passes, but a direct call used to accept an
infinite, NaN or negative indel (4, NaN and -8 for one pair), ignore a NaN substitution cost, raise RuntimeError for
1-D input, and hold the GIL for the whole call (another Python thread ran at 2.4% of its idle rate). The binding now
refuses those inputs with ValueError and releases the GIL around the kernels. It checks only what a kernel reads: a
zero indel stays valid (the dispatcher passes one for an all-zero TRATE matrix), and so does a non-finite cost
between states that do not occur, as in the dispatcher's own check.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

PB = pytest.importorskip("traminex._pb_kernels", reason="compiled _pb_kernels not built")
ROOT = Path(__file__).resolve().parents[1]
CODES = np.array([[0, 1, 0, 1], [1, 0, 1, 0], [0, 0, 1, 1]], dtype=np.int32)
MATRIX = np.array([[0.0, 2.0], [2.0, 0.0]])
METHODS = ("ham_unit", "om", "om_mat", "lcs", "lcp", "rlcp", "dhd")


def _call(method: str, codes: np.ndarray, *, indel: float = 1.0, subm=None, freq=None):
    kernel = getattr(PB, "pairwise_" + method)
    if method == "om":
        return kernel(codes, indel, 2.0 if subm is None else subm)
    if method == "om_mat":
        return kernel(codes, indel, MATRIX if subm is None else subm)
    if method == "dhd":
        width = codes.shape[-1] if codes.ndim else 1
        return kernel(codes, np.full((width, 2), 0.5) if freq is None else freq, 0)
    return kernel(codes)


@pytest.mark.parametrize("method", ["om", "om_mat"])
@pytest.mark.parametrize(
    "indel", [np.inf, -np.inf, np.nan, -1.0], ids=["inf", "-inf", "nan", "negative"]
)
def test_an_indel_that_is_not_finite_and_non_negative_is_refused(method: str, indel: float) -> None:
    with pytest.raises(ValueError, match="indel"):
        _call(method, CODES, indel=indel)


@pytest.mark.parametrize("method", ["om", "om_mat"])
def test_a_zero_indel_stays_valid(method: str) -> None:
    # Deleting and inserting is free, so every distance is 0.
    np.testing.assert_array_equal(_call(method, CODES, indel=0.0), np.zeros((3, 3)))


@pytest.mark.parametrize("dtype", [np.int32, np.uint64], ids=["signed", "unsigned"])
@pytest.mark.parametrize("subm", [np.nan, np.inf, -1.0], ids=["nan", "inf", "negative"])
def test_a_scalar_cost_that_is_not_finite_and_non_negative_is_refused(subm: float, dtype) -> None:
    # Unsigned codes are valid input (test_remediation_native_codes.py); the cost is refused with them too.
    with pytest.raises(ValueError, match="cost"):
        _call("om", CODES.astype(dtype), subm=subm)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1.0], ids=["nan", "inf", "negative"])
def test_a_read_matrix_cost_that_is_not_finite_and_non_negative_is_refused(bad: float) -> None:
    matrix = MATRIX.copy()
    matrix[0, 1] = bad
    with pytest.raises(ValueError, match="cost"):
        _call("om_mat", CODES, subm=matrix)


def test_a_cost_involving_a_state_that_does_not_occur_is_not_checked() -> None:
    matrix = np.full((3, 3), np.inf)
    matrix[:2, :2] = MATRIX
    matrix[2, 2] = 0.0
    np.testing.assert_array_equal(_call("om_mat", CODES, subm=matrix), _call("om_mat", CODES))


@pytest.mark.parametrize(
    "bad", [np.nan, np.inf, -0.5, 1.5], ids=["nan", "inf", "negative", "above-1"]
)
def test_a_read_dhd_frequency_outside_0_and_1_is_refused(bad: float) -> None:
    freq = np.full((CODES.shape[1], 2), 0.5)
    freq[0, 1] = bad
    with pytest.raises(ValueError, match="freq"):
        _call("dhd", CODES, freq=freq)


@pytest.mark.parametrize("method", METHODS)
def test_one_dimensional_input_is_a_value_error(method: str) -> None:
    with pytest.raises(ValueError, match="2D"):
        _call(method, np.array([0, 1, 0], dtype=np.int32))


_THREAD_PROBE = """
import json, sys, threading, time
import numpy as np
sys.path.insert(0, {root!r})
from traminex import _pb_kernels as pb

codes = np.random.default_rng(0).integers(0, 4, size=(260, 60)).astype(np.int32)
pb.pairwise_om(codes[:8], 1.0, 2.0)  # warm up

def counted(action):
    count, stop = [0], [False]
    def spin():
        while not stop[0]:
            count[0] += 1
    worker = threading.Thread(target=spin)
    worker.start()
    time.sleep(0.05)
    count[0] = 0
    start = time.perf_counter()
    action()
    elapsed = time.perf_counter() - start
    done = count[0]
    stop[0] = True
    worker.join()
    return done / elapsed, elapsed

idle, _ = counted(lambda: time.sleep(0.5))
busy, seconds = counted(lambda: pb.pairwise_om(codes, 1.0, 2.0))
print(json.dumps({{"ratio": busy / idle, "seconds": seconds}}))
"""


def test_another_python_thread_runs_during_a_long_call() -> None:
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    done = subprocess.run(
        [sys.executable, "-c", _THREAD_PROBE.format(root=str(ROOT))],
        capture_output=True,
        text=True,
        check=True,
        cwd=ROOT,
        env=env,
    )
    result = json.loads(done.stdout.strip().splitlines()[-1])
    assert result["seconds"] > 0.1, result  # long enough to measure
    assert result["ratio"] > 0.25, result  # held, the GIL let the thread run at about 2%
