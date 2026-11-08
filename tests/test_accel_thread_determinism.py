"""The compiled kernels give byte-identical results at every thread count (W-06).

The compiled suite had only ever run at the default OpenMP thread count, so a
reduction whose order changed with the number of workers would have gone
unnoticed: it would show up here and nowhere else.

Structurally it cannot happen. ``cpp/src/accel_kernels.h:119`` opens one
``#pragma omp parallel`` region in which each worker takes a contiguous,
disjoint block of the pair space, keeps its own ``Workspace``, and writes only
its own output cells. No value is accumulated across threads, so no summation
order depends on the worker count. This module proves that empirically rather
than resting on the reading: each thread count runs in its own process --
OpenMP fixes its team size at start-up, so it cannot be varied in-process --
and the raw kernel outputs are compared byte for byte.

``accel_max_threads()`` is asserted to follow ``OMP_NUM_THREADS`` in every
child, so the comparison cannot pass merely because OpenMP was inactive.

Skips cleanly when the compiled kernel is not built, as
``tests/test_accel_byte_identity.py`` does.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap

import pytest

try:
    from traminex import _pb_kernels as PB
except ImportError:  # pragma: no cover - the pure-Python install
    PB = None

pytestmark = pytest.mark.skipif(PB is None, reason="compiled _pb_kernels not built")

THREAD_COUNTS = (1, 2, 4, 8)

# Corpora big enough that eight workers genuinely split the pair space: the
# smallest has 200 sequences, so 19,900 pairs over 8 blocks.
_CHILD = r"""
import hashlib, json, os, sys
import numpy as np

from scripts._checkout_imports import prefer_checkout_imports
prefer_checkout_imports()

import traminex
from traminex import _pb_kernels as PB
from traminex import distances as D

rng = np.random.default_rng(20260830)
mats = []
for n, L, K in [(200, 30, 4), (150, 45, 7), (300, 12, 3), (120, 60, 10)]:
    mats.append(rng.integers(0, K, size=(n, L), dtype=np.int32))

digest = hashlib.sha256()
arrays = 0
for index, mat in enumerate(mats):
    states = int(mat.max()) + 1
    sub = np.ones((states, states)) * 2.0 - 2.0 * np.eye(states)
    outputs = {
        "ham_unit": PB.pairwise_ham_unit(mat),
        "om": PB.pairwise_om(mat, 1.0, 2.0),
        "om_mat": PB.pairwise_om_mat(mat, 1.0, sub),
        "lcs": PB.pairwise_lcs(mat),
        # the dispatcher path as well, not only the raw kernels
        "dispatch_om": np.asarray(D.pairwise(mat, method="OM", indel=1.0, subm=2.0)),
        "dispatch_ham": np.asarray(D.pairwise(mat, method="HAM")),
    }
    for key in sorted(outputs):
        digest.update(f"{index}:{key}".encode())
        digest.update(np.ascontiguousarray(outputs[key], dtype=np.float64).tobytes())
        arrays += 1

print(json.dumps({
    "sha256": digest.hexdigest(),
    "arrays": arrays,
    "max_threads": int(PB.accel_max_threads()),
    "openmp": bool(PB.accel_openmp_enabled()),
    "omp_env": os.environ.get("OMP_NUM_THREADS"),
    "package": traminex.__file__,
}, sort_keys=True))
"""


def _run(threads: int) -> dict:
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = str(threads)
    environment["TRAMINEX_BINDING"] = "pybind11"
    completed = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(_CHILD)],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert completed.returncode == 0, completed.stderr
    return json.loads(completed.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def per_thread_results() -> dict[int, dict]:
    return {threads: _run(threads) for threads in THREAD_COUNTS}


def test_openmp_is_actually_active_and_follows_the_environment(per_thread_results):
    """Without this the byte comparison could pass with OpenMP disabled."""
    for threads, result in per_thread_results.items():
        assert result["openmp"] is True, threads
        assert result["omp_env"] == str(threads)
        assert result["max_threads"] == threads, (
            f"OMP_NUM_THREADS={threads} but accel_max_threads() reported "
            f"{result['max_threads']}"
        )


def test_every_child_used_this_checkout(per_thread_results):
    """Guards against a child importing an installed copy instead."""
    expected = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for threads, result in per_thread_results.items():
        assert result["package"].startswith(expected), (threads, result["package"])


@pytest.mark.parametrize("threads", THREAD_COUNTS[1:])
def test_kernel_output_is_byte_identical_to_the_single_threaded_run(threads, per_thread_results):
    single = per_thread_results[1]
    other = per_thread_results[threads]

    assert other["arrays"] == single["arrays"] == 24
    assert other["sha256"] == single["sha256"], (
        f"compiled kernels differ between 1 and {threads} threads; a reduction "
        "is reordering with the worker count"
    )
