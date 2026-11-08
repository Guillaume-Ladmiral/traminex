import numpy as np
import pytest

import traminex.distances as dist


class _StubCompiled:
    def __init__(self):
        self.called_rlcp = False
        self.called_dhd = False
        self.called_lcs = False
        self.called_om = False
        self.last_freq_shape = None
        self.last_variant = None

    def pairwise_rlcp(self, seqmat):
        self.called_rlcp = True
        n = seqmat.shape[0]
        out = np.zeros((n, n), dtype=float)
        if n >= 2:
            out[0, 1] = out[1, 0] = 2.0
        return out

    def pairwise_dhd(self, seqmat, freq, variant):
        self.called_dhd = True
        self.last_freq_shape = tuple(freq.shape)
        self.last_variant = int(variant)
        n = seqmat.shape[0]
        return np.zeros((n, n), dtype=float)

    def pairwise_lcs(self, seqmat):
        self.called_lcs = True
        return np.zeros((seqmat.shape[0], seqmat.shape[0]), dtype=float)

    def pairwise_om(self, seqmat, indel, subm):
        self.called_om = True
        return np.zeros((seqmat.shape[0], seqmat.shape[0]), dtype=float)


def test_compiled_backend_used_for_rlcp(monkeypatch):
    stub = _StubCompiled()
    monkeypatch.setattr(dist, "_backend", "pybind11")
    monkeypatch.setattr(dist, "_pb", stub)

    seqmat = np.array([[0, 1, 2, 3], [1, 1, 2, 3], [0, 0, 2, 3]], dtype=np.int32)
    D = dist.pairwise(seqmat, method="RLCP")

    assert D.shape == (3, 3)
    assert stub.called_rlcp, "compiled RLCP kernel was not used"


def test_compiled_backend_used_for_dhd(monkeypatch):
    stub = _StubCompiled()
    monkeypatch.setattr(dist, "_backend", "pybind11")
    monkeypatch.setattr(dist, "_pb", stub)

    seqmat = np.array([[0, 1, 2, 0], [1, 1, 2, 2], [0, 2, 2, 1]], dtype=np.int32)
    D = dist.pairwise(seqmat, method="DHD", variant="sum")

    assert D.shape == (3, 3)
    assert stub.called_dhd, "compiled DHD kernel was not used"
    assert stub.last_variant == 1, "DHD variant='sum' should map to variant code 1"
    assert stub.last_freq_shape == (4, 3), "expected freq shape (L, K)"


@pytest.mark.parametrize("smooth_window", [2, 4])
def test_even_dhd_smoothing_preserves_position_count_for_python_backend(monkeypatch, smooth_window):
    monkeypatch.setattr(dist, "_backend", "python")
    monkeypatch.setattr(dist, "_pb", None)
    seqmat = np.array(
        [[0, 1, 2, 0, 1], [1, 1, 2, 2, 0], [0, 2, 2, 1, 1]],
        dtype=np.int32,
    )

    result = dist.pairwise(
        seqmat,
        method="DHD",
        variant="sum",
        smooth_window=smooth_window,
    )

    assert result.shape == (3, 3)
    assert np.all(np.isfinite(result))


def test_even_dhd_smoothing_reaches_compiled_backend_with_original_position_count(monkeypatch):
    stub = _StubCompiled()
    monkeypatch.setattr(dist, "_backend", "pybind11")
    monkeypatch.setattr(dist, "_pb", stub)
    seqmat = np.array(
        [[0, 1, 2, 0, 1], [1, 1, 2, 2, 0], [0, 2, 2, 1, 1]],
        dtype=np.int32,
    )

    result = dist.pairwise(seqmat, method="DHD", variant="mean", smooth_window=2)

    assert result.shape == (3, 3)
    assert stub.called_dhd
    assert stub.last_freq_shape == (5, 3)


@pytest.mark.parametrize("method", ["LCS", "OM"])
def test_seqdist_r_norm_uses_dense_compiled_raw_kernel(monkeypatch, method):
    stub = _StubCompiled()
    monkeypatch.setattr(dist, "_backend", "pybind11")
    monkeypatch.setattr(dist, "_pb", stub)
    seqmat = np.array([[0, 1, 2], [0, 2, 2]], dtype=np.int32)
    kwargs = {"indel": 1.0, "subm": 2.0} if method == "OM" else {}

    result = dist.seqdist(seqmat, method=method, norm="YujianBo", **kwargs)

    assert result.shape == (2, 2)
    assert stub.called_om if method == "OM" else stub.called_lcs


@pytest.mark.parametrize("sanitize", [False, True], ids=["release", "ubsan"])
def test_prefix_header_arithmetic_at_int_boundaries(tmp_path, sanitize):
    """Sparse virtual mappings touch only endpoint pages, never huge resident arrays."""
    import shutil
    import subprocess
    import sys
    from pathlib import Path

    compiler = shutil.which("g++")
    if sys.platform != "linux" or compiler is None:
        pytest.skip("bounded native header probe requires the maintained Linux compiler")
    source = tmp_path / "prefix.cpp"
    source.write_text(r"""
#include <sys/mman.h>
#include <climits>
#include <cstdint>
#include "kernels.h"
int main() {
    const int lengths[] = {1, (1 << 30) - 1, 1 << 30, (1 << 30) + 1, INT_MAX};
    for (const int length : lengths) {
        // LCS shares this scalar conversion after its recurrence. Test the
        // arithmetic domain without constructing a huge LCS workspace.
        if (traminex_core::common_length_distance(length, 0) != 2.0 * length) return 3;
        if (traminex_core::common_length_distance(length, length - 1) != 2.0) return 4;
        if (traminex_core::common_length_distance(length, length) != 0.0) return 5;
        const std::size_t size = std::size_t(2) * length * sizeof(std::int32_t);
        auto* data = static_cast<std::int32_t*>(mmap(nullptr, size,
            PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0));
        if (data == MAP_FAILED) return 77;
        data[length] = 1;
        data[std::size_t(2) * length - 1] = 1;
        double result[4];
        traminex_core::pairwise_lcp(data, 2, length, result);
        if (result[1] != 2.0 * length || result[2] != result[1]) return 1;
        traminex_core::pairwise_rlcp(data, 2, length, result);
        if (result[1] != 2.0 * length || result[2] != result[1]) return 2;
        munmap(data, size);
    }
    return 0;
}
""")
    executable = tmp_path / "prefix"
    root = Path(__file__).resolve().parents[1]
    flags = ["-fsanitize=undefined", "-fno-sanitize-recover=undefined"] if sanitize else ["-O2"]
    subprocess.run(
        [
            compiler,
            "-std=c++17",
            *flags,
            "-I",
            str(root / "cpp/src"),
            str(source),
            "-o",
            str(executable),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    result = subprocess.run([str(executable)], capture_output=True, text=True)
    if result.returncode == 77:
        pytest.skip("host refuses bounded sparse virtual mapping")
    assert result.returncode == 0, result.stdout + result.stderr
