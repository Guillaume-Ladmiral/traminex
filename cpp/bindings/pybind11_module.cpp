#include "kernels.h"
#include "accel_kernels.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace py = pybind11;

// The module is private and the dispatcher validates what it passes, but a direct
// call is checked here too (TX-AUD-20261011-045): ValueError for what a kernel
// would misread, and only for what it reads.

// Validate the original dtype/range before pybind11's forcecast can narrow it.
// Private equality-only kernels retain signed int32 labels (including negatives).
static py::array_t<int32_t, py::array::c_style | py::array::forcecast>
narrow_state_codes(const py::array& seq) {
    if (seq.ndim() != 2) throw py::value_error("seq must be 2D [n,L]");
    const char kind = seq.dtype().kind();
    if (kind != 'i' && kind != 'u') {
        throw py::type_error("seq must contain integer state codes");
    }
    if (seq.size() > 0) {
        const auto maximum = seq.attr("max")();
        bool outside;
        if (kind == 'u') {
            outside = maximum.cast<uint64_t>() > static_cast<uint64_t>(std::numeric_limits<int32_t>::max());
        } else {
            outside = maximum.cast<int64_t>() > std::numeric_limits<int32_t>::max()
                || seq.attr("min")().cast<int64_t>() < std::numeric_limits<int32_t>::min();
        }
        if (outside) throw py::value_error("seq state codes must fit a 32-bit integer");
    }
    for (auto dimension : seq.request().shape) {
        if (dimension > std::numeric_limits<int>::max()) {
            throw py::value_error("seq dimensions must fit a 32-bit integer");
        }
    }
    return py::array_t<int32_t, py::array::c_style | py::array::forcecast>::ensure(seq);
}

static void validate_state_codes(
    const py::buffer_info& seq,
    int state_count,
    const char* context
) {
    const auto* values = static_cast<const int32_t*>(seq.ptr);
    for (py::ssize_t i = 0; i < seq.size; ++i) {
        if (values[i] < 0 || values[i] >= state_count) {
            throw py::value_error(
                std::string(context) + " state codes must be in [0, "
                + std::to_string(state_count) + ")"
            );
        }
    }
}

// Zero stays valid: the dispatcher passes it for an all-zero TRATE matrix.
static void validate_indel(double indel) {
    if (!std::isfinite(indel) || indel < 0.0) {
        throw py::value_error("indel must be finite and non-negative");
    }
}

static void validate_scalar_cost(double subm) {
    if (!std::isfinite(subm) || subm < 0.0) {
        throw py::value_error("substitution cost must be finite and non-negative");
    }
}

// The costs between the states that occur, the only ones the kernel reads (as the
// dispatcher's own check); codes are already known to lie in [0, K).
static void validate_read_costs(const py::buffer_info& seq, const py::buffer_info& subm, int K) {
    const auto* codes = static_cast<const int32_t*>(seq.ptr);
    const auto* costs = static_cast<const double*>(subm.ptr);
    std::vector<char> occurs(static_cast<std::size_t>(K), 0);
    for (py::ssize_t i = 0; i < seq.size; ++i) occurs[static_cast<std::size_t>(codes[i])] = 1;
    for (int a = 0; a < K; ++a) {
        if (!occurs[static_cast<std::size_t>(a)]) continue;
        for (int b = 0; b < K; ++b) {
            if (!occurs[static_cast<std::size_t>(b)]) continue;
            const double value = costs[static_cast<std::size_t>(a) * K + b];
            if (!std::isfinite(value) || value < 0.0) {
                throw py::value_error(
                    "substitution costs between occurring states must be finite and non-negative"
                );
            }
        }
    }
}

// DHD reads freq[t, code] for the codes at position t; each cost 1 - (p + q)/2 (or
// 2 - (p + q)) is finite and non-negative when those proportions lie in [0, 1].
static void validate_read_frequencies(const py::buffer_info& seq, const py::buffer_info& freq, int K) {
    const auto* codes = static_cast<const int32_t*>(seq.ptr);
    const auto* values = static_cast<const double*>(freq.ptr);
    const py::ssize_t n = seq.shape[0];
    const py::ssize_t L = seq.shape[1];
    for (py::ssize_t i = 0; i < n; ++i) {
        for (py::ssize_t t = 0; t < L; ++t) {
            const double p = values[static_cast<std::size_t>(t) * K + codes[i * L + t]];
            if (!(p >= 0.0 && p <= 1.0)) {
                throw py::value_error("freq values read by DHD must lie in [0, 1]");
            }
        }
    }
}

static const int32_t* dense_lcs_codes(
    const py::buffer_info& seq,
    std::vector<int32_t>& remapped
) {
    const auto* values = static_cast<const int32_t*>(seq.ptr);
    bool allocation_safe = true;
    for (py::ssize_t i = 0; i < seq.size; ++i) {
        // Keep already-dense labels unchanged at the binding boundary.
        if (values[i] < 0 || static_cast<py::ssize_t>(values[i]) >= seq.size) {
            allocation_safe = false;
            break;
        }
    }
    if (allocation_safe) return values;

    // LCS depends only on equality, so coordinate-compress arbitrary int32
    // labels into a compact, stable representation without rejecting valid
    // sparse or negative codes.
    remapped.reserve(static_cast<std::size_t>(seq.size));
    std::unordered_map<int32_t, int32_t> labels;
    labels.reserve(static_cast<std::size_t>(seq.size));
    for (py::ssize_t i = 0; i < seq.size; ++i) {
        auto entry = labels.emplace(values[i], static_cast<int32_t>(labels.size())).first;
        remapped.push_back(entry->second);
    }
    return remapped.data();
}

static py::array_t<double> pb_pairwise_ham_unit(const py::array& seq) {
    auto dense = narrow_state_codes(seq);
    auto buf = dense.request();
    int n = static_cast<int>(buf.shape[0]);
    int L = static_cast<int>(buf.shape[1]);
    auto out = py::array_t<double>({n, n});
    const auto* codes = static_cast<const int32_t*>(buf.ptr);
    auto* result = static_cast<double*>(out.request().ptr);
    {
        py::gil_scoped_release release;
        traminex_accel::pairwise_ham_unit(codes, n, L, result);
    }
    return out;
}

static py::array_t<double> pb_pairwise_om(const py::array& seq,
                                          double indel, double subm) {
    auto dense = narrow_state_codes(seq);
    validate_indel(indel);
    validate_scalar_cost(subm);
    auto buf = dense.request();
    int n = static_cast<int>(buf.shape[0]);
    int L = static_cast<int>(buf.shape[1]);
    auto out = py::array_t<double>({n, n});
    const auto* codes = static_cast<const int32_t*>(buf.ptr);
    auto* result = static_cast<double*>(out.request().ptr);
    // Direct one-row OpenMP OM with repeated-addition borders, byte-identical
    // (atol=0) to the pure-Python / numba OM recurrence.
    {
        py::gil_scoped_release release;
        traminex_accel::pairwise_om(codes, n, L, indel, subm, result);
    }
    return out;
}

static py::array_t<double> pb_pairwise_om_mat(
    const py::array& seq,
    double indel,
    const py::array_t<double, py::array::c_style | py::array::forcecast>& subm
) {
    auto dense = narrow_state_codes(seq);
    validate_indel(indel);
    auto s = dense.request();
    auto m = subm.request();
    if (m.ndim != 2 || m.shape[0] != m.shape[1]) throw py::value_error("subm must be square [K,K]");
    int n = static_cast<int>(s.shape[0]);
    int L = static_cast<int>(s.shape[1]);
    int K = static_cast<int>(m.shape[0]);
    validate_state_codes(s, K, "OM");
    validate_read_costs(s, m, K);
    auto out = py::array_t<double>({n, n});
    const auto* codes = static_cast<const int32_t*>(s.ptr);
    const auto* costs = static_cast<const double*>(m.ptr);
    auto* result = static_cast<double*>(out.request().ptr);
    // RP9: OpenMP OM with matrix substitution costs (accel), byte-identical to Python/numba.
    {
        py::gil_scoped_release release;
        traminex_accel::pairwise_om_mat(codes, n, L, indel, costs, K, result);
    }
    return out;
}

static py::array_t<double> pb_pairwise_rlcp(
    const py::array& seq
) {
    auto dense = narrow_state_codes(seq);
    auto buf = dense.request();
    int n = static_cast<int>(buf.shape[0]);
    int L = static_cast<int>(buf.shape[1]);
    auto out = py::array_t<double>({n, n});
    const auto* codes = static_cast<const int32_t*>(buf.ptr);
    auto* result = static_cast<double*>(out.request().ptr);
    {
        py::gil_scoped_release release;
        traminex_core::pairwise_rlcp(codes, n, L, result);
    }
    return out;
}

static py::array_t<double> pb_pairwise_dhd(
    const py::array& seq,
    const py::array_t<double, py::array::c_style | py::array::forcecast>& freq,
    int variant
) {
    auto dense = narrow_state_codes(seq);
    auto s = dense.request();
    auto f = freq.request();
    if (f.ndim != 2) throw py::value_error("freq must be 2D [L,K]");
    int n = static_cast<int>(s.shape[0]);
    int L = static_cast<int>(s.shape[1]);
    if (static_cast<int>(f.shape[0]) != L) throw py::value_error("freq must have shape [L,K]");
    int K = static_cast<int>(f.shape[1]);
    validate_state_codes(s, K, "DHD");
    validate_read_frequencies(s, f, K);
    if (variant != 0 && variant != 1) {
        throw py::value_error("DHD variant must be 0 (mean) or 1 (sum)");
    }
    auto out = py::array_t<double>({n, n});
    const auto* codes = static_cast<const int32_t*>(s.ptr);
    const auto* frequencies = static_cast<const double*>(f.ptr);
    auto* result = static_cast<double*>(out.request().ptr);
    {
        py::gil_scoped_release release;
        traminex_core::pairwise_dhd(codes, n, L, frequencies, K, variant, result);
    }
    return out;
}

PYBIND11_MODULE(_pb_kernels, m) {
    m.doc() = "pybind11-bound fast sequence distance kernels";
    // Internal diagnostics for the accel (OpenMP) kernels; not part of the public traminex API.
    m.def("accel_max_threads", &traminex_accel::max_threads,
          "OpenMP thread count the accel kernels will use (1 if built without OpenMP)");
    m.def("accel_openmp_enabled", &traminex_accel::openmp_enabled,
          "True if the compiled kernel was built with OpenMP");
    m.def("pairwise_ham_unit", &pb_pairwise_ham_unit, py::arg("seq"),
          "Pairwise raw unit-cost Hamming (mismatch count; R-faithful HAM default, OpenMP)");
    m.def("pairwise_om", &pb_pairwise_om, py::arg("seq"), py::arg("indel")=1.0, py::arg("subm")=2.0,
          "Pairwise Optimal Matching (uniform substitution cost)");
    m.def("pairwise_om_mat", &pb_pairwise_om_mat, py::arg("seq"), py::arg("indel")=1.0, py::arg("subm"),
          "Pairwise Optimal Matching with substitution matrix [K,K]");
    // Optional: LCS/LCP when same-length sequences
    m.def("pairwise_lcs", [](const py::array& seq){
        auto dense = narrow_state_codes(seq);
        auto buf = dense.request();
        int n = static_cast<int>(buf.shape[0]);
        int L = static_cast<int>(buf.shape[1]);
        std::vector<int32_t> remapped;
        const int32_t* values = dense_lcs_codes(buf, remapped);
        auto out = py::array_t<double>({n, n});
        auto* result = static_cast<double*>(out.request().ptr);
        // Direct one-row OpenMP LCS, byte-identical to the Python recurrence.
        {
            py::gil_scoped_release release;
            traminex_accel::pairwise_lcs(values, n, L, result);
        }
        return out;
    });
    m.def("pairwise_lcp", [](const py::array& seq){
        auto dense = narrow_state_codes(seq);
        auto buf = dense.request();
        int n = static_cast<int>(buf.shape[0]);
        int L = static_cast<int>(buf.shape[1]);
        auto out = py::array_t<double>({n, n});
        const auto* codes = static_cast<const int32_t*>(buf.ptr);
        auto* result = static_cast<double*>(out.request().ptr);
        {
            py::gil_scoped_release release;
            traminex_core::pairwise_lcp(codes, n, L, result);
        }
        return out;
    });
    m.def("pairwise_rlcp", &pb_pairwise_rlcp, py::arg("seq"),
          "Pairwise RLCP distance (reverse LCP)");
    m.def("pairwise_dhd", &pb_pairwise_dhd, py::arg("seq"), py::arg("freq"), py::arg("variant") = 0,
          "Pairwise DHD with per-position frequencies [L,K] and variant code (0=mean,1=sum)");
}
