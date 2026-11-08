#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>

namespace traminex_core {

inline double common_length_distance(int length, int common) {
    return 2.0 * static_cast<double>(length) - 2.0 * static_cast<double>(common);
}

inline void pairwise_lcp(const int32_t* seq, int n, int L, double* out) {
    std::fill(out, out + static_cast<size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {
        const int32_t* ai = seq + static_cast<size_t>(i) * L;
        for (int j = i + 1; j < n; ++j) {
            const int32_t* bj = seq + static_cast<size_t>(j) * L;
            int k = 0;
            while (k < L && ai[k] == bj[k]) ++k;
            double d = common_length_distance(L, k);
            out[static_cast<size_t>(i) * n + j] = d;
            out[static_cast<size_t>(j) * n + i] = d;
        }
    }
}

inline void pairwise_rlcp(const int32_t* seq, int n, int L, double* out) {
    std::fill(out, out + static_cast<size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {
        const int32_t* ai = seq + static_cast<size_t>(i) * L;
        for (int j = i + 1; j < n; ++j) {
            const int32_t* bj = seq + static_cast<size_t>(j) * L;
            int k = 0;
            while (k < L && ai[L - 1 - k] == bj[L - 1 - k]) ++k;
            double d = common_length_distance(L, k);
            out[static_cast<size_t>(i) * n + j] = d;
            out[static_cast<size_t>(j) * n + i] = d;
        }
    }
}

inline double dhd_distance(const int32_t* a, const int32_t* b, int L, const double* freq, int K, int variant) {
    // freq is row-major [L, K], variant: 0="mean", 1="sum"
    double acc = 0.0;
    for (int t = 0; t < L; ++t) {
        int ai = a[t];
        int bj = b[t];
        if (ai == bj) continue;
        double p_ai = freq[static_cast<size_t>(t) * K + ai];
        double p_bj = freq[static_cast<size_t>(t) * K + bj];
        if (variant == 0) acc += 1.0 - 0.5 * (p_ai + p_bj);
        else acc += 2.0 - (p_ai + p_bj);
    }
    return (L > 0) ? (acc / static_cast<double>(L)) : 0.0;
}

inline void pairwise_dhd(const int32_t* seq, int n, int L, const double* freq, int K, int variant, double* out) {
    std::fill(out, out + static_cast<size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {
        const int32_t* ai = seq + static_cast<size_t>(i) * L;
        for (int j = i + 1; j < n; ++j) {
            const int32_t* bj = seq + static_cast<size_t>(j) * L;
            double d = dhd_distance(ai, bj, L, freq, K, variant);
            out[static_cast<size_t>(i) * n + j] = d;
            out[static_cast<size_t>(j) * n + i] = d;
        }
    }
}

} // namespace traminex_core
