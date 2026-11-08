#pragma once

#include "kernels.h"

// Direct HAM/LCS/OM kernels derived from TraMineR 2.2-14's seqdist contract
// and OMdistance recurrence.  The upper triangle is divided into contiguous
// pair ranges; each worker owns its dynamic-programming storage and writes both
// symmetric result cells for every pair in its range.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace traminex_accel {

inline int max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

inline bool openmp_enabled() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}

namespace detail {

struct NoWorkspace {};

struct LcsWorkspace {
    explicit LcsWorkspace(std::size_t width)
        : row(width + 1U, 0) {}

    std::vector<int> row;
};

struct OmWorkspace {
    explicit OmWorkspace(std::size_t width)
        : row(width + 1U, 0.0) {}

    std::vector<double> row;
};

inline std::uint64_t row_hash(const std::int32_t* row, int width) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (int position = 0; position < width; ++position) {
        hash ^= static_cast<std::uint32_t>(row[position]);
        hash *= 1099511628211ULL;
    }
    return hash;
}

template <typename WorkspaceFactory, typename Distance>
inline void fill_symmetric(
    const std::int32_t* sequences,
    int count,
    int width,
    double* output,
    WorkspaceFactory make_workspace,
    Distance distance
) {
    using Workspace = decltype(make_workspace());
    if (count <= 0) return;
    std::fill(
        output,
        output + static_cast<std::size_t>(count) * static_cast<std::size_t>(count),
        0.0
    );
    if (count == 1 || width == 0) return;

    std::vector<std::vector<int>> groups;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> hash_groups;
    hash_groups.reserve(static_cast<std::size_t>(count));
    for (int index = 0; index < count; ++index) {
        const std::int32_t* candidate = sequences
            + static_cast<std::size_t>(index) * static_cast<std::size_t>(width);
        const std::uint64_t hash = row_hash(candidate, width);
        auto& matching_groups = hash_groups[hash];
        std::size_t group_index = groups.size();
        for (const std::size_t possible : matching_groups) {
            const std::int32_t* representative = sequences
                + static_cast<std::size_t>(groups[possible][0])
                    * static_cast<std::size_t>(width);
            if (std::equal(candidate, candidate + width, representative)) {
                group_index = possible;
                break;
            }
        }
        if (group_index == groups.size()) {
            groups.push_back({index});
            matching_groups.push_back(group_index);
        } else {
            groups[group_index].push_back(index);
        }
    }

    const int group_count = static_cast<int>(groups.size());
    if (group_count == 1) return;

    const int worker_count = max_threads();
    std::vector<Workspace> workspaces;
    workspaces.reserve(static_cast<std::size_t>(worker_count));
    for (int worker = 0; worker < worker_count; ++worker) {
        workspaces.emplace_back(make_workspace());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
        const int worker = omp_get_thread_num();
        const int active_workers = omp_get_num_threads();
#else
        const int worker = 0;
        const int active_workers = 1;
#endif
        Workspace& workspace = workspaces[static_cast<std::size_t>(worker)];
        const std::size_t total_pairs = static_cast<std::size_t>(group_count)
            * static_cast<std::size_t>(group_count - 1) / 2U;
        const std::size_t pair_begin = total_pairs * static_cast<std::size_t>(worker)
            / static_cast<std::size_t>(active_workers);
        const std::size_t pair_end = total_pairs * static_cast<std::size_t>(worker + 1)
            / static_cast<std::size_t>(active_workers);

        int left_index = 0;
        std::size_t row_start = 0U;
        std::size_t row_pairs = static_cast<std::size_t>(group_count - 1);
        while (pair_begin >= row_start + row_pairs && row_pairs > 0U) {
            row_start += row_pairs;
            --row_pairs;
            ++left_index;
        }
        int right_index = left_index + 1
            + static_cast<int>(pair_begin - row_start);

        for (std::size_t pair = pair_begin; pair < pair_end; ++pair) {
            const auto& left_group = groups[static_cast<std::size_t>(left_index)];
            const auto& right_group = groups[static_cast<std::size_t>(right_index)];
            const std::int32_t* left = sequences
                + static_cast<std::size_t>(left_group[0]) * static_cast<std::size_t>(width);
            const std::int32_t* right = sequences
                + static_cast<std::size_t>(right_group[0]) * static_cast<std::size_t>(width);
            const double value = distance(left, right, workspace);
            for (const int left_member : left_group) {
                for (const int right_member : right_group) {
                    const std::size_t left_row = static_cast<std::size_t>(left_member);
                    const std::size_t right_row = static_cast<std::size_t>(right_member);
                    output[left_row * static_cast<std::size_t>(count) + right_row] = value;
                    output[right_row * static_cast<std::size_t>(count) + left_row] = value;
                }
            }

            ++right_index;
            if (right_index == group_count) {
                ++left_index;
                right_index = left_index + 1;
            }
        }
    }
}

inline double hamming(
    const std::int32_t* left,
    const std::int32_t* right,
    int width,
    NoWorkspace&
) {
    int mismatches = 0;
    for (int position = 0; position < width; ++position) {
        mismatches += left[position] != right[position];
    }
    return static_cast<double>(mismatches);
}

inline double lcs(
    const std::int32_t* left,
    const std::int32_t* right,
    int width,
    LcsWorkspace& workspace
) {
    std::fill(workspace.row.begin(), workspace.row.end(), 0);
    for (int left_position = 0; left_position < width; ++left_position) {
        int diagonal = 0;
        for (int column = 0; column < width; ++column) {
            const std::size_t cell = static_cast<std::size_t>(column) + 1U;
            const int above = workspace.row[cell];
            if (left[left_position] == right[column]) {
                workspace.row[cell] = diagonal + 1;
            } else {
                workspace.row[cell] = std::max(above, workspace.row[cell - 1U]);
            }
            diagonal = above;
        }
    }
    const int common = workspace.row[static_cast<std::size_t>(width)];
    return traminex_core::common_length_distance(width, common);
}

template <typename SubstitutionCost>
inline double optimal_matching(
    const std::int32_t* left,
    const std::int32_t* right,
    int width,
    double indel,
    SubstitutionCost substitution_cost,
    OmWorkspace& workspace
) {
    const std::size_t columns = static_cast<std::size_t>(width);
    workspace.row[0] = 0.0;
    // width may equal INT_MAX: an inclusive signed counter would overflow
    // when incremented after its final iteration.
    for (std::size_t column = 1U; column <= columns; ++column) {
        workspace.row[column] = workspace.row[column - 1U] + indel;
    }
    for (std::size_t row = 1U; row <= columns; ++row) {
        double diagonal = workspace.row[0];
        workspace.row[0] += indel;
        for (std::size_t column = 1U; column <= columns; ++column) {
            const std::size_t cell = column;
            const double above = workspace.row[cell];
            const double deletion =
                above + indel;
            const double insertion =
                workspace.row[cell - 1U] + indel;
            const double substitution =
                diagonal
                + substitution_cost(left[row - 1], right[column - 1]);
            workspace.row[cell] =
                std::min(std::min(deletion, insertion), substitution);
            diagonal = above;
        }
    }
    return workspace.row[static_cast<std::size_t>(width)];
}

}  // namespace detail

inline void pairwise_ham_unit(
    const std::int32_t* sequences,
    int count,
    int width,
    double* output
) {
    detail::fill_symmetric(
        sequences,
        count,
        width,
        output,
        [] { return detail::NoWorkspace{}; },
        [width](const std::int32_t* left, const std::int32_t* right, detail::NoWorkspace& ws) {
            return detail::hamming(left, right, width, ws);
        }
    );
}

inline void pairwise_lcs(
    const std::int32_t* sequences,
    int count,
    int width,
    double* output
) {
    detail::fill_symmetric(
        sequences,
        count,
        width,
        output,
        [width] { return detail::LcsWorkspace(static_cast<std::size_t>(width)); },
        [width](const std::int32_t* left, const std::int32_t* right, detail::LcsWorkspace& ws) {
            return detail::lcs(left, right, width, ws);
        }
    );
}

inline void pairwise_om(
    const std::int32_t* sequences,
    int count,
    int width,
    double indel,
    double substitution,
    double* output
) {
    detail::fill_symmetric(
        sequences,
        count,
        width,
        output,
        [width] { return detail::OmWorkspace(static_cast<std::size_t>(width)); },
        [width, indel, substitution](
            const std::int32_t* left,
            const std::int32_t* right,
            detail::OmWorkspace& ws
        ) {
            const auto cost = [substitution](std::int32_t a, std::int32_t b) {
                return a == b ? 0.0 : substitution;
            };
            return detail::optimal_matching(left, right, width, indel, cost, ws);
        }
    );
}

inline void pairwise_om_mat(
    const std::int32_t* sequences,
    int count,
    int width,
    double indel,
    const double* substitution,
    int state_count,
    double* output
) {
    detail::fill_symmetric(
        sequences,
        count,
        width,
        output,
        [width] { return detail::OmWorkspace(static_cast<std::size_t>(width)); },
        [width, indel, substitution, state_count](
            const std::int32_t* left,
            const std::int32_t* right,
            detail::OmWorkspace& ws
        ) {
            const auto cost = [substitution, state_count](std::int32_t a, std::int32_t b) {
                if (a == b) return 0.0;
                const std::size_t index = static_cast<std::size_t>(a)
                    * static_cast<std::size_t>(state_count)
                    + static_cast<std::size_t>(b);
                return substitution[index];
            };
            return detail::optimal_matching(left, right, width, indel, cost, ws);
        }
    );
}

}  // namespace traminex_accel
