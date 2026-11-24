#pragma once

#include "dlmst/types.hpp"

#include <algorithm>
#include <numeric>
#include <vector>

namespace dlmst {

/**
 * Ground-truth Kruskal MSF baseline using the exact same tie-breaking rule
 * (weight, orig_id) and optional comparison counting.
 */
inline std::vector<int> kruskal_msf(const Graph& g, uint64_t* cmp_count = nullptr) {
    std::vector<int> order;
    order.reserve(g.edges.size());
    for (int i = 0; i < static_cast<int>(g.edges.size()); ++i) {
        if (g.edges[i].u != g.edges[i].v) {
            order.push_back(i);
        }
    }

    uint64_t cmps = 0;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        ++cmps;
        const int ea = g.edges[a].orig_id;
        const int eb = g.edges[b].orig_id;
        const double wa = g.weights[ea];
        const double wb = g.weights[eb];
        if (wa != wb) return wa < wb;
        return ea < eb;
    });

    if (cmp_count) *cmp_count = cmps;

    DSU dsu(g.n);
    std::vector<int> result;
    result.reserve(g.n > 0 ? g.n - 1 : 0);
    for (int idx : order) {
        const Edge& e = g.edges[idx];
        if (dsu.unite(e.u, e.v)) {
            result.push_back(e.orig_id);
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

}  // namespace dlmst
