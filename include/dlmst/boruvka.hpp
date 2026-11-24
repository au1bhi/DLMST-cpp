#pragma once

#include "dlmst/types.hpp"

#include <utility>
#include <vector>

namespace dlmst {

struct BoruvkaResult {
    std::vector<int> contracted_orig_ids;
    int contracted_n = 0;
    std::vector<Edge> contracted_edges;
};

/**
 * Runs up to `max_phases` Boruvka contraction steps (Definition 2.8 / Algorithm 2 line 1),
 * discarding self-loops and isolated vertices after contraction.
 */
inline BoruvkaResult run_boruvka_phases(
    int n,
    const std::vector<Edge>& edges,
    int max_phases,
    WeightOracle& oracle) {
    BoruvkaResult res;
    int cur_n = n;
    std::vector<Edge> cur_edges;
    cur_edges.reserve(edges.size());
    for (const auto& e : edges) {
        if (e.u != e.v) cur_edges.push_back(e);
    }

    for (int phase = 0; phase < max_phases; ++phase) {
        if (cur_n <= 1 || cur_edges.empty()) break;

        std::vector<int> best_idx(cur_n, -1);
        for (int i = 0; i < static_cast<int>(cur_edges.size()); ++i) {
            const auto& e = cur_edges[i];
            if (best_idx[e.u] == -1 ||
                oracle.less(e.orig_id, cur_edges[best_idx[e.u]].orig_id, CmpCategory::Boruvka)) {
                best_idx[e.u] = i;
            }
            if (best_idx[e.v] == -1 ||
                oracle.less(e.orig_id, cur_edges[best_idx[e.v]].orig_id, CmpCategory::Boruvka)) {
                best_idx[e.v] = i;
            }
        }

        DSU dsu(cur_n);
        std::vector<bool> picked(cur_edges.size(), false);
        bool any_picked = false;
        for (int v = 0; v < cur_n; ++v) {
            const int idx = best_idx[v];
            if (idx != -1 && !picked[idx]) {
                picked[idx] = true;
                any_picked = true;
                res.contracted_orig_ids.push_back(cur_edges[idx].orig_id);
                dsu.unite(cur_edges[idx].u, cur_edges[idx].v);
            }
        }
        if (!any_picked) break;

        std::vector<int> comp(cur_n);
        for (int v = 0; v < cur_n; ++v) {
            comp[v] = dsu.find(v);
        }

        std::vector<int> comp_remap(cur_n, -1);
        int new_n = 0;
        std::vector<Edge> next_edges;
        next_edges.reserve(cur_edges.size());
        for (const auto& e : cur_edges) {
            const int cu = comp[e.u];
            const int cv = comp[e.v];
            if (cu == cv) continue;
            if (comp_remap[cu] == -1) comp_remap[cu] = new_n++;
            if (comp_remap[cv] == -1) comp_remap[cv] = new_n++;
            next_edges.push_back({comp_remap[cu], comp_remap[cv], e.orig_id});
        }

        cur_n = new_n;
        cur_edges = std::move(next_edges);
    }

    res.contracted_n = cur_n;
    res.contracted_edges = std::move(cur_edges);
    return res;
}

}  // namespace dlmst
