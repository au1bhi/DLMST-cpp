#pragma once

#include "dlmst/types.hpp"

#include <bit>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

namespace dlmst {

struct ForestEdge {
    int u = -1;
    int v = -1;
    int edge_id = -1;
};

namespace detail {

/**
 * Hagerup (2009) bitmask filter:
 * Keeps exactly those bits k in S that serve as the minimum set bit >= t
 * for at least one query top-depth t in D.
 * Bits of S and D must lie in [1, 61].
 */
inline uint64_t hagerup_filter_mask(uint64_t S, uint64_t D) {
    return (((~S + (D & ~S)) & S) | (D & S));
}

}  // namespace detail

/**
 * Deterministic O(n + q) path-maximum queries on a forest using:
 *   - King (1997) Boruvka tree reduction (depth <= log2(n) < 60),
 *   - Tarjan's offline LCA on the Boruvka tree (0 comparisons),
 *   - Komlos (1985) / Hagerup (2009) bitmask + binary search traversal.
 *
 * `cmp(edge_id_a, edge_id_b)` must return true iff edge_id_a is strictly lighter
 * than edge_id_b under a strict total order.
 *
 * Returns a vector of size `queries.size()`, where entry i is the edge_id of the
 * maximum-weight edge on the unique forest path between queries[i].first and
 * queries[i].second, or -1 if the endpoints are equal or disconnected in `forest`.
 */
template <typename Cmp>
std::vector<int> path_max_queries(
    int n,
    const std::vector<ForestEdge>& forest,
    const std::vector<std::pair<int, int>>& queries,
    Cmp&& cmp) {
    std::vector<int> answers(queries.size(), -1);
    if (n <= 1 || forest.empty() || queries.empty()) {
        return answers;
    }

    // 1. Identify connected components of `forest` via BFS (0 comparisons).
    std::vector<std::vector<std::pair<int, int>>> adj(n);
    for (int i = 0; i < static_cast<int>(forest.size()); ++i) {
        const auto& fe = forest[i];
        if (fe.u == fe.v) continue;
        adj[fe.u].push_back({fe.v, fe.edge_id});
        adj[fe.v].push_back({fe.u, fe.edge_id});
    }

    std::vector<int> comp_id(n, -1);
    std::vector<std::vector<int>> comp_verts;
    std::vector<int> queue;
    for (int s = 0; s < n; ++s) {
        if (comp_id[s] != -1 || adj[s].empty()) continue;
        const int cid = static_cast<int>(comp_verts.size());
        comp_verts.emplace_back();
        queue.clear();
        queue.push_back(s);
        comp_id[s] = cid;
        for (size_t head = 0; head < queue.size(); ++head) {
            int u = queue[head];
            comp_verts[cid].push_back(u);
            for (const auto& [v, eid] : adj[u]) {
                (void)eid;
                if (comp_id[v] == -1) {
                    comp_id[v] = cid;
                    queue.push_back(v);
                }
            }
        }
    }

    // 2. Group queries by connected component.
    std::vector<std::vector<int>> comp_queries(comp_verts.size());
    for (int qi = 0; qi < static_cast<int>(queries.size()); ++qi) {
        const int u = queries[qi].first;
        const int v = queries[qi].second;
        if (u < 0 || u >= n || v < 0 || v >= n || u == v) continue;
        if (comp_id[u] != -1 && comp_id[u] == comp_id[v]) {
            comp_queries[comp_id[u]].push_back(qi);
        }
    }

    std::vector<int> local_v(n, -1);

    // 3. Process each component with >= 1 query independently.
    for (size_t cid = 0; cid < comp_verts.size(); ++cid) {
        if (comp_queries[cid].empty()) continue;

        const auto& verts = comp_verts[cid];
        const int sz = static_cast<int>(verts.size());
        assert(sz >= 2);

        for (int i = 0; i < sz; ++i) {
            local_v[verts[i]] = i;
        }

        std::vector<ForestEdge> cur_edges;
        cur_edges.reserve(sz - 1);
        for (int u_glob : verts) {
            for (const auto& [v_glob, eid] : adj[u_glob]) {
                if (u_glob < v_glob) {
                    cur_edges.push_back({local_v[u_glob], local_v[v_glob], eid});
                }
            }
        }

        // Build King's Boruvka tree B. Leaves are 0..sz-1.
        std::vector<int> cur_nodes(sz);
        for (int i = 0; i < sz; ++i) cur_nodes[i] = i;

        // Upper bound on B nodes is 2 * sz.
        std::vector<std::vector<int>> b_children(sz);
        std::vector<int> b_edge(sz, -1);
        std::vector<int> node_depth(sz, 0);

        std::vector<int> node_to_pos(sz, -1);
        int rounds = 0;

        while (cur_nodes.size() > 1) {
            const int k = static_cast<int>(cur_nodes.size());
            if (static_cast<int>(node_to_pos.size()) < static_cast<int>(b_children.size())) {
                node_to_pos.resize(b_children.size(), -1);
            }
            for (int i = 0; i < k; ++i) {
                node_to_pos[cur_nodes[i]] = i;
            }

            // Find minimum incident edge for each active Boruvka node.
            std::vector<int> best_e_idx(k, -1);
            for (int ei = 0; ei < static_cast<int>(cur_edges.size()); ++ei) {
                const auto& e = cur_edges[ei];
                const int pu = node_to_pos[e.u];
                const int pv = node_to_pos[e.v];
                if (best_e_idx[pu] == -1 || cmp(e.edge_id, cur_edges[best_e_idx[pu]].edge_id)) {
                    best_e_idx[pu] = ei;
                }
                if (best_e_idx[pv] == -1 || cmp(e.edge_id, cur_edges[best_e_idx[pv]].edge_id)) {
                    best_e_idx[pv] = ei;
                }
            }

            // Group nodes connected by selected edges using DSU.
            DSU dsu(k);
            for (int i = 0; i < k; ++i) {
                const int ei = best_e_idx[i];
                assert(ei != -1);
                b_edge[cur_nodes[i]] = cur_edges[ei].edge_id;
                const int pu = node_to_pos[cur_edges[ei].u];
                const int pv = node_to_pos[cur_edges[ei].v];
                dsu.unite(pu, pv);
            }

            std::vector<int> rep_to_parent(k, -1);
            std::vector<int> next_nodes;
            next_nodes.reserve(k / 2 + 1);

            for (int i = 0; i < k; ++i) {
                const int r = dsu.find(i);
                if (rep_to_parent[r] == -1) {
                    const int pid = static_cast<int>(b_children.size());
                    rep_to_parent[r] = pid;
                    b_children.emplace_back();
                    b_edge.push_back(-1);
                    node_depth.push_back(0);
                    next_nodes.push_back(pid);
                }
                b_children[rep_to_parent[r]].push_back(cur_nodes[i]);
            }

            std::vector<ForestEdge> next_edges;
            next_edges.reserve(cur_edges.size());
            for (const auto& e : cur_edges) {
                const int pu = node_to_pos[e.u];
                const int pv = node_to_pos[e.v];
                const int p1 = rep_to_parent[dsu.find(pu)];
                const int p2 = rep_to_parent[dsu.find(pv)];
                if (p1 != p2) {
                    next_edges.push_back({p1, p2, e.edge_id});
                }
            }

            cur_nodes = std::move(next_nodes);
            cur_edges = std::move(next_edges);
            ++rounds;
        }

        const int total_b = static_cast<int>(b_children.size());
        const int root = cur_nodes[0];
        const int h = rounds;
        assert(h >= 1 && h <= 61);

        // Assign depths top-down in B (0 comparisons).
        {
            std::vector<int> stack = {root};
            node_depth[root] = 0;
            while (!stack.empty()) {
                int u = stack.back();
                stack.pop_back();
                for (int ch : b_children[u]) {
                    node_depth[ch] = node_depth[u] + 1;
                    stack.push_back(ch);
                }
            }
        }

        // Tarjan's offline LCA on B (0 comparisons).
        const auto& cqueries = comp_queries[cid];
        const int num_cq = static_cast<int>(cqueries.size());

        struct LeafQueryRef {
            int cq_idx;
            int other_leaf;
        };
        std::vector<std::vector<LeafQueryRef>> leaf_q(sz);
        for (int i = 0; i < num_cq; ++i) {
            const int qi = cqueries[i];
            const int lu = local_v[queries[qi].first];
            const int lv = local_v[queries[qi].second];
            leaf_q[lu].push_back({i, lv});
            leaf_q[lv].push_back({i, lu});
        }

        std::vector<int> lca_node(num_cq, -1);
        {
            DSU lca_dsu(total_b);
            std::vector<int> ancestor(total_b);
            for (int i = 0; i < total_b; ++i) ancestor[i] = i;
            std::vector<bool> visited(total_b, false);

            auto dfs_lca = [&](auto& self, int u) -> void {
                ancestor[u] = u;
                for (int ch : b_children[u]) {
                    self(self, ch);
                    lca_dsu.unite(u, ch);
                    ancestor[lca_dsu.find(u)] = u;
                }
                visited[u] = true;
                if (u < sz) {
                    for (const auto& ref : leaf_q[u]) {
                        if (visited[ref.other_leaf]) {
                            lca_node[ref.cq_idx] = ancestor[lca_dsu.find(ref.other_leaf)];
                        }
                    }
                }
            };
            dfs_lca(dfs_lca, root);
        }

        // Attach vertical sub-queries to leaves and compute D_mask bottom-up.
        struct VerticalQuery {
            int cq_idx;
            int side;
            int top_depth;
        };
        std::vector<std::vector<VerticalQuery>> leaf_vqueries(sz);
        std::vector<uint64_t> d_mask(total_b, 0ULL);

        for (int i = 0; i < num_cq; ++i) {
            const int qi = cqueries[i];
            const int lu = local_v[queries[qi].first];
            const int lv = local_v[queries[qi].second];
            const int lca = lca_node[i];
            assert(lca >= sz);
            const int top_d = node_depth[lca] + 1;
            assert(top_d >= 1 && top_d <= h);
            leaf_vqueries[lu].push_back({i, 0, top_d});
            leaf_vqueries[lv].push_back({i, 1, top_d});
            const uint64_t bit = 1ULL << top_d;
            d_mask[lu] |= bit;
            d_mask[lv] |= bit;
        }

        // Internal nodes were created in strictly increasing depth-from-bottom order (sz..total_b-1).
        for (int u = sz; u < total_b; ++u) {
            uint64_t m = 0ULL;
            for (int ch : b_children[u]) {
                m |= d_mask[ch];
            }
            d_mask[u] = m;
        }

        // Top-down DFS with Hagerup bitmask filter + binary search.
        std::vector<int> vert_ans_0(num_cq, -1);
        std::vector<int> vert_ans_1(num_cq, -1);
        int path_edge[64] = {};

        auto dfs_solve = [&](auto& self, int p, uint64_t S_p) -> void {
            for (int v : b_children[p]) {
                if (d_mask[v] == 0ULL) continue;
                const int d = node_depth[v];
                path_edge[d] = b_edge[v];

                const uint64_t D_above = d_mask[v] & ((1ULL << d) - 1ULL);
                uint64_t S_cur = detail::hagerup_filter_mask(S_p, D_above);

                int bits[64];
                int len = 0;
                for (uint64_t tmp = S_cur; tmp != 0ULL; tmp &= (tmp - 1ULL)) {
                    bits[len++] = std::countr_zero(tmp);
                }

                // Binary search for the shallowest (heaviest) suffix-max depth in S_cur
                // that is strictly lighter than path_edge[d].
                int lo = 0;
                int hi = len;
                while (lo < hi) {
                    const int mid = lo + (hi - lo) / 2;
                    if (cmp(path_edge[bits[mid]], path_edge[d])) {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                if (lo < len) {
                    S_cur &= (1ULL << bits[lo]) - 1ULL;
                }

                uint64_t S_v = S_cur | (1ULL << d);
                const uint64_t D_upto_d =
                    (d == 63) ? d_mask[v] : (d_mask[v] & ((1ULL << (d + 1)) - 1ULL));
                S_v = detail::hagerup_filter_mask(S_v, D_upto_d);

                if (v < sz) {
                    for (const auto& vq : leaf_vqueries[v]) {
                        const uint64_t active = S_v & (~0ULL << vq.top_depth);
                        assert(active != 0ULL);
                        const int best_d = std::countr_zero(active);
                        if (vq.side == 0) {
                            vert_ans_0[vq.cq_idx] = path_edge[best_d];
                        } else {
                            vert_ans_1[vq.cq_idx] = path_edge[best_d];
                        }
                    }
                } else {
                    self(self, v, S_v);
                }
            }
        };

        dfs_solve(dfs_solve, root, 0ULL);

        // Combine the two vertical halves for each query (at most 1 comparison).
        for (int i = 0; i < num_cq; ++i) {
            const int a = vert_ans_0[i];
            const int b = vert_ans_1[i];
            assert(a != -1 && b != -1);
            if (a == b) {
                answers[cqueries[i]] = a;
            } else {
                answers[cqueries[i]] = cmp(a, b) ? b : a;
            }
        }

        for (int u_glob : verts) {
            local_v[u_glob] = -1;
        }
    }

    return answers;
}

}  // namespace dlmst
