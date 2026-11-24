#pragma once

#include "dlmst/path_max.hpp"
#include "dlmst/types.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

namespace dlmst {

enum class EdgeAction : uint8_t {
    Find = 0,
    Project = 1
};

/**
 * Estimates preferred actions for edges in A_i (Section 3.5 / Algorithm 1 line 10).
 *
 * True rule:
 *   - Find if sigma[u] == -1 || sigma[v] == -1 || pi[u] == pi[v]
 *   - Project if u, v are disconnected in F
 *   - Otherwise Find iff F_{<= p(e)} connects u and v, where p(e) = max(w(sigma[u]), w(sigma[v])).
 *
 * In DAGTopological mode, queries a deterministic O(t * sqrt(m_i * n)) sample S*
 * on F, records every comparison edge into a DAG on ground set T = F U {sigma[v]},
 * topologically sorts T (0 comparisons) into a consistent total order, and resolves
 * all remaining edges in A_i using non-counting path_max_queries over topo_rank.
 */
inline std::vector<EdgeAction> estimate_preferred_actions(
    int n,
    const std::vector<CertEdge>& edges,
    const std::vector<int>& F,
    const std::vector<int>& sigma,
    const std::vector<int>& pi,
    const std::vector<int>& A_i,
    int t,
    uint64_t round_seed,
    const DLMSTConfig& config,
    WeightOracle& oracle,
    DLMSTStats* stats = nullptr) {
    std::vector<EdgeAction> actions(A_i.size(), EdgeAction::Find);
    if (A_i.empty()) return actions;

    // 1. Identify connected components of F (0 weight comparisons).
    DSU f_dsu(n);
    for (int f_idx : F) {
        f_dsu.unite(edges[f_idx].u, edges[f_idx].v);
    }

    std::vector<int> ambig_pos;
    ambig_pos.reserve(A_i.size());
    for (int pos = 0; pos < static_cast<int>(A_i.size()); ++pos) {
        const auto& e = edges[A_i[pos]];
        if (sigma[e.u] == -1 || sigma[e.v] == -1 || pi[e.u] == pi[e.v]) {
            actions[pos] = EdgeAction::Find;
        } else if (f_dsu.find(e.u) != f_dsu.find(e.v)) {
            actions[pos] = EdgeAction::Project;
        } else {
            ambig_pos.push_back(pos);
        }
    }

    const int m_a = static_cast<int>(ambig_pos.size());
    if (m_a == 0) return actions;

    if (stats) stats->estimator_ambig_edges += static_cast<uint64_t>(m_a);

    int s_budget = m_a;
    if (config.estimator_mode == EstimatorMode::DAGTopological) {
        const double raw = config.estimator_sample_c * static_cast<double>(t) *
                           std::sqrt(static_cast<double>(m_a) * static_cast<double>(n));
        const int cap = std::max(
            1, static_cast<int>(std::floor(config.max_sample_fraction * static_cast<double>(m_a))));
        if (raw < static_cast<double>(m_a) || cap < m_a) {
            s_budget = std::min(m_a, std::min(cap, static_cast<int>(std::ceil(raw))));
            s_budget = std::max(1, s_budget);
        }
    }

    // 2. Exact path when budget covers all ambiguous edges.
    if (s_budget >= m_a) {
        std::vector<ForestEdge> forest_edges;
        forest_edges.reserve(F.size());
        for (int f_idx : F) {
            forest_edges.push_back({edges[f_idx].u, edges[f_idx].v, edges[f_idx].orig_id});
        }

        std::vector<std::pair<int, int>> queries;
        queries.reserve(m_a);
        for (int pos : ambig_pos) {
            const auto& e = edges[A_i[pos]];
            queries.push_back({e.u, e.v});
        }

        auto f_max = path_max_queries(
            n, forest_edges, queries,
            [&](int a, int b) { return oracle.less(a, b, CmpCategory::Estimator); });

        for (int j = 0; j < m_a; ++j) {
            const int pos = ambig_pos[j];
            const auto& e = edges[A_i[pos]];
            const int su = edges[sigma[e.u]].orig_id;
            const int sv = edges[sigma[e.v]].orig_id;
            const int s_max = oracle.less(su, sv, CmpCategory::Estimator) ? sv : su;
            const int fe = f_max[j];
            assert(fe != -1);
            const bool f_le_s = (fe == s_max) || oracle.less(fe, s_max, CmpCategory::Estimator);
            actions[pos] = f_le_s ? EdgeAction::Find : EdgeAction::Project;
        }
        if (stats) stats->estimator_queried_edges += static_cast<uint64_t>(m_a);
        return actions;
    }

    // 3. Sub-sampled Comparison-DAG topological-sort estimator (Section 3.5).
    if (stats) {
        ++stats->estimator_dag_invocations;
        stats->estimator_queried_edges += static_cast<uint64_t>(s_budget);
    }

    // Deterministic sample S* without replacement (0 weight comparisons).
    std::vector<int> perm(m_a);
    std::iota(perm.begin(), perm.end(), 0);
    for (int i = 0; i < s_budget; ++i) {
        const uint64_t r = splitmix64(round_seed ^ static_cast<uint64_t>(i));
        const int j = i + static_cast<int>(r % static_cast<uint64_t>(m_a - i));
        std::swap(perm[i], perm[j]);
    }

    // Build ground set T = F U {active spokes}, |T| <= 2n.
    std::vector<int> T_orig_ids;
    T_orig_ids.reserve(F.size() + n);
    for (int f_idx : F) T_orig_ids.push_back(edges[f_idx].orig_id);
    for (int v = 0; v < n; ++v) {
        if (sigma[v] != -1) T_orig_ids.push_back(edges[sigma[v]].orig_id);
    }
    std::sort(T_orig_ids.begin(), T_orig_ids.end());
    T_orig_ids.erase(std::unique(T_orig_ids.begin(), T_orig_ids.end()), T_orig_ids.end());

    const int num_T = static_cast<int>(T_orig_ids.size());
    auto tid_of = [&](int orig_id) -> int {
        auto it = std::lower_bound(T_orig_ids.begin(), T_orig_ids.end(), orig_id);
        assert(it != T_orig_ids.end() && *it == orig_id);
        return static_cast<int>(it - T_orig_ids.begin());
    };

    std::vector<ForestEdge> forest_tid;
    forest_tid.reserve(F.size());
    for (int f_idx : F) {
        forest_tid.push_back({edges[f_idx].u, edges[f_idx].v, tid_of(edges[f_idx].orig_id)});
    }

    std::vector<std::vector<int>> dag_adj(num_T);
    std::vector<int> indeg(num_T, 0);

    auto cmp_record = [&](int tid_a, int tid_b) -> bool {
        if (tid_a == tid_b) return false;
        const bool a_lt_b =
            oracle.less(T_orig_ids[tid_a], T_orig_ids[tid_b], CmpCategory::Estimator);
        if (a_lt_b) {
            dag_adj[tid_a].push_back(tid_b);
            ++indeg[tid_b];
        } else {
            dag_adj[tid_b].push_back(tid_a);
            ++indeg[tid_a];
        }
        return a_lt_b;
    };

    std::vector<std::pair<int, int>> sample_queries;
    sample_queries.reserve(s_budget);
    for (int i = 0; i < s_budget; ++i) {
        const auto& e = edges[A_i[ambig_pos[perm[i]]]];
        sample_queries.push_back({e.u, e.v});
    }

    auto f_max_sample = path_max_queries(n, forest_tid, sample_queries, cmp_record);

    std::vector<EdgeAction> sample_exact_actions;
    if (config.verify_invariants) sample_exact_actions.resize(s_budget);

    for (int i = 0; i < s_budget; ++i) {
        const auto& e = edges[A_i[ambig_pos[perm[i]]]];
        const int tid_su = tid_of(edges[sigma[e.u]].orig_id);
        const int tid_sv = tid_of(edges[sigma[e.v]].orig_id);
        const int tid_smax = cmp_record(tid_su, tid_sv) ? tid_sv : tid_su;
        const int tid_fmax = f_max_sample[i];
        assert(tid_fmax != -1);
        const bool f_le_s = (tid_fmax == tid_smax) || cmp_record(tid_fmax, tid_smax);
        if (config.verify_invariants) {
            sample_exact_actions[i] = f_le_s ? EdgeAction::Find : EdgeAction::Project;
        }
    }

    // Kahn's topological sort on comparison DAG (0 weight comparisons).
    std::vector<int> q;
    q.reserve(num_T);
    for (int i = 0; i < num_T; ++i) {
        if (indeg[i] == 0) q.push_back(i);
    }
    std::vector<int> topo_rank(num_T, -1);
    int next_rank = 0;
    for (size_t head = 0; head < q.size(); ++head) {
        const int u = q[head];
        topo_rank[u] = next_rank++;
        for (int v : dag_adj[u]) {
            if (--indeg[v] == 0) q.push_back(v);
        }
    }
    assert(next_rank == num_T);

    // Resolve ALL ambiguous edges using topo_rank (0 real weight comparisons).
    std::vector<std::pair<int, int>> all_queries;
    all_queries.reserve(m_a);
    for (int pos : ambig_pos) {
        const auto& e = edges[A_i[pos]];
        all_queries.push_back({e.u, e.v});
    }

    auto f_max_all = path_max_queries(
        n, forest_tid, all_queries,
        [&](int tid_a, int tid_b) { return topo_rank[tid_a] < topo_rank[tid_b]; });

    for (int j = 0; j < m_a; ++j) {
        const int pos = ambig_pos[j];
        const auto& e = edges[A_i[pos]];
        const int tid_su = tid_of(edges[sigma[e.u]].orig_id);
        const int tid_sv = tid_of(edges[sigma[e.v]].orig_id);
        const int r_smax = std::max(topo_rank[tid_su], topo_rank[tid_sv]);
        const int r_fmax = topo_rank[f_max_all[j]];
        actions[pos] = (r_smax >= r_fmax) ? EdgeAction::Find : EdgeAction::Project;
    }

    if (config.verify_invariants) {
        for (int i = 0; i < s_budget; ++i) {
            assert(actions[ambig_pos[perm[i]]] == sample_exact_actions[i]);
        }
    }

    return actions;
}

}  // namespace dlmst
