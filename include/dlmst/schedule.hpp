#pragma once

#include "dlmst/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace dlmst {

/**
 * Deterministic comparison-free schedule (Section 3.4 / Algorithm 1 line 5).
 * Assigns each edge in `edges` to a round in [0, num_rounds).
 */
inline std::vector<int> build_schedule(
    int num_rounds,
    uint64_t call_seed,
    const std::vector<CertEdge>& edges) {
    std::vector<int> sched(edges.size(), 0);
    if (num_rounds <= 1) return sched;
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        const uint64_t key =
            call_seed ^
            (static_cast<uint64_t>(i) * 0x9e3779b97f4a7c15ULL) ^
            (static_cast<uint64_t>(edges[i].orig_id) * 0xbf58476d1ce4e5b9ULL);
        sched[i] = static_cast<int>(splitmix64(key) % static_cast<uint64_t>(num_rounds));
    }
    return sched;
}

struct ExactScheduleResult {
    std::vector<int> best_schedule;
    int min_worst_bad_edges = 0;
    double lemma_3_11_bound = 0.0;
    size_t num_schedules_checked = 0;
    size_t num_weight_orders_checked = 0;
};

/**
 * Constructive small-instance verification of Lemma 3.11:
 * For a small unweighted graph (e.g., m <= 6, t = 2), enumerates all (t^2)^m
 * schedules and all m! weight permutations, simulating exact round preferences
 * to verify that there exists a deterministic schedule rho* achieving small
 * worst-case bad-round count across all weight orderings.
 */
inline ExactScheduleResult find_exact_schedule_small(
    int n,
    const std::vector<Edge>& edges,
    int t) {
    const int m = static_cast<int>(edges.size());
    const int num_rounds = t * t;
    const int r = 2 * t + 2;
    const int k_rep = r * r;

    // Enumerate all m! edge weight rank permutations.
    std::vector<std::vector<int>> perms;
    {
        std::vector<int> p(m);
        std::iota(p.begin(), p.end(), 0);
        do {
            perms.push_back(p);
        } while (std::next_permutation(p.begin(), p.end()));
    }

    auto compute_exact_msf_small = [&](const std::vector<int>& edge_indices,
                                       const std::vector<int>& rank) {
        std::vector<int> sorted = edge_indices;
        std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
            return rank[a] < rank[b];
        });
        DSU dsu(n);
        std::vector<int> tree;
        for (int idx : sorted) {
            if (dsu.unite(edges[idx].u, edges[idx].v)) {
                tree.push_back(idx);
            }
        }
        return tree;
    };

    auto path_max_in_forest = [&](const std::vector<int>& forest,
                                  int u,
                                  int v,
                                  const std::vector<int>& rank) -> int {
        if (u == v) return -1;
        std::vector<std::vector<std::pair<int, int>>> adj(n);
        for (int idx : forest) {
            adj[edges[idx].u].push_back({edges[idx].v, idx});
            adj[edges[idx].v].push_back({edges[idx].u, idx});
        }
        std::vector<int> parent_v(n, -1), parent_e(n, -1);
        std::vector<int> q = {u};
        parent_v[u] = u;
        for (size_t head = 0; head < q.size(); ++head) {
            int cur = q[head];
            if (cur == v) break;
            for (const auto& [nxt, eid] : adj[cur]) {
                if (parent_v[nxt] == -1) {
                    parent_v[nxt] = cur;
                    parent_e[nxt] = eid;
                    q.push_back(nxt);
                }
            }
        }
        if (parent_v[v] == -1) return -1;
        int best = -1;
        for (int cur = v; cur != u; cur = parent_v[cur]) {
            const int eid = parent_e[cur];
            if (best == -1 || rank[eid] > rank[best]) best = eid;
        }
        return best;
    };

    ExactScheduleResult res;
    res.min_worst_bad_edges = m + 1;
    res.num_weight_orders_checked = perms.size();

    size_t total_sched = 1;
    for (int i = 0; i < m; ++i) total_sched *= static_cast<size_t>(num_rounds);
    res.num_schedules_checked = total_sched;

    std::vector<int> sched(m, 0);
    for (size_t s_code = 0; s_code < total_sched; ++s_code) {
        size_t tmp = s_code;
        for (int i = 0; i < m; ++i) {
            sched[i] = static_cast<int>(tmp % static_cast<size_t>(num_rounds));
            tmp /= static_cast<size_t>(num_rounds);
        }

        int worst_bad = 0;
        for (const auto& rank : perms) {
            // Simulate t blocks of t rounds on this small graph under `rank` and `sched`.
            std::vector<int> F;
            std::vector<int> sigma(n, -1);
            std::vector<int> pi(n, -1);
            std::vector<std::vector<int>> pref_history(num_rounds, std::vector<int>(m, 0));

            for (int blk = 0; blk < t; ++blk) {
                std::vector<int> B;
                for (int rb = 0; rb < t; ++rb) {
                    const int rnd = blk * t + rb;
                    for (int e = 0; e < m; ++e) {
                        const int u = edges[e].u;
                        const int v = edges[e].v;
                        int pref = 0;  // 0 = Find, 1 = Project
                        if (sigma[u] != -1 && sigma[v] != -1 && pi[u] != pi[v]) {
                            const int fmax = path_max_in_forest(F, u, v, rank);
                            const int smax = (rank[sigma[u]] > rank[sigma[v]]) ? sigma[u] : sigma[v];
                            if (fmax == -1 || rank[fmax] > rank[smax]) {
                                pref = 1;
                            }
                        }
                        pref_history[rnd][e] = pref;
                    }
                    std::vector<int> projected;
                    for (int e = 0; e < m; ++e) {
                        if (sched[e] == rnd) {
                            if (pref_history[rnd][e] == 0) {
                                B.push_back(e);
                            } else {
                                projected.push_back(e);
                            }
                        }
                    }
                    if (!projected.empty()) {
                        // In a single level without deeper portal recursion, valid projected
                        // cycles resolve back into F via exact MSF on portals; here simulate
                        // exact MSF incorporation of projected edges whose projection is valid.
                        std::vector<int> cand = F;
                        for (int e : projected) {
                            const int smax =
                                (rank[sigma[edges[e].u]] > rank[sigma[edges[e].v]])
                                    ? sigma[edges[e].u]
                                    : sigma[edges[e].v];
                            if (rank[smax] <= rank[e]) cand.push_back(e);
                        }
                        F = compute_exact_msf_small(cand, rank);
                    }
                }

                // Batch processing on B: each vertex picks up to k_rep lightest neighbors.
                std::vector<int> D;
                std::vector<std::vector<int>> Sv(n);
                for (int v = 0; v < n; ++v) {
                    std::vector<int> inc;
                    for (int e : B) {
                        if (edges[e].u == v || edges[e].v == v) inc.push_back(e);
                    }
                    std::sort(inc.begin(), inc.end(), [&](int a, int b) {
                        return rank[a] < rank[b];
                    });
                    for (int idx : inc) {
                        const int other = (edges[idx].u == v) ? edges[idx].v : edges[idx].u;
                        if (std::find(Sv[v].begin(), Sv[v].end(), other) == Sv[v].end()) {
                            Sv[v].push_back(other);
                            D.push_back(idx);
                            if (static_cast<int>(Sv[v].size()) == k_rep) break;
                        }
                    }
                }
                if (!D.empty()) {
                    std::vector<int> cand = F;
                    cand.insert(cand.end(), D.begin(), D.end());
                    F = compute_exact_msf_small(cand, rank);
                }
                for (int v = 0; v < n; ++v) {
                    if (static_cast<int>(Sv[v].size()) == k_rep) {
                        const int x = Sv[v][0];
                        for (int idx : D) {
                            if ((edges[idx].u == v && edges[idx].v == x) ||
                                (edges[idx].v == v && edges[idx].u == x)) {
                                if (sigma[v] == -1 || rank[idx] < rank[sigma[v]]) {
                                    sigma[v] = idx;
                                    pi[v] = x;
                                }
                                break;
                            }
                        }
                    }
                }
            }

            // Count bad-round edges per Definition 3.4.
            int bad = 0;
            for (int e = 0; e < m; ++e) {
                const int r_assigned = sched[e];
                const int b_assigned = r_assigned / t;
                const int pref = pref_history[r_assigned][e];
                bool is_bad = false;
                if (pref == 1) {
                    // Last round in its block with Project preference.
                    int last_p = -1;
                    for (int rb = 0; rb < t; ++rb) {
                        if (pref_history[b_assigned * t + rb][e] == 1) last_p = b_assigned * t + rb;
                    }
                    if (r_assigned == last_p) is_bad = true;
                } else {
                    // Find preference in the last block containing any Find preference.
                    int last_f_blk = -1;
                    for (int blk = 0; blk < t; ++blk) {
                        for (int rb = 0; rb < t; ++rb) {
                            if (pref_history[blk * t + rb][e] == 0) last_f_blk = blk;
                        }
                    }
                    if (b_assigned == last_f_blk) is_bad = true;
                }
                if (is_bad) ++bad;
            }
            if (bad > worst_bad) worst_bad = bad;
            if (worst_bad >= res.min_worst_bad_edges) break;
        }

        if (worst_bad < res.min_worst_bad_edges) {
            res.min_worst_bad_edges = worst_bad;
            res.best_schedule = sched;
        }
    }

    // Lemma 3.11 bound: (log(|H_approx|) + 2m/t) / log(2), where |H| >= perms.size()
    // gives the concrete comparison-permutation lower bound; full paper upper bound uses
    // Lemma 3.10's bound on log|H|.
    const double log_h_est = std::log(static_cast<double>(perms.size()));
    res.lemma_3_11_bound =
        (log_h_est + 2.0 * static_cast<double>(m) / static_cast<double>(t)) / std::log(2.0);
    return res;
}

}  // namespace dlmst
