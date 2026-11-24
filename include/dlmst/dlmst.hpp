#pragma once

#include "dlmst/action_estimator.hpp"
#include "dlmst/batch_processor.hpp"
#include "dlmst/boruvka.hpp"
#include "dlmst/path_max.hpp"
#include "dlmst/schedule.hpp"
#include "dlmst/types.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

namespace dlmst {

namespace detail {

struct CertificateResult {
    std::vector<int> forest_indices;        // F: local indices into cert_edges
    std::vector<RejectedEdge> rejected;     // R: edges whose first invalid projection is < depth
};

// Forward declaration for mutual recursion between Algorithm 1 and Algorithm 2.
inline std::vector<int> msf_recursive(
    int n,
    const std::vector<Edge>& edges,
    WeightOracle& oracle,
    const DLMSTConfig& config,
    DLMSTStats& stats,
    int msf_depth);

inline CertificateResult certificate_recursive(
    int n,
    const std::vector<CertEdge>& cert_edges,
    int depth,
    WeightOracle& oracle,
    const DLMSTConfig& config,
    DLMSTStats& stats) {
    CertificateResult result;
    if (cert_edges.empty() || n <= 1) return result;

    ++stats.cert_calls;
    if (depth > stats.max_cert_depth) stats.max_cert_depth = depth;

    int t = 2;
    int r = 5;
    int64_t k_rep = 25;
    config.compute_level_params(depth, t, r, k_rep);

    const int m = static_cast<int>(cert_edges.size());
    const int num_rounds = t * t;

    // Map orig_id -> local index in cert_edges (0 weight comparisons).
    std::vector<std::pair<int, int>> orig_to_local(m);
    for (int i = 0; i < m; ++i) {
        orig_to_local[i] = {cert_edges[i].orig_id, i};
    }
    std::sort(orig_to_local.begin(), orig_to_local.end());

    auto lookup_local = [&](int orig_id) -> int {
        auto it = std::lower_bound(
            orig_to_local.begin(),
            orig_to_local.end(),
            std::make_pair(orig_id, -1));
        assert(it != orig_to_local.end() && it->first == orig_id);
        return it->second;
    };

    // Deterministic schedule assigning each edge to one of the t^2 rounds.
    const uint64_t call_seed =
        config.schedule_seed ^
        (static_cast<uint64_t>(depth) * 0x9e3779b97f4a7c15ULL) ^
        (static_cast<uint64_t>(m) * 0xbf58476d1ce4e5b9ULL) ^
        (static_cast<uint64_t>(n) * 0x94d049bb133111ebULL);
    const std::vector<int> sched = build_schedule(num_rounds, call_seed, cert_edges);

    std::vector<std::vector<int>> round_edges(num_rounds);
    for (int i = 0; i < m; ++i) {
        round_edges[sched[i]].push_back(i);
    }

    std::vector<int> F;
    std::vector<int> sigma(n, -1);
    std::vector<int> pi(n, -1);
    std::vector<bool> in_P(n, false);

    // Helper implementing F <- MSF(G[F U extra]) via an exact recursive MSF call.
    std::vector<int> v_remap(n, -1);
    std::vector<bool> edge_seen(m, false);
    auto insert_into_F = [&](const std::vector<int>& extra) {
        if (extra.empty()) return;

        std::vector<int> union_local;
        union_local.reserve(F.size() + extra.size());
        for (int idx : F) {
            if (!edge_seen[idx]) {
                edge_seen[idx] = true;
                union_local.push_back(idx);
            }
        }
        for (int idx : extra) {
            if (!edge_seen[idx]) {
                edge_seen[idx] = true;
                union_local.push_back(idx);
            }
        }

        int sub_n = 0;
        std::vector<Edge> sub_edges;
        sub_edges.reserve(union_local.size());
        std::vector<int> touched_verts;
        touched_verts.reserve(union_local.size() * 2);

        for (int idx : union_local) {
            edge_seen[idx] = false;
            const auto& ce = cert_edges[idx];
            if (ce.u == ce.v) continue;
            if (v_remap[ce.u] == -1) {
                v_remap[ce.u] = sub_n++;
                touched_verts.push_back(ce.u);
            }
            if (v_remap[ce.v] == -1) {
                v_remap[ce.v] = sub_n++;
                touched_verts.push_back(ce.v);
            }
            sub_edges.push_back({v_remap[ce.u], v_remap[ce.v], ce.orig_id});
        }
        for (int v : touched_verts) v_remap[v] = -1;

        if (sub_edges.empty()) return;

        const std::vector<int> kept_orig =
            msf_recursive(sub_n, sub_edges, oracle, config, stats, /*msf_depth=*/1);
        F.clear();
        F.reserve(kept_orig.size());
        for (int oid : kept_orig) {
            F.push_back(lookup_local(oid));
        }
    };

    std::vector<int> portal_remap(n, -1);
    std::vector<bool> batch_rej_mark;
    if (config.verify_invariants) batch_rej_mark.assign(m, false);

    // Main loop: t blocks of t rounds each (Algorithm 1 lines 6-21).
    for (int blk = 0; blk < t; ++blk) {
        std::vector<int> B;

        for (int rb = 0; rb < t; ++rb) {
            const int rnd = blk * t + rb;
            const auto& A_i = round_edges[rnd];
            if (A_i.empty()) continue;

            const uint64_t rnd_seed =
                call_seed ^ (static_cast<uint64_t>(rnd + 1) * 0x94d049bb133111ebULL);
            const auto actions = estimate_preferred_actions(
                n, cert_edges, F, sigma, pi, A_i, t, rnd_seed, config, oracle, &stats);

            std::vector<int> proj_indices;
            for (size_t pos = 0; pos < A_i.size(); ++pos) {
                if (actions[pos] == EdgeAction::Find) {
                    B.push_back(A_i[pos]);
                } else {
                    proj_indices.push_back(A_i[pos]);
                }
            }

            if (proj_indices.empty()) continue;

            // Build projected child graph G' on active portals.
            int child_n = 0;
            std::vector<int> used_portals;
            std::vector<CertEdge> child_edges;
            child_edges.reserve(proj_indices.size());

            for (int e_idx : proj_indices) {
                const auto& ce = cert_edges[e_idx];
                const int pu = pi[ce.u];
                const int pv = pi[ce.v];
                assert(pu != -1 && pv != -1 && pu != pv);

                if (portal_remap[pu] == -1) {
                    portal_remap[pu] = child_n++;
                    used_portals.push_back(pu);
                }
                if (portal_remap[pv] == -1) {
                    portal_remap[pv] = child_n++;
                    used_portals.push_back(pv);
                }

                CertEdge pe;
                pe.u = portal_remap[pu];
                pe.v = portal_remap[pv];
                pe.orig_id = ce.orig_id;
                pe.parent_edge_idx = e_idx;
                pe.proj_record = ce.proj_record;
                pe.proj_record.push_back({
                    depth,
                    cert_edges[sigma[ce.u]].orig_id,
                    cert_edges[sigma[ce.v]].orig_id
                });
                child_edges.push_back(std::move(pe));
            }
            for (int p : used_portals) portal_remap[p] = -1;

            stats.total_projected_edges += static_cast<uint64_t>(child_edges.size());

            auto child_res = certificate_recursive(
                child_n, child_edges, depth + 1, oracle, config, stats);

            std::vector<int> to_insert;
            to_insert.reserve(child_res.forest_indices.size() + child_res.rejected.size());
            for (int cf_idx : child_res.forest_indices) {
                to_insert.push_back(child_edges[cf_idx].parent_edge_idx);
            }
            for (const auto& rref : child_res.rejected) {
                const int parent_idx = child_edges[rref.edge_idx].parent_edge_idx;
                if (rref.first_invalid_level == depth) {
                    to_insert.push_back(parent_idx);  // R''
                } else {
                    assert(rref.first_invalid_level < depth);
                    result.rejected.push_back({parent_idx, rref.first_invalid_level});
                }
            }
            insert_into_F(to_insert);
        }

        // Process cumulative batch B at the end of the block (Algorithm 1 lines 20-21).
        if (!B.empty()) {
            stats.total_batch_edges += static_cast<uint64_t>(B.size());
            auto bres = process_batch(n, cert_edges, B, r, k_rep, sigma, pi, in_P, oracle);

            stats.total_portals += static_cast<uint64_t>(bres.new_portals.size());
            stats.total_declared_edges += static_cast<uint64_t>(bres.declared_edges.size());
            stats.total_rejected_edges += static_cast<uint64_t>(bres.rejected.size());

            for (const auto& rref : bres.rejected) {
                assert(rref.first_invalid_level < depth);
                result.rejected.push_back(rref);
            }

            insert_into_F(bres.declared_edges);

            // Invariant check: Lemma 3.2 guarantee on B \ R'.
            if (config.verify_invariants) {
                for (const auto& rref : bres.rejected) {
                    batch_rej_mark[rref.edge_idx] = true;
                }
                std::vector<ForestEdge> f_edges;
                f_edges.reserve(F.size());
                for (int f_idx : F) {
                    f_edges.push_back({cert_edges[f_idx].u, cert_edges[f_idx].v, cert_edges[f_idx].orig_id});
                }
                std::vector<int> need_f_check;
                std::vector<std::pair<int, int>> f_queries;
                for (int e_idx : B) {
                    if (batch_rej_mark[e_idx]) continue;
                    const auto& ce = cert_edges[e_idx];
                    const bool valid_proj =
                        sigma[ce.u] != -1 && sigma[ce.v] != -1 &&
                        !oracle.less_unchecked(ce.orig_id, cert_edges[sigma[ce.u]].orig_id) &&
                        !oracle.less_unchecked(ce.orig_id, cert_edges[sigma[ce.v]].orig_id);
                    if (!valid_proj) {
                        need_f_check.push_back(e_idx);
                        f_queries.push_back({ce.u, ce.v});
                    }
                }
                if (!f_queries.empty()) {
                    auto f_max = path_max_queries(
                        n, f_edges, f_queries,
                        [&](int a, int b) { return oracle.less_unchecked(a, b); });
                    for (size_t i = 0; i < need_f_check.size(); ++i) {
                        const int e_idx = need_f_check[i];
                        assert(f_max[i] != -1);
                        assert(!oracle.less_unchecked(cert_edges[e_idx].orig_id, f_max[i]));
                    }
                }
                for (const auto& rref : bres.rejected) {
                    batch_rej_mark[rref.edge_idx] = false;
                }
            }
        }
    }

    // Invariant check: Claim 3.5 (every edge in F must be admissible).
    if (config.verify_invariants) {
        for (int f_idx : F) {
            const auto& ce = cert_edges[f_idx];
            for (const auto& step : ce.proj_record) {
                assert(!oracle.less_unchecked(ce.orig_id, step.spoke_u));
                assert(!oracle.less_unchecked(ce.orig_id, step.spoke_v));
            }
        }
    }

    if (depth == 1) {
        assert(result.rejected.empty());
    }

    result.forest_indices = std::move(F);
    return result;
}

inline std::vector<int> msf_recursive(
    int n,
    const std::vector<Edge>& edges,
    WeightOracle& oracle,
    const DLMSTConfig& config,
    DLMSTStats& stats,
    int msf_depth) {
    ++stats.msf_calls;
    if (msf_depth > stats.max_msf_depth) stats.max_msf_depth = msf_depth;

    if (n <= 1 || edges.empty()) return {};

    // Line 1-5: Perform k >= 1 Boruvka contraction phases.
    const int k_boruvka = std::max(1, config.boruvka_steps);
    auto bres = run_boruvka_phases(n, edges, k_boruvka, oracle);
    std::vector<int> result = std::move(bres.contracted_orig_ids);
    if (bres.contracted_edges.empty() || bres.contracted_n <= 1) {
        return result;
    }

    // Line 6: Run Certificate(G') to obtain forest F.
    const int gp_n = bres.contracted_n;
    const int gp_m = static_cast<int>(bres.contracted_edges.size());

    std::vector<CertEdge> cert_edges(gp_m);
    for (int i = 0; i < gp_m; ++i) {
        cert_edges[i].u = bres.contracted_edges[i].u;
        cert_edges[i].v = bres.contracted_edges[i].v;
        cert_edges[i].orig_id = bres.contracted_edges[i].orig_id;
    }

    auto cres = certificate_recursive(gp_n, cert_edges, /*depth=*/1, oracle, config, stats);

    // Line 7: Compute U(F) using Komlos/King O(n + q) path maximum queries.
    std::vector<bool> in_F(gp_m, false);
    std::vector<ForestEdge> forest_edges;
    forest_edges.reserve(cres.forest_indices.size());
    std::vector<Edge> next_edges;
    next_edges.reserve(gp_m);

    for (int f_idx : cres.forest_indices) {
        in_F[f_idx] = true;
        const auto& ce = cert_edges[f_idx];
        forest_edges.push_back({ce.u, ce.v, ce.orig_id});
        next_edges.push_back({ce.u, ce.v, ce.orig_id});
    }

    std::vector<int> check_indices;
    std::vector<std::pair<int, int>> check_queries;
    check_indices.reserve(gp_m - cres.forest_indices.size());
    check_queries.reserve(gp_m - cres.forest_indices.size());
    for (int i = 0; i < gp_m; ++i) {
        if (!in_F[i]) {
            check_indices.push_back(i);
            check_queries.push_back({cert_edges[i].u, cert_edges[i].v});
        }
    }

    if (!check_queries.empty()) {
        auto f_max = path_max_queries(
            gp_n, forest_edges, check_queries,
            [&](int a, int b) { return oracle.less(a, b, CmpCategory::FilterVerify); });

        for (size_t qi = 0; qi < check_indices.size(); ++qi) {
            const int idx = check_indices[qi];
            const int oid = cert_edges[idx].orig_id;
            const int fmax_oid = f_max[qi];
            // Supported iff endpoints are connected in F and max path edge <= w(e).
            // Otherwise unsupported -> retain in U(F).
            if (fmax_oid == -1 || oracle.less(oid, fmax_oid, CmpCategory::FilterVerify)) {
                next_edges.push_back({cert_edges[idx].u, cert_edges[idx].v, oid});
            }
        }
    }

    if (msf_depth == 0) {
        stats.root_edges_before_filter = static_cast<uint64_t>(gp_m);
        stats.root_edges_after_filter = static_cast<uint64_t>(next_edges.size());
    }

    // Line 8: Recurse on G'[F U U(F)].
    auto rec_ids = msf_recursive(gp_n, next_edges, oracle, config, stats, msf_depth + 1);
    result.insert(result.end(), rec_ids.begin(), rec_ids.end());
    return result;
}

}  // namespace detail

/**
 * Computes the exact Minimum Spanning Forest (MSF) of an undirected weighted
 * multigraph `g` using our C++20 reproduction of Khanna-Putterman-Song (2026) DLMST.
 */
inline MSFResult dlmst_msf(const Graph& g, const DLMSTConfig& config = {}) {
    MSFResult res;
    if (g.n <= 1 || g.edges.empty()) return res;

    WeightOracle oracle(g.weights, &res.stats);
    res.edge_ids = detail::msf_recursive(g.n, g.edges, oracle, config, res.stats, /*msf_depth=*/0);
    std::sort(res.edge_ids.begin(), res.edge_ids.end());
    return res;
}

}  // namespace dlmst
