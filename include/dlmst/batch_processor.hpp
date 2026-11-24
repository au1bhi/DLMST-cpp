#pragma once

#include "dlmst/bfprt.hpp"
#include "dlmst/types.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

namespace dlmst {

struct RejectedEdge {
    int edge_idx = -1;
    int first_invalid_level = -1;
};

struct BatchProcessResult {
    std::vector<int> declared_edges;       // D: local edge indices into CertEdge vector
    std::vector<int> new_portals;          // P': newly selected portal vertex IDs
    std::vector<RejectedEdge> rejected;    // R': edges failing admissibility during batch
};

namespace detail {

struct IncidentEntry {
    int edge_idx = -1;
    int orig_id = -1;
    int group = -1;  // Other endpoint vertex
};

struct BatchWorkspace {
    std::vector<int8_t> adm_state;      // 0 = untested, 1 = admissible, -1 = rejected
    std::vector<int> first_invalid_lvl;
    std::vector<int> group_stamp;
    std::vector<int> group_bucket;
    int cur_stamp = 0;

    void init(int n, size_t num_edges) {
        adm_state.assign(num_edges, 0);
        first_invalid_lvl.assign(num_edges, -1);
        group_stamp.assign(n, 0);
        group_bucket.assign(n, -1);
        cur_stamp = 0;
    }
};

template <typename AdmFn>
void select_k_representatives_rec(
    std::vector<IncidentEntry>& items,
    int k,
    BatchWorkspace& ws,
    WeightOracle& oracle,
    AdmFn&& ensure_admissible,
    std::vector<IncidentEntry>& out) {
    if (k <= 0 || items.empty()) return;

    // Drop entries already known to be rejected (0 comparisons).
    size_t write_pos = 0;
    for (size_t i = 0; i < items.size(); ++i) {
        if (ws.adm_state[items[i].edge_idx] != -1) {
            items[write_pos++] = items[i];
        }
    }
    items.resize(write_pos);
    if (items.empty()) return;

    // Base case: m <= k -> test every remaining edge and take lightest admissible per group.
    if (static_cast<int>(items.size()) <= k) {
        const int stamp = ++ws.cur_stamp;
        std::vector<int> seen_groups;
        for (const auto& entry : items) {
            if (!ensure_admissible(entry.edge_idx)) continue;
            if (ws.group_stamp[entry.group] != stamp) {
                ws.group_stamp[entry.group] = stamp;
                ws.group_bucket[entry.group] = static_cast<int>(seen_groups.size());
                seen_groups.push_back(entry.group);
                out.push_back(entry);
            } else {
                const int pos = static_cast<int>(out.size()) -
                                static_cast<int>(seen_groups.size()) +
                                ws.group_bucket[entry.group];
                if (oracle.less(entry.orig_id, out[pos].orig_id, CmpCategory::BatchSelect)) {
                    out[pos] = entry;
                }
            }
        }
        for (int g : seen_groups) ws.group_bucket[g] = -1;
        return;
    }

    // Partition items into L (floor(m/2) lightest) and H (remaining) via BFPRT.
    const size_t half = items.size() / 2;
    assert(half >= 1 && half < items.size());

    auto cmp_entry = [&](const IncidentEntry& a, const IncidentEntry& b) {
        return oracle.less(a.orig_id, b.orig_id, CmpCategory::BatchSelect);
    };
    bfprt_nth_element(items, 0, items.size(), half - 1, cmp_entry);

    // Scan L, testing edges until k distinct groups pass or L is exhausted.
    const int stamp = ++ws.cur_stamp;
    std::vector<int> passed_groups;
    passed_groups.reserve(k);

    for (size_t i = 0; i < half; ++i) {
        const auto& entry = items[i];
        if (ws.group_stamp[entry.group] == stamp) continue;
        if (ensure_admissible(entry.edge_idx)) {
            ws.group_stamp[entry.group] = stamp;
            passed_groups.push_back(entry.group);
            if (static_cast<int>(passed_groups.size()) == k) break;
        }
    }

    const int k_prime = static_cast<int>(passed_groups.size());
    if (k_prime == k) {
        // All k desired representatives lie in L.
        items.resize(half);
        select_k_representatives_rec(items, k, ws, oracle, ensure_admissible, out);
        return;
    }

    // k' < k:
    // Mark passed groups so we can bucket L and filter H without scanning passed_groups.
    for (int bi = 0; bi < k_prime; ++bi) {
        ws.group_bucket[passed_groups[bi]] = bi;
    }

    std::vector<std::vector<IncidentEntry>> group_buckets(k_prime);
    for (size_t i = 0; i < half; ++i) {
        const auto& entry = items[i];
        if (ws.adm_state[entry.edge_idx] == -1) continue;
        const int bi = ws.group_bucket[entry.group];
        if (bi != -1) {
            group_buckets[bi].push_back(entry);
        }
    }

    std::vector<IncidentEntry> h_filtered;
    h_filtered.reserve(items.size() - half);
    for (size_t i = half; i < items.size(); ++i) {
        const auto& entry = items[i];
        if (ws.adm_state[entry.edge_idx] == -1) continue;
        if (ws.group_bucket[entry.group] == -1) {
            h_filtered.push_back(entry);
        }
    }

    // Clean group_bucket before making recursive calls.
    for (int g : passed_groups) {
        ws.group_bucket[g] = -1;
    }

    // Quota-1 search per passed group on L.
    for (int bi = 0; bi < k_prime; ++bi) {
        select_k_representatives_rec(
            group_buckets[bi], 1, ws, oracle, ensure_admissible, out);
    }

    // Recurse on remaining H with quota k - k'.
    select_k_representatives_rec(
        h_filtered, k - k_prime, ws, oracle, ensure_admissible, out);
}

}  // namespace detail

/**
 * Lemma 3.3 standalone wrapper:
 * Selects up to k lightest admissible group representatives from `items`
 * (modifying `items` order in place).
 */
template <typename AdmFn>
std::vector<detail::IncidentEntry> select_k_representatives(
    int num_groups,
    size_t num_edges,
    std::vector<detail::IncidentEntry>& items,
    int k,
    WeightOracle& oracle,
    AdmFn&& is_edge_admissible) {
    detail::BatchWorkspace ws;
    ws.init(num_groups, num_edges);

    auto ensure_adm = [&](int edge_idx) -> bool {
        if (ws.adm_state[edge_idx] == 1) return true;
        if (ws.adm_state[edge_idx] == -1) return false;
        const bool ok = is_edge_admissible(edge_idx);
        ws.adm_state[edge_idx] = ok ? 1 : -1;
        return ok;
    };

    std::vector<detail::IncidentEntry> out;
    detail::select_k_representatives_rec(items, k, ws, oracle, ensure_adm, out);
    return out;
}

/**
 * Batch processing procedure (Section 3.2 / Algorithm 1 line 20):
 *  1. For each vertex v, select up to k_rep lightest admissible representatives
 *     grouped by neighbor endpoint via Lemma 3.3.
 *  2. Declare all selected edges D and track neighbors S_v.
 *  3. Greedily hit every S_v with |S_v| == k_rep using <= 2n/r portals P'.
 *  4. Update spokes sigma[v], pi[v] for vertices with |S_v| == k_rep.
 */
inline BatchProcessResult process_batch(
    int n,
    const std::vector<CertEdge>& edges,
    const std::vector<int>& batch_indices,
    int r,
    int64_t k_rep,
    std::vector<int>& sigma,
    std::vector<int>& pi,
    std::vector<bool>& in_P,
    WeightOracle& oracle) {
    BatchProcessResult result;
    if (n <= 0 || batch_indices.empty()) return result;

    detail::BatchWorkspace ws;
    ws.init(n, edges.size());

    auto ensure_admissible = [&](int edge_idx) -> bool {
        if (ws.adm_state[edge_idx] == 1) return true;
        if (ws.adm_state[edge_idx] == -1) return false;

        const auto& e = edges[edge_idx];
        uint64_t local_cmps = 0;
        auto lt = [&](int a, int b) -> bool {
            if (a == b) return false;
            ++local_cmps;
            return oracle.less_unchecked(a, b);
        };

        for (const auto& step : e.proj_record) {
            if (lt(e.orig_id, step.spoke_u) || lt(e.orig_id, step.spoke_v)) {
                ws.adm_state[edge_idx] = -1;
                ws.first_invalid_lvl[edge_idx] = step.level;
                result.rejected.push_back({edge_idx, step.level});
                oracle.record_cmp(CmpCategory::AdmissibilityFail, local_cmps);
                return false;
            }
        }

        ws.adm_state[edge_idx] = 1;
        oracle.record_cmp(CmpCategory::AdmissibilityPass, local_cmps);
        return true;
    };

    // Build incident lists for B (ignoring self-loops if any).
    std::vector<std::vector<detail::IncidentEntry>> inc(n);
    for (int edge_idx : batch_indices) {
        const auto& e = edges[edge_idx];
        if (e.u == e.v) continue;
        inc[e.u].push_back({edge_idx, e.orig_id, e.v});
        inc[e.v].push_back({edge_idx, e.orig_id, e.u});
    }

    // Per-vertex k-lightest admissible representative selection.
    std::vector<std::vector<int>> Sv(n);
    std::vector<std::vector<int>> Sv_edge(n);
    std::vector<bool> declared_mark(edges.size(), false);

    for (int v = 0; v < n; ++v) {
        if (inc[v].empty()) continue;
        const int k_eff = static_cast<int>(
            std::min<int64_t>(k_rep, static_cast<int64_t>(inc[v].size())));
        std::vector<detail::IncidentEntry> selected;
        selected.reserve(k_eff);
        detail::select_k_representatives_rec(
            inc[v], k_eff, ws, oracle, ensure_admissible, selected);

        Sv[v].reserve(selected.size());
        Sv_edge[v].reserve(selected.size());
        for (const auto& entry : selected) {
            Sv[v].push_back(entry.group);
            Sv_edge[v].push_back(entry.edge_idx);
            if (!declared_mark[entry.edge_idx]) {
                declared_mark[entry.edge_idx] = true;
                result.declared_edges.push_back(entry.edge_idx);
            }
        }
    }

    // Select portal hitting set P' covering every S_v of size == k_rep.
    std::vector<bool> uncovered(n, false);
    std::vector<int> first_portal(n, -1);
    std::vector<int> first_portal_edge(n, -1);
    std::vector<std::vector<std::pair<int, int>>> in_sets(n);  // x -> list of (v, edge_idx)
    std::vector<int> cnt(n, 0);

    bool any_full = false;
    for (int v = 0; v < n; ++v) {
        if (static_cast<int64_t>(Sv[v].size()) == k_rep) {
            uncovered[v] = true;
            any_full = true;
            for (size_t i = 0; i < Sv[v].size(); ++i) {
                const int x = Sv[v][i];
                const int eidx = Sv_edge[v][i];
                in_sets[x].push_back({v, eidx});
                ++cnt[x];
            }
        }
    }

    if (any_full) {
        std::vector<int> heavy_queue;
        for (int x = 0; x < n; ++x) {
            if (cnt[x] >= r) heavy_queue.push_back(x);
        }

        auto select_portal = [&](int x) {
            if (!in_P[x]) {
                in_P[x] = true;
                result.new_portals.push_back(x);
            }
            for (const auto& [v, eidx] : in_sets[x]) {
                if (!uncovered[v]) continue;
                uncovered[v] = false;
                first_portal[v] = x;
                first_portal_edge[v] = eidx;
                for (int y : Sv[v]) {
                    --cnt[y];
                }
            }
        };

        for (size_t head = 0; head < heavy_queue.size(); ++head) {
            const int x = heavy_queue[head];
            if (cnt[x] >= r) {
                select_portal(x);
            }
        }

        for (int v = 0; v < n; ++v) {
            if (uncovered[v]) {
                assert(!Sv[v].empty());
                select_portal(Sv[v][0]);
            }
        }

        // Update spokes for each v with |S_v| == k_rep.
        for (int v = 0; v < n; ++v) {
            if (first_portal[v] == -1) continue;
            const int cand_idx = first_portal_edge[v];
            const int portal_v = first_portal[v];
            if (sigma[v] == -1) {
                sigma[v] = cand_idx;
                pi[v] = portal_v;
            } else if (oracle.less(
                           edges[cand_idx].orig_id,
                           edges[sigma[v]].orig_id,
                           CmpCategory::BatchSpoke)) {
                sigma[v] = cand_idx;
                pi[v] = portal_v;
            }
        }
    }

    return result;
}

}  // namespace dlmst
