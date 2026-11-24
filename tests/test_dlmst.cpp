#include "dlmst/action_estimator.hpp"
#include "dlmst/batch_processor.hpp"
#include "dlmst/bfprt.hpp"
#include "dlmst/boruvka.hpp"
#include "dlmst/dlmst.hpp"
#include "dlmst/kruskal.hpp"
#include "dlmst/path_max.hpp"
#include "dlmst/schedule.hpp"
#include "dlmst/types.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <utility>
#include <vector>

using namespace dlmst;

namespace {

int naive_path_max(
    int n,
    const std::vector<ForestEdge>& forest,
    int u,
    int v,
    const std::vector<double>& weights) {
    if (u == v) return -1;
    std::vector<std::vector<std::pair<int, int>>> adj(n);
    for (const auto& fe : forest) {
        adj[fe.u].push_back({fe.v, fe.edge_id});
        adj[fe.v].push_back({fe.u, fe.edge_id});
    }
    std::vector<int> pv(n, -1), pe(n, -1);
    std::vector<int> q = {u};
    pv[u] = u;
    for (size_t head = 0; head < q.size(); ++head) {
        int cur = q[head];
        if (cur == v) break;
        for (const auto& [nxt, eid] : adj[cur]) {
            if (pv[nxt] == -1) {
                pv[nxt] = cur;
                pe[nxt] = eid;
                q.push_back(nxt);
            }
        }
    }
    if (pv[v] == -1) return -1;
    int best = -1;
    for (int cur = v; cur != u; cur = pv[cur]) {
        const int eid = pe[cur];
        if (best == -1) {
            best = eid;
        } else if (weights[best] < weights[eid] ||
                   (weights[best] == weights[eid] && best < eid)) {
            best = eid;
        }
    }
    return best;
}

Graph make_random_graph(
    int n,
    int m,
    uint64_t seed,
    bool allow_parallel = false,
    bool equal_weights = false) {
    std::mt19937_64 rng(seed);
    Graph g;
    g.n = n;
    std::uniform_int_distribution<int> vdist(0, std::max(0, n - 1));
    std::uniform_real_distribution<double> wdist(-1000.0, 1000.0);

    std::set<std::pair<int, int>> seen;
    if (n >= 2) {
        // Connected spanning backbone first if m >= n - 1
        int backbone = std::min(m, n - 1);
        for (int i = 0; i < backbone; ++i) {
            int p = static_cast<int>(rng() % static_cast<uint64_t>(i + 1));
            int u = i + 1;
            int v = p;
            if (!allow_parallel) {
                seen.insert({std::min(u, v), std::max(u, v)});
            }
            int oid = static_cast<int>(g.edges.size());
            g.edges.push_back({u, v, oid});
            g.weights.push_back(equal_weights ? 42.0 : wdist(rng));
        }
    }

    int tries = 0;
    while (static_cast<int>(g.edges.size()) < m && tries < m * 20 + 100) {
        ++tries;
        if (n < 2) break;
        int u = vdist(rng);
        int v = vdist(rng);
        if (u == v) continue;
        if (!allow_parallel) {
            auto key = std::make_pair(std::min(u, v), std::max(u, v));
            if (seen.count(key)) continue;
            seen.insert(key);
        }
        int oid = static_cast<int>(g.edges.size());
        g.edges.push_back({u, v, oid});
        g.weights.push_back(equal_weights ? 42.0 : wdist(rng));
    }
    return g;
}

Graph make_complete_graph(int n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> wdist(0.0, 100.0);
    Graph g;
    g.n = n;
    for (int u = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v) {
            int oid = static_cast<int>(g.edges.size());
            g.edges.push_back({u, v, oid});
            g.weights.push_back(wdist(rng));
        }
    }
    return g;
}

Graph make_grid_graph(int rows, int cols, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> wdist(0.0, 1000.0);
    Graph g;
    g.n = rows * cols;
    auto id = [&](int r, int c) { return r * cols + c; };
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            if (r + 1 < rows) {
                int oid = static_cast<int>(g.edges.size());
                g.edges.push_back({id(r, c), id(r + 1, c), oid});
                g.weights.push_back(wdist(rng));
            }
            if (c + 1 < cols) {
                int oid = static_cast<int>(g.edges.size());
                g.edges.push_back({id(r, c), id(r, c + 1), oid});
                g.weights.push_back(wdist(rng));
            }
        }
    }
    return g;
}

Graph make_ba_hub_graph(int n, int attach_k, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> wdist(0.0, 100.0);
    Graph g;
    g.n = n;
    int init_k = std::min(n, attach_k + 1);
    std::vector<int> endpoints;
    for (int u = 0; u < init_k; ++u) {
        for (int v = u + 1; v < init_k; ++v) {
            int oid = static_cast<int>(g.edges.size());
            g.edges.push_back({u, v, oid});
            g.weights.push_back(wdist(rng));
            endpoints.push_back(u);
            endpoints.push_back(v);
        }
    }
    for (int u = init_k; u < n; ++u) {
        std::set<int> chosen;
        int target = std::min(attach_k, u);
        while (static_cast<int>(chosen.size()) < target) {
            int pick = endpoints[rng() % endpoints.size()];
            chosen.insert(pick);
        }
        for (int v : chosen) {
            int oid = static_cast<int>(g.edges.size());
            g.edges.push_back({u, v, oid});
            g.weights.push_back(wdist(rng));
            endpoints.push_back(u);
            endpoints.push_back(v);
        }
    }
    return g;
}

void assert_matches_kruskal(const Graph& g, const DLMSTConfig& cfg, const char* label,
                            DLMSTStats* out_stats = nullptr) {
    auto expected = kruskal_msf(g);
    auto actual = dlmst_msf(g, cfg);
    if (actual.edge_ids != expected) {
        std::cerr << "FAIL [" << label << "]: expected " << expected.size()
                  << " edges, got " << actual.edge_ids.size() << "\n";
        std::abort();
    }
    if (out_stats) *out_stats = actual.stats;
}

void test_bfprt() {
    std::cout << "[1/6] Testing BFPRT deterministic selection...\n";
    std::mt19937_64 rng(12345);
    for (int trial = 0; trial < 200; ++trial) {
        int n = 1 + (trial % 300);
        std::vector<std::pair<int, int>> arr(n);
        for (int i = 0; i < n; ++i) {
            int val = (trial % 4 == 0) ? i :
                      (trial % 4 == 1) ? (n - i) :
                      (trial % 4 == 2) ? static_cast<int>(rng() % 20) :
                                         static_cast<int>(rng() % 10000);
            arr[i] = {val, i};
        }
        size_t nth = rng() % static_cast<size_t>(n);
        auto ref = arr;
        auto cmp = [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return a < b;
        };
        std::sort(ref.begin(), ref.end(), cmp);
        bfprt_nth_element(arr, 0, arr.size(), nth, cmp);
        assert(arr[nth] == ref[nth]);
        for (size_t i = 0; i < nth; ++i) assert(arr[i] < arr[nth]);
        for (size_t i = nth + 1; i < arr.size(); ++i) assert(arr[nth] < arr[i]);
    }

    // Comparison bound check on n = 10,000.
    std::vector<int> big(10000);
    std::iota(big.begin(), big.end(), 0);
    std::shuffle(big.begin(), big.end(), rng);
    uint64_t cmps = 0;
    bfprt_nth_element(big, 0, big.size(), 5000, [&](int a, int b) {
        ++cmps;
        return a < b;
    });
    assert(big[5000] == 5000);
    assert(cmps < 15ULL * 10000ULL);
    std::cout << "  PASS (n=10000 BFPRT comparisons = " << cmps
              << ", ratio = " << static_cast<double>(cmps) / 10000.0 << "x)\n";
}

void test_path_max() {
    std::cout << "[2/6] Testing King + Komlos/Hagerup path_max_queries...\n";
    std::mt19937_64 rng(777);

    auto run_forest_case = [&](int n, const std::vector<ForestEdge>& forest,
                               const std::vector<double>& weights,
                               const std::vector<std::pair<int, int>>& queries,
                               bool check_dag) {
        std::vector<std::vector<int>> dag_adj(weights.size());
        std::vector<int> indeg(weights.size(), 0);
        uint64_t cmps = 0;
        auto ans = path_max_queries(n, forest, queries, [&](int a, int b) {
            if (a == b) return false;
            ++cmps;
            bool res = (weights[a] < weights[b]) || (weights[a] == weights[b] && a < b);
            if (check_dag) {
                if (res) {
                    dag_adj[a].push_back(b);
                    ++indeg[b];
                } else {
                    dag_adj[b].push_back(a);
                    ++indeg[a];
                }
            }
            return res;
        });

        for (size_t i = 0; i < queries.size(); ++i) {
            int expected = naive_path_max(n, forest, queries[i].first, queries[i].second, weights);
            if (ans[i] != expected) {
                std::cerr << "PathMax mismatch at q=" << i << ": got " << ans[i]
                          << ", expected " << expected << "\n";
                std::abort();
            }
        }

        if (check_dag) {
            std::vector<int> q;
            for (size_t i = 0; i < weights.size(); ++i) {
                if (indeg[i] == 0) q.push_back(static_cast<int>(i));
            }
            std::vector<int> topo_rank(weights.size(), -1);
            int nxt = 0;
            for (size_t head = 0; head < q.size(); ++head) {
                int u = q[head];
                topo_rank[u] = nxt++;
                for (int v : dag_adj[u]) {
                    if (--indeg[v] == 0) q.push_back(v);
                }
            }
            assert(nxt == static_cast<int>(weights.size()));

            auto ans_dag = path_max_queries(n, forest, queries, [&](int a, int b) {
                return topo_rank[a] < topo_rank[b];
            });
            assert(ans_dag == ans);
        }
        return cmps;
    };

    // Random trees, paths, stars, caterpillars, disconnected forests.
    for (int trial = 0; trial < 40; ++trial) {
        int n = 20 + trial * 5;
        std::vector<ForestEdge> forest;
        std::vector<double> weights;
        std::uniform_real_distribution<double> wdist(-500.0, 500.0);

        for (int i = 1; i < n; ++i) {
            if (trial % 5 == 0 && (i % 7 == 0)) continue;  // Disconnected forest
            int p = (trial % 4 == 1) ? (i - 1) :             // Path
                    (trial % 4 == 2) ? 0 :                   // Star
                    (trial % 4 == 3) ? (i % 3 == 0 ? i - 1 : (i / 3) * 3) :  // Caterpillar-like
                                       static_cast<int>(rng() % static_cast<uint64_t>(i));
            int eid = static_cast<int>(forest.size());
            forest.push_back({p, i, eid});
            weights.push_back(wdist(rng));
        }

        std::vector<std::pair<int, int>> queries(300);
        for (auto& qr : queries) {
            qr = {static_cast<int>(rng() % static_cast<uint64_t>(n)),
                  static_cast<int>(rng() % static_cast<uint64_t>(n))};
        }
        run_forest_case(n, forest, weights, queries, /*check_dag=*/true);
    }

    // Large n=2000, q=10000 linear comparison check.
    {
        int n = 2000;
        std::vector<ForestEdge> forest;
        std::vector<double> weights;
        std::uniform_real_distribution<double> wdist(0.0, 1e6);
        for (int i = 1; i < n; ++i) {
            int p = static_cast<int>(rng() % static_cast<uint64_t>(i));
            int eid = static_cast<int>(forest.size());
            forest.push_back({p, i, eid});
            weights.push_back(wdist(rng));
        }
        std::vector<std::pair<int, int>> queries(10000);
        for (auto& qr : queries) {
            qr = {static_cast<int>(rng() % static_cast<uint64_t>(n)),
                  static_cast<int>(rng() % static_cast<uint64_t>(n))};
        }
        uint64_t cmps = run_forest_case(n, forest, weights, queries, /*check_dag=*/true);
        double ratio = static_cast<double>(cmps) / static_cast<double>(n + queries.size());
        assert(ratio < 6.0);
        std::cout << "  PASS (n=2000, q=10000 comparisons = " << cmps
                  << ", C/(n+q) = " << ratio << ", DAG dominance verified)\n";
    }
}

void test_batch_processor() {
    std::cout << "[3/6] Testing Batch Processor (Lemma 3.2 & Lemma 3.3)...\n";
    std::mt19937_64 rng(999);
    std::uniform_real_distribution<double> wdist(0.0, 100.0);

    // Lemma 3.3 test against brute-force over parallel groups with mixed admissibility.
    for (int trial = 0; trial < 30; ++trial) {
        int num_groups = 12;
        int m = 60;
        std::vector<double> weights(m);
        std::vector<bool> adm(m);
        std::vector<detail::IncidentEntry> items(m);
        for (int i = 0; i < m; ++i) {
            weights[i] = wdist(rng);
            adm[i] = (rng() % 3 != 0);
            items[i] = {i, i, static_cast<int>(rng() % static_cast<uint64_t>(num_groups))};
        }

        for (int k : {0, 1, 4, 8, 20}) {
            std::vector<detail::IncidentEntry> copy = items;
            DLMSTStats st;
            WeightOracle oracle(weights, &st);
            auto got = select_k_representatives(
                num_groups, m, copy, k, oracle, [&](int idx) { return adm[idx]; });

            // Brute force reference.
            std::vector<detail::IncidentEntry> best_per_group(num_groups, {-1, -1, -1});
            for (const auto& it : items) {
                if (!adm[it.edge_idx]) continue;
                auto& b = best_per_group[it.group];
                if (b.edge_idx == -1 || oracle.less_unchecked(it.orig_id, b.orig_id)) {
                    b = it;
                }
            }
            std::vector<detail::IncidentEntry> expected;
            for (const auto& b : best_per_group) {
                if (b.edge_idx != -1) expected.push_back(b);
            }
            std::sort(expected.begin(), expected.end(), [&](const auto& a, const auto& b) {
                return oracle.less_unchecked(a.orig_id, b.orig_id);
            });
            if (static_cast<int>(expected.size()) > k) expected.resize(k);

            std::vector<int> got_ids, exp_ids;
            for (const auto& x : got) got_ids.push_back(x.orig_id);
            for (const auto& x : expected) exp_ids.push_back(x.orig_id);
            std::sort(got_ids.begin(), got_ids.end());
            std::sort(exp_ids.begin(), exp_ids.end());
            assert(got_ids == exp_ids);
        }
    }

    // Lemma 3.2 full batch test: |P'| <= 2n/r, |D| <= k_rep * n, and invariant holds.
    {
        Graph g = make_random_graph(80, 1800, 2026, /*allow_parallel=*/true);
        std::vector<CertEdge> cedges(g.edges.size());
        std::vector<int> batch_idx(g.edges.size());
        for (size_t i = 0; i < g.edges.size(); ++i) {
            cedges[i] = {g.edges[i].u, g.edges[i].v, g.edges[i].orig_id, -1, {}};
            batch_idx[i] = static_cast<int>(i);
        }
        int r = 5;
        int64_t k_rep = 25;
        std::vector<int> sigma(g.n, -1), pi(g.n, -1);
        std::vector<bool> in_P(g.n, false);
        DLMSTStats st;
        WeightOracle oracle(g.weights, &st);
        auto bres = process_batch(g.n, cedges, batch_idx, r, k_rep, sigma, pi, in_P, oracle);

        assert(static_cast<int>(bres.new_portals.size()) <= (2 * g.n) / r);
        assert(static_cast<int64_t>(bres.declared_edges.size()) <= k_rep * g.n);

        // Compute F = MSF(D) using Kruskal and verify Lemma 3.2 on all of B \ R'.
        Graph dg;
        dg.n = g.n;
        dg.weights = g.weights;
        for (int idx : bres.declared_edges) dg.edges.push_back(g.edges[idx]);
        auto msf_d = kruskal_msf(dg);
        std::vector<ForestEdge> f_edges;
        for (int oid : msf_d) f_edges.push_back({g.edges[oid].u, g.edges[oid].v, oid});

        std::vector<std::pair<int, int>> queries;
        for (const auto& ce : cedges) queries.push_back({ce.u, ce.v});
        auto fmax = path_max_queries(
            g.n, f_edges, queries, [&](int a, int b) { return oracle.less_unchecked(a, b); });

        for (size_t i = 0; i < cedges.size(); ++i) {
            const auto& ce = cedges[i];
            bool supported = (fmax[i] != -1) && !oracle.less_unchecked(ce.orig_id, fmax[i]);
            bool valid_proj = (sigma[ce.u] != -1 && sigma[ce.v] != -1 &&
                               !oracle.less_unchecked(ce.orig_id, cedges[sigma[ce.u]].orig_id) &&
                               !oracle.less_unchecked(ce.orig_id, cedges[sigma[ce.v]].orig_id));
            assert(supported || valid_proj);
        }
        std::cout << "  PASS (|P'|=" << bres.new_portals.size() << " <= " << (2 * g.n) / r
                  << ", |D|=" << bres.declared_edges.size() << ", Lemma 3.2 & 3.3 verified)\n";
    }
}

void test_small_schedule_lemma_3_11() {
    std::cout << "[4/6] Testing constructive Lemma 3.11 exact schedule search on K4...\n";
    Graph k4 = make_complete_graph(4, 1);
    auto res = find_exact_schedule_small(k4.n, k4.edges, /*t=*/2);
    assert(res.num_schedules_checked == 4096);
    assert(res.num_weight_orders_checked == 720);
    assert(res.min_worst_bad_edges < static_cast<int>(k4.edges.size()));
    std::cout << "  PASS (checked " << res.num_schedules_checked << " schedules x "
              << res.num_weight_orders_checked << " permutations; min worst bad edges = "
              << res.min_worst_bad_edges << "/" << k4.edges.size() << ")\n";
}

void test_action_estimator_dag() {
    std::cout << "[5/6] Testing Action Estimator (DAGTopological sub-sampling)...\n";
    Graph g = make_random_graph(180, 3600, 5555);
    std::vector<CertEdge> cedges(g.edges.size());
    for (size_t i = 0; i < g.edges.size(); ++i) {
        cedges[i] = {g.edges[i].u, g.edges[i].v, g.edges[i].orig_id, -1, {}};
    }

    // Build a spanning forest F via Kruskal on first 400 edges and assign spokes.
    Graph sub;
    sub.n = g.n;
    sub.weights = g.weights;
    sub.edges.assign(g.edges.begin(), g.edges.begin() + 400);
    auto f_oids = kruskal_msf(sub);
    std::vector<int> F = f_oids;
    std::vector<int> sigma(g.n, -1), pi(g.n, -1);
    for (int f_oid : F) {
        const auto& ce = cedges[f_oid];
        if (sigma[ce.u] == -1) { sigma[ce.u] = f_oid; pi[ce.u] = ce.v % 10; }
        if (sigma[ce.v] == -1) { sigma[ce.v] = f_oid; pi[ce.v] = ce.u % 10; }
    }

    std::vector<int> A_i;
    for (int i = 400; i < static_cast<int>(cedges.size()); ++i) A_i.push_back(i);

    DLMSTConfig cfg_exact;
    cfg_exact.estimator_mode = EstimatorMode::Exact;
    DLMSTStats st_exact;
    WeightOracle oracle_exact(g.weights, &st_exact);
    auto act_exact = estimate_preferred_actions(
        g.n, cedges, F, sigma, pi, A_i, 3, 1234ULL, cfg_exact, oracle_exact, &st_exact);

    DLMSTConfig cfg_dag;
    cfg_dag.estimator_mode = EstimatorMode::DAGTopological;
    cfg_dag.max_sample_fraction = 0.25;
    cfg_dag.verify_invariants = true;
    DLMSTStats st_dag;
    WeightOracle oracle_dag(g.weights, &st_dag);
    auto act_dag = estimate_preferred_actions(
        g.n, cedges, F, sigma, pi, A_i, 3, 1234ULL, cfg_dag, oracle_dag, &st_dag);

    assert(st_dag.estimator_dag_invocations == 1);
    assert(st_dag.cmp_estimator < st_exact.cmp_estimator);
    int agree = 0;
    for (size_t i = 0; i < A_i.size(); ++i) {
        if (act_dag[i] == act_exact[i]) ++agree;
    }
    double acc = static_cast<double>(agree) / static_cast<double>(A_i.size());
    assert(acc >= 0.70);
    std::cout << "  PASS (DAG queried " << st_dag.estimator_queried_edges << "/"
              << st_dag.estimator_ambig_edges << " ambig edges, saved "
              << (st_exact.cmp_estimator - st_dag.cmp_estimator)
              << " cmps, full-set agreement = " << (acc * 100.0) << "%)\n";
}

void test_end_to_end_dlmst() {
    std::cout << "[6/6] Testing End-to-End DLMST across edge cases, families & configs...\n";

    DLMSTConfig cfg_practical;
    cfg_practical.boruvka_steps = 4;
    cfg_practical.param_schedule = ParamSchedule::Practical;
    cfg_practical.verify_invariants = true;

    // Edge cases.
    {
        Graph empty{0, {}, {}};
        assert_matches_kruskal(empty, cfg_practical, "empty");

        Graph isolated{10, {}, {}};
        assert_matches_kruskal(isolated, cfg_practical, "isolated_only");

        Graph single_v{1, {}, {}};
        assert_matches_kruskal(single_v, cfg_practical, "single_vertex");

        Graph single_e{2, {{0, 1, 0}}, {3.14}};
        assert_matches_kruskal(single_e, cfg_practical, "single_edge");

        Graph loops_only{5, {{0, 0, 0}, {2, 2, 1}, {4, 4, 2}}, {1.0, 2.0, 3.0}};
        assert_matches_kruskal(loops_only, cfg_practical, "loops_only");

        Graph equal_w = make_random_graph(300, 2000, 111, /*allow_parallel=*/true, /*equal_weights=*/true);
        assert_matches_kruskal(equal_w, cfg_practical, "all_equal_weights");

        // Disconnected multi-component graph with ascending/descending weights.
        Graph disconn;
        disconn.n = 250;
        for (int c = 0; c < 5; ++c) {
            int base = c * 50;
            for (int u = 0; u < 45; ++u) {
                for (int v = u + 1; v < 45; v += 3) {
                    int oid = static_cast<int>(disconn.edges.size());
                    disconn.edges.push_back({base + u, base + v, oid});
                    disconn.weights.push_back(
                        (c % 2 == 0) ? static_cast<double>(oid) : -static_cast<double>(oid));
                }
            }
        }
        assert_matches_kruskal(disconn, cfg_practical, "disconnected_adversarial");
    }

    // Graph families across Practical and PaperAsymptotic.
    {
        Graph er_sparse = make_random_graph(1200, 4800, 2001);
        Graph er_dense = make_random_graph(400, 8000, 2002);
        Graph complete = make_complete_graph(80, 2003);
        Graph grid = make_grid_graph(30, 30, 2004);
        Graph ba_hub = make_ba_hub_graph(800, 6, 2005);
        Graph multi = make_random_graph(300, 6000, 2006, /*allow_parallel=*/true);

        DLMSTConfig cfg_paper;
        cfg_paper.boruvka_steps = 4;
        cfg_paper.param_schedule = ParamSchedule::PaperAsymptotic;
        cfg_paper.gamma = 1;
        cfg_paper.verify_invariants = true;

        DLMSTConfig cfg_exact = cfg_practical;
        cfg_exact.estimator_mode = EstimatorMode::Exact;

        for (const auto& [name, g] : std::vector<std::pair<const char*, Graph>>{
                 {"ER_sparse", er_sparse},
                 {"ER_dense", er_dense},
                 {"Complete_K80", complete},
                 {"Grid_30x30", grid},
                 {"BA_hub", ba_hub},
                 {"Multigraph", multi}}) {
            assert_matches_kruskal(g, cfg_practical, name);
            assert_matches_kruskal(g, cfg_paper, name);
            assert_matches_kruskal(g, cfg_exact, name);
        }
    }

    // DeepStress test explicitly proving multi-level portal recursion, projections,
    // admissibility rejections, and DAG estimator sub-sampling all fire and stay exact.
    {
        Graph stress = make_random_graph(140, 4800, 314159, /*allow_parallel=*/true);
        DLMSTConfig cfg_deep;
        cfg_deep.boruvka_steps = 1;
        cfg_deep.param_schedule = ParamSchedule::DeepStress;
        cfg_deep.estimator_mode = EstimatorMode::DAGTopological;
        cfg_deep.max_sample_fraction = 0.35;
        cfg_deep.verify_invariants = true;

        DLMSTStats st;
        assert_matches_kruskal(stress, cfg_deep, "DeepStress_multigraph", &st);
        std::cout << "  DeepStress stats: max_cert_depth=" << st.max_cert_depth
                  << ", portals=" << st.total_portals
                  << ", projected=" << st.total_projected_edges
                  << ", rejected=" << st.total_rejected_edges
                  << ", dag_invocations=" << st.estimator_dag_invocations << "\n";
        assert(st.max_cert_depth >= 2);
        assert(st.total_portals > 0);
        assert(st.total_projected_edges > 0);
        assert(st.total_rejected_edges > 0);
        assert(st.estimator_dag_invocations > 0);
    }

    std::cout << "  PASS (all end-to-end tests matched Kruskal MSF 100%)\n";
}

}  // namespace

int main() {
    test_bfprt();
    test_path_max();
    test_batch_processor();
    test_small_schedule_lemma_3_11();
    test_action_estimator_dag();
    test_end_to_end_dlmst();
    std::cout << "\nALL DLMST TESTS PASSED.\n";
    return 0;
}
