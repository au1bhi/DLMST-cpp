#include "dlmst/boruvka.hpp"
#include "dlmst/dlmst.hpp"
#include "dlmst/kruskal.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace dlmst;

namespace {

Graph gen_graph(int n, int m, uint64_t seed, bool allow_parallel) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> wdist(0.0, 1e6);
    Graph g;
    g.n = n;
    g.edges.reserve(m);
    g.weights.reserve(m);

    std::set<std::pair<int, int>> seen;
    // Connected spanning tree backbone.
    for (int i = 1; i < n && static_cast<int>(g.edges.size()) < m; ++i) {
        int p = static_cast<int>(rng() % static_cast<uint64_t>(i));
        if (!allow_parallel) {
            seen.insert({std::min(i, p), std::max(i, p)});
        }
        int oid = static_cast<int>(g.edges.size());
        g.edges.push_back({p, i, oid});
        g.weights.push_back(wdist(rng));
    }

    while (static_cast<int>(g.edges.size()) < m) {
        int u = static_cast<int>(rng() % static_cast<uint64_t>(n));
        int v = static_cast<int>(rng() % static_cast<uint64_t>(n));
        if (u == v) continue;
        if (!allow_parallel) {
            auto key = std::make_pair(std::min(u, v), std::max(u, v));
            if (seen.count(key)) continue;
            seen.insert(key);
        }
        int oid = static_cast<int>(g.edges.size());
        g.edges.push_back({u, v, oid});
        g.weights.push_back(wdist(rng));
    }
    return g;
}

uint64_t pure_boruvka_cmps(const Graph& g, std::vector<int>& out_ids) {
    DLMSTStats st;
    WeightOracle oracle(g.weights, &st);
    auto bres = run_boruvka_phases(g.n, g.edges, 64, oracle);
    out_ids = std::move(bres.contracted_orig_ids);
    std::sort(out_ids.begin(), out_ids.end());
    return st.cmp_boruvka;
}

struct BenchCase {
    std::string family;
    int n;
    int m;
    bool parallel;
    DLMSTConfig cfg;
};

}  // namespace

int main() {
    std::cout << "=== DLMST C++20 Reproduction Benchmark & Scaling Verification ===\n";
    std::cout << "Paper: Khanna, Putterman, Song (2026) - https://github.com/Louiealbp/DLMST\n\n";

    DLMSTConfig cfg_prac;
    cfg_prac.boruvka_steps = 4;
    cfg_prac.param_schedule = ParamSchedule::Practical;
    cfg_prac.estimator_mode = EstimatorMode::DAGTopological;

    DLMSTConfig cfg_dense_rec;
    cfg_dense_rec.boruvka_steps = 2;
    cfg_dense_rec.param_schedule = ParamSchedule::DeepStress;
    cfg_dense_rec.estimator_mode = EstimatorMode::DAGTopological;

    std::vector<BenchCase> cases = {
        {"Sparse (m=4n)",    2000,   8000, false, cfg_prac},
        {"Sparse (m=4n)",   10000,  40000, false, cfg_prac},
        {"Sparse (m=4n)",   50000, 200000, false, cfg_prac},
        {"Medium (m=16n)",   2000,  32000, false, cfg_prac},
        {"Medium (m=16n)",  10000, 160000, false, cfg_prac},
        {"Dense (m=50n)",    1000,  50000, false, cfg_dense_rec},
        {"Dense (m=50n)",    3000, 150000, false, cfg_dense_rec},
        {"Multi (m=40n)",    2000,  80000, true,  cfg_dense_rec},
    };

    std::cout << std::left
              << std::setw(16) << "Family"
              << std::right
              << std::setw(8) << "n"
              << std::setw(9) << "m"
              << std::setw(11) << "DLMST_cmp"
              << std::setw(8) << "C/m"
              << std::setw(11) << "Krusk_cmp"
              << std::setw(8) << "K/m"
              << std::setw(11) << "Boruv_cmp"
              << std::setw(8) << "B/m"
              << std::setw(9) << "Time(ms)"
              << std::setw(7) << "Depth"
              << std::setw(9) << "Portals"
              << std::setw(9) << "Reject"
              << "  Status\n";
    std::cout << std::string(126, '-') << "\n";

    uint64_t seed = 424242;
    for (const auto& bc : cases) {
        Graph g = gen_graph(bc.n, bc.m, seed++, bc.parallel);

        uint64_t k_cmp = 0;
        auto expected = kruskal_msf(g, &k_cmp);

        std::vector<int> b_ids;
        uint64_t b_cmp = pure_boruvka_cmps(g, b_ids);
        assert(b_ids == expected);

        auto t0 = std::chrono::steady_clock::now();
        auto res = dlmst_msf(g, bc.cfg);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        bool ok = (res.edge_ids == expected);
        assert(ok);

        const auto& st = res.stats;
        const double m_d = static_cast<double>(g.edges.size());
        std::cout << std::left << std::setw(16) << bc.family
                  << std::right
                  << std::setw(8) << g.n
                  << std::setw(9) << g.edges.size()
                  << std::setw(11) << st.total_comparisons()
                  << std::setw(8) << std::fixed << std::setprecision(2)
                  << (static_cast<double>(st.total_comparisons()) / m_d)
                  << std::setw(11) << k_cmp
                  << std::setw(8) << (static_cast<double>(k_cmp) / m_d)
                  << std::setw(11) << b_cmp
                  << std::setw(8) << (static_cast<double>(b_cmp) / m_d)
                  << std::setw(9) << std::setprecision(1) << ms
                  << std::setw(7) << st.max_cert_depth
                  << std::setw(9) << st.total_portals
                  << std::setw(9) << st.total_rejected_edges
                  << (ok ? "  PASS" : "  FAIL") << "\n";

        std::cout << "    Breakdown: Boruvka=" << st.cmp_boruvka
                  << " BatchSel=" << st.cmp_batch_select
                  << " Spoke=" << st.cmp_batch_spokes
                  << " AdmPass=" << st.cmp_admissibility_pass
                  << " AdmFail=" << st.cmp_admissibility_fail
                  << " Est=" << st.cmp_estimator
                  << " Filter=" << st.cmp_filter_verify
                  << " | RootFilter: " << st.root_edges_before_filter
                  << " -> " << st.root_edges_after_filter
                  << " | Proj=" << st.total_projected_edges
                  << " DAG_inv=" << st.estimator_dag_invocations << "\n";
    }

    return 0;
}
