#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace dlmst {

struct Edge {
    int u = -1;
    int v = -1;
    int orig_id = -1;
};

struct Graph {
    int n = 0;
    std::vector<Edge> edges;
    std::vector<double> weights;  // Indexed by orig_id in [0, weights.size())
};

enum class CmpCategory {
    Boruvka,
    BatchSelect,
    BatchSpoke,
    AdmissibilityPass,
    AdmissibilityFail,
    Estimator,
    FilterVerify,
    ExactReference
};

struct DLMSTStats {
    uint64_t cmp_boruvka = 0;
    uint64_t cmp_batch_select = 0;
    uint64_t cmp_batch_spokes = 0;
    uint64_t cmp_admissibility_pass = 0;
    uint64_t cmp_admissibility_fail = 0;
    uint64_t cmp_estimator = 0;
    uint64_t cmp_filter_verify = 0;
    uint64_t cmp_exact_ref = 0;

    uint64_t msf_calls = 0;
    uint64_t cert_calls = 0;
    int max_msf_depth = 0;
    int max_cert_depth = 0;

    uint64_t total_portals = 0;
    uint64_t total_projected_edges = 0;
    uint64_t total_rejected_edges = 0;
    uint64_t total_batch_edges = 0;
    uint64_t total_declared_edges = 0;

    uint64_t estimator_ambig_edges = 0;
    uint64_t estimator_queried_edges = 0;
    uint64_t estimator_dag_invocations = 0;

    uint64_t root_edges_before_filter = 0;
    uint64_t root_edges_after_filter = 0;

    uint64_t total_comparisons() const {
        return cmp_boruvka + cmp_batch_select + cmp_batch_spokes +
               cmp_admissibility_pass + cmp_admissibility_fail +
               cmp_estimator + cmp_filter_verify;
    }

    void add_cmp(CmpCategory cat, uint64_t delta = 1) {
        switch (cat) {
            case CmpCategory::Boruvka:           cmp_boruvka += delta; break;
            case CmpCategory::BatchSelect:       cmp_batch_select += delta; break;
            case CmpCategory::BatchSpoke:        cmp_batch_spokes += delta; break;
            case CmpCategory::AdmissibilityPass: cmp_admissibility_pass += delta; break;
            case CmpCategory::AdmissibilityFail: cmp_admissibility_fail += delta; break;
            case CmpCategory::Estimator:         cmp_estimator += delta; break;
            case CmpCategory::FilterVerify:      cmp_filter_verify += delta; break;
            case CmpCategory::ExactReference:    cmp_exact_ref += delta; break;
        }
    }
};

/**
 * Opaque comparison oracle implementing strict edge weight ordering
 * with deterministic tie-breaking by original edge index.
 */
class WeightOracle {
public:
    WeightOracle(const std::vector<double>& weights, DLMSTStats* stats = nullptr)
        : weights_(weights), stats_(stats) {}

    bool less(int e1, int e2, CmpCategory cat) {
        if (e1 == e2) return false;
        if (stats_) stats_->add_cmp(cat, 1);
        return less_unchecked(e1, e2);
    }

    bool less_unchecked(int e1, int e2) const {
        if (e1 == e2) return false;
        const double w1 = weights_[e1];
        const double w2 = weights_[e2];
        if (w1 != w2) return w1 < w2;
        return e1 < e2;
    }

    void record_cmp(CmpCategory cat, uint64_t count = 1) {
        if (stats_) stats_->add_cmp(cat, count);
    }

    const std::vector<double>& weights() const { return weights_; }
    DLMSTStats* stats() const { return stats_; }

private:
    const std::vector<double>& weights_;
    DLMSTStats* stats_;
};

struct ProjectionStep {
    int level = 0;
    int spoke_u = -1;  // orig_id of spoke at u
    int spoke_v = -1;  // orig_id of spoke at v
};

struct CertEdge {
    int u = -1;
    int v = -1;
    int orig_id = -1;
    int parent_edge_idx = -1;
    std::vector<ProjectionStep> proj_record;
};

enum class ParamSchedule {
    PaperAsymptotic,  // t = max(2, gamma * d^3), r = 8 * t^3, k = r^2
    Practical,        // t = 2 + d, r = 2*t + 2, k = r^2
    DeepStress        // t = 2, r = 5, k = 25 (triggers multi-level projections even on small graphs)
};

enum class EstimatorMode {
    DAGTopological,
    Exact
};

struct DLMSTConfig {
    int boruvka_steps = 6;
    ParamSchedule param_schedule = ParamSchedule::Practical;
    int gamma = 1;
    EstimatorMode estimator_mode = EstimatorMode::DAGTopological;
    double estimator_sample_c = 1.5;
    double max_sample_fraction = 1.0;
    uint64_t schedule_seed = 0x9e3779b97f4a7c15ULL;
    bool verify_invariants = false;

    void compute_level_params(int depth, int& t, int& r, int64_t& k_rep) const {
        const int d = std::max(1, depth);
        switch (param_schedule) {
            case ParamSchedule::PaperAsymptotic: {
                const int64_t g = std::max(1, gamma);
                const int64_t d3 = static_cast<int64_t>(d) * d * d;
                const int64_t t64 = std::min<int64_t>(1000, std::max<int64_t>(2, g * d3));
                t = static_cast<int>(t64);
                const int64_t r64 = std::min<int64_t>(1000000, 8LL * t64 * t64 * t64);
                r = static_cast<int>(std::max<int64_t>(2LL * t + 1, r64));
                k_rep = static_cast<int64_t>(r) * r;
                break;
            }
            case ParamSchedule::Practical: {
                t = 2 + d;
                r = 2 * t + 2;
                k_rep = static_cast<int64_t>(r) * r;
                break;
            }
            case ParamSchedule::DeepStress: {
                t = 2;
                r = 5;
                k_rep = 25;
                break;
            }
        }
    }
};

struct MSFResult {
    std::vector<int> edge_ids;  // Sorted orig_ids in the minimum spanning forest
    DLMSTStats stats;
};

/**
 * Deterministic SplitMix64 mixing function (no weight comparisons).
 */
inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

/**
 * Simple Disjoint Set Union with path compression and union by rank.
 */
struct DSU {
    std::vector<int> parent;
    std::vector<int> rank;

    explicit DSU(int n = 0) { reset(n); }

    void reset(int n) {
        parent.resize(n);
        rank.assign(n, 0);
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int find(int x) {
        int root = x;
        while (parent[root] != root) root = parent[root];
        while (parent[x] != root) {
            int nxt = parent[x];
            parent[x] = root;
            x = nxt;
        }
        return root;
    }

    bool unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
        return true;
    }
};

}  // namespace dlmst
