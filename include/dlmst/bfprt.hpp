#pragma once

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace dlmst {

namespace detail {

template <typename T, typename Cmp>
void insertion_sort_range(std::vector<T>& a, size_t lo, size_t hi, Cmp& cmp) {
    for (size_t i = lo + 1; i < hi; ++i) {
        T key = std::move(a[i]);
        size_t j = i;
        while (j > lo && cmp(key, a[j - 1])) {
            a[j] = std::move(a[j - 1]);
            --j;
        }
        a[j] = std::move(key);
    }
}

/**
 * Partition a[lo..hi) around pivot initially stored at pivot_idx,
 * returning its final sorted position p.
 * Assumes elements have a strict total order (no ties under cmp).
 */
template <typename T, typename Cmp>
size_t partition_around(std::vector<T>& a, size_t lo, size_t hi, size_t pivot_idx, Cmp& cmp) {
    std::swap(a[pivot_idx], a[hi - 1]);
    size_t store = lo;
    for (size_t i = lo; i + 1 < hi; ++i) {
        if (cmp(a[i], a[hi - 1])) {
            std::swap(a[i], a[store]);
            ++store;
        }
    }
    std::swap(a[store], a[hi - 1]);
    return store;
}

}  // namespace detail

/**
 * Deterministic Blum-Floyd-Pratt-Rivest-Tarjan (BFPRT) median-of-medians
 * selection in O(m) worst-case comparisons and zero dynamic allocations.
 *
 * Places the nth element (absolute index in [lo, hi)) in its sorted position,
 * with all elements in [lo, nth) strictly smaller and [nth + 1, hi) strictly larger.
 */
template <typename T, typename Cmp>
void bfprt_nth_element(std::vector<T>& a, size_t lo, size_t hi, size_t nth, Cmp cmp) {
    while (true) {
        if (hi - lo <= 10) {
            detail::insertion_sort_range(a, lo, hi, cmp);
            return;
        }

        const size_t n = hi - lo;
        const size_t num_chunks = (n + 4) / 5;
        for (size_t c = 0; c < num_chunks; ++c) {
            const size_t sub_lo = lo + c * 5;
            const size_t sub_hi = std::min(hi, sub_lo + 5);
            detail::insertion_sort_range(a, sub_lo, sub_hi, cmp);
            const size_t median_idx = sub_lo + (sub_hi - sub_lo - 1) / 2;
            std::swap(a[lo + c], a[median_idx]);
        }

        const size_t mom_idx = lo + num_chunks / 2;
        bfprt_nth_element(a, lo, lo + num_chunks, mom_idx, cmp);

        const size_t p = detail::partition_around(a, lo, hi, mom_idx, cmp);
        if (nth == p) {
            return;
        } else if (nth < p) {
            hi = p;
        } else {
            lo = p + 1;
        }
    }
}

}  // namespace dlmst
