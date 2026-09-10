#pragma once
#include "uncertain_graph.hpp"
#include "dbl_index.hpp"
#include "hoeffding.hpp"
#include <random>
#include <queue>
#include <vector>

// Plain BFS on an already-materialised world. Only batchEstimate() needs this,
// because it indexes a real graph; estimate() never builds a world at all.
inline bool bfsReachable(const Graph& g, int s, int t) {
    if (s == t) return true;
    std::vector<char> visited(g.numVertices(), 0);
    std::queue<int> q;
    q.push(s);
    visited[s] = 1;
    while (!q.empty()) {
        int p = q.front(); q.pop();
        for (int x : g.Suc(p)) {
            if (x == t) return true;
            if (!visited[x]) { visited[x] = 1; q.push(x); }
        }
    }
    return false;
}

struct ReachabilityEstimate {
    double probability;
    uint64_t samples;
    double epsilon; 
    double delta;
};

class MonteCarloSampler {
public:
    explicit MonteCarloSampler(const UncertainGraph& ug, uint64_t seed = 1234)
        : ug_(ug), rng_(seed), adj_(ug.numVertices()), visited_(ug.numVertices(), 0) {
        for (auto& [u, v, p] : ug.edges()) adj_[u].emplace_back(v, p);
    }

    // Draws N worlds (N from Hoeffding) and counts how many connect s to t
    ReachabilityEstimate estimate(int s, int t, double epsilon = 0.01, double delta = 0.05) {
        uint64_t n = hoeffdingSampleSize(epsilon, delta);
        uint64_t successes = 0;
        for (uint64_t i = 0; i < n; ++i) {
            if (lazyReachable(s, t)) ++successes;
        }
        return { static_cast<double>(successes) / static_cast<double>(n), n, epsilon, delta };
    }

    // Many pairs against the same worlds: build one DBL index per world and
    // query it, instead of a BFS per (world, pair).
    //
    // Worth it only in bulk. Indexing a world costs 2k + k' traversals (~64 at
    // the sizes below) against one BFS per pair, so it does not pay until the
    // query list runs into the thousands.
    std::vector<double> batchEstimate(
        const std::vector<std::pair<int,int>>& queries,
        double epsilon = 0.01,
        double delta = 0.05) {

        uint64_t n = hoeffdingSampleSize(epsilon, delta);
        std::vector<uint64_t> successCounts(queries.size(), 0);

        for (uint64_t i = 0; i < n; ++i) {
            Graph world = ug_.sampleWorld(rng_);
            DBLIndex idx(world, /*k=*/16, /*kp=*/32);
            idx.build();
            for (size_t qi = 0; qi < queries.size(); ++qi) {
                if (idx.query(queries[qi].first, queries[qi].second)) {
                    ++successCounts[qi];
                }
            }
        }

        std::vector<double> results(queries.size());
        for (size_t qi = 0; qi < queries.size(); ++qi) {
            results[qi] = static_cast<double>(successCounts[qi]) / static_cast<double>(n);
        }
        return results;
    }

private:
    // One world's s-t check, flipping an edge's coin only when the search
    // arrives at it. Each vertex is dequeued once, so no edge is drawn twice,
    // and an edge never reached cannot sit on an s-t path -- so this comes out
    // distributed exactly like sampling the whole world and then walking it.
    // On wiki-Vote that skips about 99% of the coin flips.
    bool lazyReachable(int s, int t) {
        if (s == t) return true;
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        std::queue<int> q;
        q.push(s);
        visited_[s] = 1;
        touched_.push_back(s);
        bool found = false;
        while (!q.empty() && !found) {
            int p = q.front(); q.pop();
            for (auto& [x, prob] : adj_[p]) {
                // Reached already, so this edge cannot change the answer. Needs
                // no special case for t, which returns below before it is ever
                // marked visited.
                if (visited_[x]) continue;
                if (unit(rng_) < prob) {
                    if (x == t) { found = true; break; }
                    visited_[x] = 1;
                    touched_.push_back(x);
                    q.push(x);
                }
            }
        }
        // Clearing by touched list, not by refilling the whole array -- a
        // sample usually reaches a handful of vertices out of thousands
        for (int v : touched_) visited_[v] = 0;
        touched_.clear();
        return found;
    }

    const UncertainGraph& ug_;
    std::mt19937 rng_;
    std::vector<std::vector<std::pair<int,double>>> adj_; // (target, probability)
    std::vector<char> visited_;                           // reused across samples
    std::vector<int> touched_;
};
