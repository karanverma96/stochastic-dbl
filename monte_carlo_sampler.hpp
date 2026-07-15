#pragma once
#include "uncertain_graph.hpp"
#include "dbl_index.hpp"
#include "hoeffding.hpp"
#include <random>
#include <queue>
#include <vector>

// Plain BFS reachability check on a single sampled possible-world graph.
// Used as the per-sample check: since each Monte Carlo sample produces a
// *different* graph that is only queried once (for the single s-t pair of
// interest), a one-off BFS is cheaper than building a full DBL index for
// each sampled world. DBL is the better tool when many queries are run
// against the *same* graph snapshot -- see batchEstimate() below, which
// uses exactly that scenario.
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
        : ug_(ug), rng_(seed) {}

    // Single-query estimator: draws N possible worlds (N from Hoeffding's
    // inequality) and checks s-t reachability once per world via BFS.
    ReachabilityEstimate estimate(int s, int t, double epsilon = 0.01, double delta = 0.05) {
        uint64_t n = hoeffdingSampleSize(epsilon, delta);
        uint64_t successes = 0;
        for (uint64_t i = 0; i < n; ++i) {
            Graph world = ug_.sampleWorld(rng_);
            if (bfsReachable(world, s, t)) ++successes;
        }
        return { static_cast<double>(successes) / static_cast<double>(n), n, epsilon, delta };
    }

    // Batch estimator: when many (s, t) pairs need probabilities against
    // the *same* set of sampled worlds, it pays off to build a DBL index
    // per world once and answer every query against it, rather than
    // re-running BFS per (world, query) pair.
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
    const UncertainGraph& ug_;
    std::mt19937 rng_;
};
