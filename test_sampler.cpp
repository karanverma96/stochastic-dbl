#include "uncertain_graph.hpp"
#include "monte_carlo_sampler.hpp"
#include "hoeffding.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

// Exact reachability probability via possible-world enumeration.
// Only tractable for small edge counts (2^|E| possible worlds).
double exactReachability(const UncertainGraph& ug, int s, int t) {
    auto& edges = ug.edges();
    int m = static_cast<int>(edges.size());
    int n = ug.numVertices();
    double total = 0.0;

    for (uint64_t mask = 0; mask < (1ULL << m); ++mask) {
        Graph g(n);
        double worldProb = 1.0;
        for (int i = 0; i < m; ++i) {
            auto& [u, v, p] = edges[i];
            bool present = (mask >> i) & 1ULL;
            worldProb *= present ? p : (1.0 - p);
            if (present) g.addEdge(u, v);
        }
        if (worldProb > 0.0 && bfsReachable(g, s, t)) {
            total += worldProb;
        }
    }
    return total;
}

int main() {
    std::cout << std::fixed << std::setprecision(4);

    // --- Sanity check: Hoeffding sample size behaves as expected ---
    std::cout << "Hoeffding sample sizes:\n";
    for (double eps : {0.05, 0.01, 0.005}) {
        uint64_t n = hoeffdingSampleSize(eps, 0.05);
        std::cout << "  epsilon=" << eps << " delta=0.05 -> N=" << n << "\n";
    }
    std::cout << "\n";

    // --- Build a small uncertain graph: 0->1->3, 0->2->3, 1->2 ---
    // Deterministic edge list (small enough for exact enumeration: 2^6 = 64 worlds)
    std::vector<std::pair<int,int>> rawEdges = {
        {0, 1}, {1, 3}, {0, 2}, {2, 3}, {1, 2}, {2, 1}
    };
    int n = 4;
    std::mt19937 rng(7);

    for (auto scheme : {ProbabilityScheme::Uniform, ProbabilityScheme::Trivalency, ProbabilityScheme::WeightedCascade}) {
        std::string name = scheme == ProbabilityScheme::Uniform ? "Uniform"
                          : scheme == ProbabilityScheme::Trivalency ? "Trivalency"
                          : "WeightedCascade";
        UncertainGraph ug = UncertainGraph::fromEdgeList(n, rawEdges, scheme, rng, /*uniformP=*/0.3);

        double exact = exactReachability(ug, 0, 3);

        MonteCarloSampler sampler(ug, /*seed=*/99);
        auto est = sampler.estimate(0, 3, /*epsilon=*/0.01, /*delta=*/0.05);

        double diff = std::abs(exact - est.probability);
        std::cout << "[" << name << "] exact=" << exact
                  << "  monte_carlo=" << est.probability
                  << "  samples=" << est.samples
                  << "  |diff|=" << diff
                  << (diff <= est.epsilon * 3 ? "  OK" : "  ** OUT OF TOLERANCE **")
                  << "\n";
    }

    // --- Batch estimation demo across multiple query pairs on the same worlds ---
    std::cout << "\nBatch estimation (Weighted Cascade, multiple query pairs):\n";
    UncertainGraph ugBatch = UncertainGraph::fromEdgeList(n, rawEdges, ProbabilityScheme::WeightedCascade, rng);
    MonteCarloSampler batchSampler(ugBatch, 555);
    std::vector<std::pair<int,int>> queries = {{0,3},{0,2},{1,3},{2,1}};
    auto results = batchSampler.batchEstimate(queries, 0.02, 0.05);
    for (size_t i = 0; i < queries.size(); ++i) {
        std::cout << "  P(" << queries[i].first << " -> " << queries[i].second
                  << ") ~= " << results[i] << "\n";
    }

    return 0;
}
