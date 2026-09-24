#include "uncertain_graph.hpp"
#include "monte_carlo_sampler.hpp"
#include "hoeffding.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

// Exact answer by walking all 2^|E| possible worlds. The thing Monte Carlo is
// approximating, and only affordable on a toy graph.
double exactReachability(const UncertainGraph& ug, int s, int t) {
    auto& edges = ug.edges();
    int m = static_cast<int>(edges.size());
    int n = ug.numVertices();
    // Past ~25 edges this is hopeless anyway, and at m >= 64 the shift below
    // is undefined behaviour
    if (m > 25) {
        throw std::invalid_argument("exactReachability: too many edges to enumerate");
    }
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
    int failures = 0;

    // --- Sanity check: Hoeffding sample size behaves as expected ---
    std::cout << "Hoeffding sample sizes:\n";
    for (double eps : {0.05, 0.01, 0.005}) {
        uint64_t n = hoeffdingSampleSize(eps, 0.05);
        std::cout << "  epsilon=" << eps << " delta=0.05 -> N=" << n << "\n";
    }
    std::cout << "\n";

    // --- Rejected arguments, NaN included ---
    // hoeffdingSampleSize() ends in static_cast<uint64_t>(std::ceil(n)), which
    // is undefined behaviour for anything not representable. A NaN argument
    // reaches it unless the guards are written as a positive range: NaN fails
    // every comparison, so the negated form `epsilon <= 0.0 || epsilon >= 1.0`
    // lets it straight through. Before the guards were rewritten this returned
    // 2^63 instead of throwing, which would have made estimate() loop for
    // effectively ever.
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        struct Bad { double eps, delta; const char* what; };
        std::vector<Bad> bad = {
            { nan,  0.05, "epsilon = NaN"      },
            { 0.02, nan,  "delta = NaN"        },
            { 0.0,  0.05, "epsilon = 0"        },
            { 1.0,  0.05, "epsilon = 1"        },
            {-0.1,  0.05, "epsilon negative"   },
            { 0.02, 0.0,  "delta = 0"          },
            { 0.02, 1.0,  "delta = 1"          },
        };
        int rejected = 0;
        for (auto& b : bad) {
            try {
                uint64_t n = hoeffdingSampleSize(b.eps, b.delta);
                std::cout << "  ** " << b.what << " was accepted, N=" << n << " **\n";
                ++failures;
            } catch (const std::invalid_argument&) {
                ++rejected;
            }
        }
        std::cout << "invalid argument rejection: " << rejected << "/" << bad.size()
                  << " rejected\n\n";
    }

    // --- A graph small enough to enumerate: 6 edges, 64 worlds ---
    // 1 and 2 point at each other, so there are cycles to get wrong
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
        bool ok = diff <= est.epsilon * 3;
        if (!ok) ++failures;
        std::cout << "[" << name << "] exact=" << exact
                  << "  monte_carlo=" << est.probability
                  << "  samples=" << est.samples
                  << "  |diff|=" << diff
                  << (ok ? "  OK" : "  ** OUT OF TOLERANCE **")
                  << "\n";
    }

    // --- Batch estimation across several pairs on the same worlds ---
    // Worth checking separately: batchEstimate() answers through a per-world
    // DBL index, a different code path from the lazy sampler estimate() uses.
    std::cout << "\nBatch estimation (Weighted Cascade, multiple query pairs):\n";
    UncertainGraph ugBatch = UncertainGraph::fromEdgeList(n, rawEdges, ProbabilityScheme::WeightedCascade, rng);
    MonteCarloSampler batchSampler(ugBatch, 555);
    std::vector<std::pair<int,int>> queries = {{0,3},{0,2},{1,3},{2,1},{3,0}};
    const double batchEpsilon = 0.02;
    auto results = batchSampler.batchEstimate(queries, batchEpsilon, 0.05);
    for (size_t i = 0; i < queries.size(); ++i) {
        double ex = exactReachability(ugBatch, queries[i].first, queries[i].second);
        double diff = std::abs(ex - results[i]);
        bool ok = diff <= batchEpsilon * 3;
        if (!ok) ++failures;
        std::cout << "  P(" << queries[i].first << " -> " << queries[i].second
                  << ") ~= " << results[i]
                  << "  exact=" << ex
                  << "  |diff|=" << diff
                  << (ok ? "  OK" : "  ** OUT OF TOLERANCE **") << "\n";
    }

    if (failures == 0) {
        std::cout << "\nALL SAMPLER TESTS PASSED\n";
    } else {
        std::cout << "\nSAMPLER TESTS FAILED  out_of_tolerance=" << failures << "\n";
        return 1;
    }
    return 0;
}
