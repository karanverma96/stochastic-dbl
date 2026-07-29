#include "dataset_loader.hpp"
#include "dbl_index.hpp"
#include "uncertain_graph.hpp"
#include "monte_carlo_sampler.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>

using Clock = std::chrono::high_resolution_clock;
double ms(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

int main(int argc, char** argv) {
    std::string path = argc > 1 ? argv[1] : "data/wiki-Vote.txt";
    std::cout << std::fixed << std::setprecision(4);

    // --- Load dataset ---
    auto t0 = Clock::now();
    LoadedEdgeList loaded = loadSnapEdgeList(path);
    auto t1 = Clock::now();
    std::cout << "Dataset: " << path << "\n";
    std::cout << "  nodes=" << loaded.numVertices
               << "  edges=" << loaded.edges.size()
               << "  load_time_ms=" << ms(t0, t1) << "\n\n";

    Graph g(loaded.numVertices);
    for (auto& e : loaded.edges) g.addEdge(e.first, e.second);

    // --- Build DBL index ---
    t0 = Clock::now();
    DBLIndex idx(g, /*k=*/64, /*kp=*/64);
    idx.build();
    t1 = Clock::now();
    double buildMs = ms(t0, t1);
    std::cout << "DBL index build_time_ms=" << buildMs << "\n\n";

    // --- Query throughput ---
    std::mt19937 rng(2026);
    std::uniform_int_distribution<int> vdist(0, loaded.numVertices - 1);
    int numQueries = 1000000;
    std::vector<std::pair<int,int>> queries(numQueries);
    for (auto& q : queries) q = {vdist(rng), vdist(rng)};

    t0 = Clock::now();
    long long reachableCount = 0;
    for (auto& q : queries) reachableCount += idx.query(q.first, q.second);
    t1 = Clock::now();
    double queryMs = ms(t0, t1);
    std::cout << "Query throughput:\n";
    std::cout << "  queries=" << numQueries
               << "  total_time_ms=" << queryMs
               << "  queries_per_sec=" << (numQueries / (queryMs / 1000.0))
               << "  avg_latency_us=" << (queryMs * 1000.0 / numQueries)
               << "  reachable_fraction=" << (double)reachableCount / numQueries << "\n\n";

    // --- Edge insertion cost ---
    int numInserts = 2000;
    std::vector<std::pair<int,int>> newEdges(numInserts);
    for (auto& e : newEdges) e = {vdist(rng), vdist(rng)};

    t0 = Clock::now();
    for (auto& e : newEdges) idx.insertEdge(e.first, e.second);
    t1 = Clock::now();
    double insertMs = ms(t0, t1);
    std::cout << "Edge insertion:\n";
    std::cout << "  inserts=" << numInserts
               << "  total_time_ms=" << insertMs
               << "  avg_latency_us=" << (insertMs * 1000.0 / numInserts) << "\n\n";

    // --- Uncertain graph / Monte Carlo pipeline across all 3 schemes ---
    std::cout << "Uncertain-graph Monte Carlo reachability (epsilon=0.02, delta=0.05):\n";
    std::vector<std::pair<int,int>> sampleQueries;
    for (int i = 0; i < 5; ++i) sampleQueries.push_back({vdist(rng), vdist(rng)});

    for (auto scheme : {ProbabilityScheme::Uniform, ProbabilityScheme::Trivalency, ProbabilityScheme::WeightedCascade}) {
        std::string name = scheme == ProbabilityScheme::Uniform ? "Uniform(p=0.1)"
                          : scheme == ProbabilityScheme::Trivalency ? "Trivalency"
                          : "WeightedCascade";
        UncertainGraph ug = UncertainGraph::fromEdgeList(loaded.numVertices, loaded.edges, scheme, rng, 0.1);
        MonteCarloSampler sampler(ug, 42);

        t0 = Clock::now();
        double totalP = 0.0;
        for (auto& q : sampleQueries) {
            auto est = sampler.estimate(q.first, q.second, 0.02, 0.05);
            totalP += est.probability;
        }
        t1 = Clock::now();
        double mcMs = ms(t0, t1);
        std::cout << "  [" << name << "] "
                  << "avg_estimate=" << (totalP / sampleQueries.size())
                  << "  total_time_ms=" << mcMs
                  << "  avg_time_per_query_ms=" << (mcMs / sampleQueries.size()) << "\n";
    }

    return 0;
}
