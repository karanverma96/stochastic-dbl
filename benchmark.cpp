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
               << "  duplicates_dropped=" << loaded.duplicatesDropped
               << "  malformed_lines=" << loaded.malformedLines
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

    // Pairs for the Monte Carlo section, picked before the insertions so they
    // match the edge list the UncertainGraph is built from. Only pairs that are
    // reachable with every edge present: about 76% of random pairs are not, and
    // those estimate zero under any scheme, which says nothing about the scheme.
    std::vector<std::pair<int,int>> sampleQueries;
    for (int tries = 0; tries < 1000000 && sampleQueries.size() < 5; ++tries) {
        int s = vdist(rng), t = vdist(rng);
        if (s != t && idx.query(s, t)) sampleQueries.emplace_back(s, t);
    }
    if (sampleQueries.empty()) {
        std::cerr << "no reachable vertex pair found; skipping the insertion, "
                     "second query and Monte Carlo sections\n";
        return 1;
    }

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

    // --- Query throughput again, against the updated graph ---
    // Same query set as the first round, so the two rows line up:
    // reachable_fraction shows what the insertions changed, and the timing
    // shows the index still answering without a rebuild.
    t0 = Clock::now();
    long long reachableAfter = 0;
    for (auto& q : queries) reachableAfter += idx.query(q.first, q.second);
    t1 = Clock::now();
    double queryMsAfter = ms(t0, t1);
    std::cout << "Query throughput after updates:\n";
    std::cout << "  queries=" << numQueries
               << "  total_time_ms=" << queryMsAfter
               << "  queries_per_sec=" << (numQueries / (queryMsAfter / 1000.0))
               << "  avg_latency_us=" << (queryMsAfter * 1000.0 / numQueries)
               << "  reachable_fraction=" << (double)reachableAfter / numQueries
               << "  (before=" << (double)reachableCount / numQueries << ")\n\n";

    // --- Uncertain graph / Monte Carlo pipeline across all 3 schemes ---
    std::cout << "Uncertain-graph Monte Carlo reachability (epsilon=0.02, delta=0.05):\n";
    std::cout << "  query pairs (all reachable in the deterministic graph):";
    for (auto& q : sampleQueries) std::cout << " " << q.first << "->" << q.second;
    std::cout << "\n";

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
