#pragma once
#include "graph.hpp"
#include <vector>
#include <tuple>
#include <random>
#include <algorithm>
#include <stdexcept>

enum class ProbabilityScheme {
    Uniform,
    Trivalency,
    WeightedCascade
};

// Directed graph where every edge carries an existence probability in (0, 1].
// Stores only the edge list; callers that want adjacency build it from edges().
class UncertainGraph {
public:
    explicit UncertainGraph(int n = 0) : n_(n) {}

    int numVertices() const { return n_; }

    void addEdge(int u, int v, double p) {
        // Positive range, not its negation: NaN fails every comparison, so
        // `p <= 0.0 || p > 1.0` would wave it through
        if (!(p > 0.0 && p <= 1.0)) {
            throw std::invalid_argument("probability must be in (0, 1]");
        }
        n_ = std::max(n_, std::max(u, v) + 1);
        edges_.emplace_back(u, v, p);
    }

    const std::vector<std::tuple<int,int,double>>& edges() const { return edges_; }

    // Real datasets ship without probabilities, so they get assigned by one of
    // three schemes borrowed from the influence-maximisation literature.
    static UncertainGraph fromEdgeList(
        int n,
        const std::vector<std::pair<int,int>>& rawEdges,
        ProbabilityScheme scheme,
        std::mt19937& rng,
        double uniformP = 0.1) {

        UncertainGraph ug(n);

        if (scheme == ProbabilityScheme::WeightedCascade) {
            // p(u, v) = 1 / indeg(v), so every in-degree has to be final before
            // the first probability is assigned -- hence two passes
            Graph tmp(n);
            for (auto& e : rawEdges) tmp.addEdge(e.first, e.second);
            for (auto& e : rawEdges) {
                ug.addEdge(e.first, e.second, 1.0 / tmp.inDegree(e.second));
            }
            return ug;
        }

        std::uniform_int_distribution<int> triDist(0, 2);
        static const double triValues[3] = {0.1, 0.01, 0.001};

        // No `default:` -- a fourth scheme should draw a compiler warning here
        // rather than silently fall back to uniformP
        for (auto& e : rawEdges) {
            double p = uniformP;
            switch (scheme) {
                case ProbabilityScheme::Uniform:
                    p = uniformP;
                    break;
                case ProbabilityScheme::Trivalency:
                    p = triValues[triDist(rng)];
                    break;
                case ProbabilityScheme::WeightedCascade:
                    break; // unreachable; handled above
            }
            ug.addEdge(e.first, e.second, p);
        }
        return ug;
    }

    // Draw one possible world: each edge flips its own coin, independently
    Graph sampleWorld(std::mt19937& rng) const {
        Graph g(numVertices());
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        for (auto& [u, v, p] : edges_) {
            if (unit(rng) < p) g.addEdge(u, v);
        }
        return g;
    }

private:
    int n_ = 0;
    std::vector<std::tuple<int,int,double>> edges_;
};
