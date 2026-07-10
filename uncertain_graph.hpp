#pragma once
#include "graph.hpp"
#include <vector>
#include <tuple>
#include <random>
#include <stdexcept>

enum class ProbabilityScheme {
    Uniform,
    Trivalency,
    WeightedCascade
};

// A directed graph where every edge additionally carries an existence
// probability p(u, v) in (0, 1]. Built on top of the deterministic Graph
// structure used by the DBL index.
class UncertainGraph {
public:
    explicit UncertainGraph(int n = 0) : base_(n) {}

    int numVertices() const { return base_.numVertices(); }
    const Graph& base() const { return base_; }

    void addEdge(int u, int v, double p) {
        if (p <= 0.0 || p > 1.0) throw std::invalid_argument("probability must be in (0, 1]");
        base_.addEdge(u, v);
        edges_.emplace_back(u, v, p);
    }

    const std::vector<std::tuple<int,int,double>>& edges() const { return edges_; }

    // Assign probabilities to a deterministic edge list according to one of
    // the three standard schemes, producing an UncertainGraph.
    static UncertainGraph fromEdgeList(
        int n,
        const std::vector<std::pair<int,int>>& rawEdges,
        ProbabilityScheme scheme,
        std::mt19937& rng,
        double uniformP = 0.1) {

        UncertainGraph ug(n);

        if (scheme == ProbabilityScheme::WeightedCascade) {
            // Need in-degree to compute deg(u) + deg(v); build a plain
            // graph first to count degrees, then assign probabilities.
            Graph tmp(n);
            for (auto& e : rawEdges) tmp.addEdge(e.first, e.second);
            for (auto& e : rawEdges) {
                int u = e.first, v = e.second;
                int deg = tmp.inDegree(u) + tmp.outDegree(u)
                        + tmp.inDegree(v) + tmp.outDegree(v);
                double p = deg > 0 ? 1.0 / deg : 1.0;
                ug.addEdge(u, v, p);
            }
            return ug;
        }

        std::uniform_int_distribution<int> triDist(0, 2);
        static const double triValues[3] = {0.1, 0.01, 0.001};

        for (auto& e : rawEdges) {
            double p;
            switch (scheme) {
                case ProbabilityScheme::Uniform:
                    p = uniformP;
                    break;
                case ProbabilityScheme::Trivalency:
                    p = triValues[triDist(rng)];
                    break;
                default:
                    p = uniformP;
            }
            ug.addEdge(e.first, e.second, p);
        }
        return ug;
    }

    // Draw one "possible world": a plain Graph containing each edge
    // independently with its assigned probability.
    Graph sampleWorld(std::mt19937& rng) const {
        Graph g(numVertices());
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        for (auto& [u, v, p] : edges_) {
            if (unit(rng) < p) g.addEdge(u, v);
        }
        return g;
    }

private:
    Graph base_; // deterministic skeleton (ignores probabilities), used for degree lookups etc.
    std::vector<std::tuple<int,int,double>> edges_;
};
