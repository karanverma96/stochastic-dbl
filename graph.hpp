#pragma once
#include <vector>
#include <cstdint>
#include <algorithm>

// Simple directed graph stored as forward and reverse adjacency lists.
// Vertices are identified by contiguous indices [0, n).
class Graph {
public:
    explicit Graph(int n = 0) : n_(n), suc_(n), pre_(n) {}

    int numVertices() const { return n_; }

    // Grow the vertex set to at least `n` vertices (used when inserting
    // edges that reference new vertex ids).
    void ensureVertex(int v) {
        if (v >= n_) {
            n_ = v + 1;
            suc_.resize(n_);
            pre_.resize(n_);
        }
    }

    void addEdge(int u, int v) {
        ensureVertex(std::max(u, v));
        suc_[u].push_back(v);
        pre_[v].push_back(u);
    }

    const std::vector<int>& Suc(int u) const { return suc_[u]; }
    const std::vector<int>& Pre(int u) const { return pre_[u]; }

    int outDegree(int u) const { return static_cast<int>(suc_[u].size()); }
    int inDegree(int u) const { return static_cast<int>(pre_[u].size()); }

private:
    int n_;
    std::vector<std::vector<int>> suc_; // forward adjacency: Suc(u)
    std::vector<std::vector<int>> pre_; // reverse adjacency: Pre(u)
};
