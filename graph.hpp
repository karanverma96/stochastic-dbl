#pragma once
#include <vector>
#include <algorithm>

// Directed graph over contiguous vertex ids [0, n), kept in both directions so
// that backward traversals are as cheap as forward ones -- DBL needs both.
class Graph {
public:
    explicit Graph(int n = 0) : n_(n), suc_(n), pre_(n) {}

    int numVertices() const { return n_; }

    // Only ever grows, so addEdge can call it with max(u, v) and cover both
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
    // Named after the paper's Suc(u) / Pre(u)
    int n_;
    std::vector<std::vector<int>> suc_;
    std::vector<std::vector<int>> pre_;
};
