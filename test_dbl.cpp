#include "graph.hpp"
#include "dbl_index.hpp"
#include <iostream>
#include <random>
#include <cassert>
#include <chrono>

// Brute-force ground truth: plain BFS reachability check.
bool bruteForceReachable(const Graph& g, int u, int v) {
    if (u == v) return true;
    std::vector<char> visited(g.numVertices(), 0);
    std::queue<int> q;
    q.push(u);
    visited[u] = 1;
    while (!q.empty()) {
        int p = q.front(); q.pop();
        for (int x : g.Suc(p)) {
            if (x == v) return true;
            if (!visited[x]) { visited[x] = 1; q.push(x); }
        }
    }
    return false;
}

Graph makeRandomGraph(int n, int m, std::mt19937& rng) {
    Graph g(n);
    std::uniform_int_distribution<int> vertexDist(0, n - 1);
    for (int i = 0; i < m; ++i) {
        int u = vertexDist(rng), v = vertexDist(rng);
        if (u != v) g.addEdge(u, v);
    }
    return g;
}

void checkAllPairs(const Graph& g, const DBLIndex& idx, int n, int& checked, int& mismatches) {
    for (int u = 0; u < n; ++u) {
        for (int v = 0; v < n; ++v) {
            bool expected = bruteForceReachable(g, u, v);
            bool actual = idx.query(u, v);
            ++checked;
            if (expected != actual) {
                ++mismatches;
                std::cerr << "MISMATCH q(" << u << "," << v << "): expected="
                          << expected << " actual=" << actual << "\n";
            }
        }
    }
}

int main() {
    std::mt19937 rng(42);

    // Test 1: static correctness on a handful of random graphs of varying density.
    struct Case { int n, m; };
    std::vector<Case> cases = {{20, 30}, {50, 150}, {100, 400}, {200, 1000}};

    int totalChecked = 0, totalMismatches = 0;
    for (auto& c : cases) {
        Graph g = makeRandomGraph(c.n, c.m, rng);
        DBLIndex idx(g, /*k=*/16, /*kp=*/32);
        idx.build();
        int checked = 0, mismatches = 0;
        checkAllPairs(g, idx, c.n, checked, mismatches);
        std::cout << "n=" << c.n << " m=" << c.m
                  << "  pairs_checked=" << checked
                  << "  mismatches=" << mismatches << "\n";
        totalChecked += checked;
        totalMismatches += mismatches;
    }

    // Test 2: dynamic correctness — insert random edges one at a time and
    // re-verify all-pairs reachability against brute force after each batch.
    {
        int n = 60;
        Graph g(n);
        DBLIndex idx(g, /*k=*/16, /*kp=*/32);
        idx.build();

        std::uniform_int_distribution<int> vertexDist(0, n - 1);
        int checked = 0, mismatches = 0;
        for (int step = 0; step < 300; ++step) {
            int u = vertexDist(rng), v = vertexDist(rng);
            if (u == v) continue;
            idx.insertEdge(u, v);
            if (step % 20 == 0) {
                checkAllPairs(g, idx, n, checked, mismatches);
            }
        }
        std::cout << "dynamic insertion test  pairs_checked=" << checked
                  << "  mismatches=" << mismatches << "\n";
        totalChecked += checked;
        totalMismatches += mismatches;
    }

    // Test 3: quick throughput measurement on a larger graph.
    {
        int n = 2000, m = 12000;
        Graph g = makeRandomGraph(n, m, rng);
        DBLIndex idx(g, 64, 64);
        idx.build();

        std::uniform_int_distribution<int> vertexDist(0, n - 1);
        int numQueries = 1000000;
        std::vector<std::pair<int,int>> queries(numQueries);
        for (auto& q : queries) q = {vertexDist(rng), vertexDist(rng)};

        auto start = std::chrono::high_resolution_clock::now();
        long long trueCount = 0;
        for (auto& q : queries) trueCount += idx.query(q.first, q.second);
        auto end = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        std::cout << "throughput test  n=" << n << " m=" << m
                  << "  queries=" << numQueries
                  << "  time_ms=" << ms
                  << "  reachable_fraction=" << (double)trueCount / numQueries << "\n";
    }

    std::cout << "\nTOTAL pairs_checked=" << totalChecked
              << " mismatches=" << totalMismatches << "\n";
    if (totalMismatches == 0) {
        std::cout << "ALL CORRECTNESS TESTS PASSED\n";
    } else {
        std::cout << "CORRECTNESS TESTS FAILED\n";
        return 1;
    }
    return 0;
}
