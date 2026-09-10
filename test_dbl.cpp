#include "graph.hpp"
#include "dbl_index.hpp"
#include <iostream>
#include <random>
#include <queue>
#include <vector>

// Ground truth. Deliberately dumb -- a check is only worth as much as the
// thing it is checked against.
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

// `m` is attempts, not edges: self-loop draws are thrown away
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
    // Each block seeds its own generator. Sharing one meant that adding a test
    // reshuffled the data of every test after it.
    int totalChecked = 0, totalMismatches = 0;

    // Test 1: static correctness on a handful of random graphs of varying density.
    {
        std::mt19937 rng(42);
        struct Case { int n, m; };
        std::vector<Case> cases = {{20, 30}, {50, 150}, {100, 400}, {200, 1000}};

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
    }

    // Test 2: insertions into an index built on an empty graph, so every label
    // starts blank.
    {
        std::mt19937 rng(1042);
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

    // Test 2b: the real case -- insertions into a populated index, where they
    // have to merge into existing labels and vertices can stop being leaves.
    {
        // The cases are picked to give BL something to say. BL_in carries
        // source-leaf reachability, so a dense graph has few leaves, those
        // labels stay near-empty, and breaking the BL update changes nothing
        // you could observe. Sparse graphs leave plenty of leaves, and a small
        // k stops DL from answering before the BL check is reached.
        //
        // Checked by deleting each of the four propagation calls in
        // insertEdge: all four are caught here, from the first insertion.
        std::mt19937 rng(2042);
        struct DynCase { int n, m, inserts, k, kp; };
        std::vector<DynCase> dynCases = {
            {100,  50, 40,  2, 16},  // sparse, minimal DL: BL does the work
            {150, 150, 40,  2, 16},  // sparse
            {150, 150, 40, 16, 64},  // sparse, ordinary label sizes
            {150, 900, 40, 16, 32},  // dense: general coverage
            {120, 400, 40,  2,  4},  // tiny labels: maximum hash collisions
        };
        int checked = 0, mismatches = 0;
        for (auto& c : dynCases) {
            Graph g = makeRandomGraph(c.n, c.m, rng);
            DBLIndex idx(g, c.k, c.kp);
            idx.build();
            std::uniform_int_distribution<int> vertexDist(0, c.n - 1);
            for (int step = 0; step < c.inserts; ++step) {
                int u = vertexDist(rng), v = vertexDist(rng);
                if (u == v) continue;
                idx.insertEdge(u, v);
                checkAllPairs(g, idx, c.n, checked, mismatches);
            }
        }
        std::cout << "populated-index insertion test  pairs_checked=" << checked
                  << "  mismatches=" << mismatches << "\n";
        totalChecked += checked;
        totalMismatches += mismatches;
    }

    // Test 3: a scale the all-pairs tests cannot reach -- 2000 vertices is 4M
    // pairs, so a random sample stands in for exhaustive checking. Timing is
    // benchmark.cpp's job; this file only asserts.
    {
        std::mt19937 rng(3042);
        int n = 2000, m = 12000;
        Graph g = makeRandomGraph(n, m, rng);
        DBLIndex idx(g, 64, 64);
        idx.build();

        std::uniform_int_distribution<int> vertexDist(0, n - 1);
        int numQueries = 5000;
        int checked = 0, mismatches = 0;
        long long reachable = 0;
        for (int i = 0; i < numQueries; ++i) {
            int u = vertexDist(rng), v = vertexDist(rng);
            bool expected = bruteForceReachable(g, u, v);
            bool actual = idx.query(u, v);
            ++checked;
            reachable += expected;
            if (expected != actual) {
                ++mismatches;
                std::cerr << "MISMATCH q(" << u << "," << v << "): expected="
                          << expected << " actual=" << actual << "\n";
            }
        }
        std::cout << "large-graph sampled test  n=" << n << " m=" << m
                  << "  pairs_checked=" << checked
                  << "  mismatches=" << mismatches
                  << "  reachable_fraction=" << (double)reachable / numQueries << "\n";
        totalChecked += checked;
        totalMismatches += mismatches;
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
