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

    // Test 4: the two theorem-based early exits, which fire zero times on
    // wiki-Vote because BL always rules the pair out first. Giving the graph no
    // sources and no sinks leaves every BL label empty, so blContain() can never
    // rule anything out and a query falls through to the theorems instead.
    // Instrumented, this graph reaches Theorem 1 four times and Theorem 2
    // sixteen times, with the pruned BFS never running.
    //
    // What this can and cannot catch, measured by mutating the two lines:
    //   Theorem 1 moved above the DL check  -> 6 mismatches, caught
    //   Theorem 2 changed to `return true`  -> 16 mismatches, caught
    //   either theorem deleted outright     -> 0 mismatches, NOT caught
    // The last one is not a hole in the test: both theorems are shortcuts, so
    // deleting one just hands the pair to the pruned BFS, which answers the
    // same. No output-based test can pin them; only their correctness is
    // testable here.
    {
        // Two 2-cycles joined one way, plus a detached 2-cycle. Every vertex
        // has indeg >= 1 and outdeg >= 1, so selectLeaves() finds nothing.
        int n = 6;
        Graph g(n);
        int edges[][2] = {{0,1},{1,0},{2,3},{3,2},{1,2},{4,5},{5,4}};
        for (auto& e : edges) g.addEdge(e[0], e[1]);

        DBLIndex idx(g, /*k=*/64, /*kp=*/64);
        idx.build();
        int checked = 0, mismatches = 0;
        checkAllPairs(g, idx, n, checked, mismatches);
        std::cout << "theorem-path test  pairs_checked=" << checked
                  << "  mismatches=" << mismatches << "\n";
        totalChecked += checked;
        totalMismatches += mismatches;
    }

    // Test 5: self-loop insertion, which every other insertion loop skips and
    // makeRandomGraph throws away. A self-loop can never change reachability --
    // query(u, u) is already true before it is added -- so the whole point is
    // that the labels come out unchanged and nothing downstream breaks.
    {
        std::mt19937 rng(4042);
        int n = 40;
        Graph g = makeRandomGraph(n, 80, rng);
        DBLIndex idx(g, /*k=*/4, /*kp=*/8);
        idx.build();

        int checked = 0, mismatches = 0;
        std::uniform_int_distribution<int> vertexDist(0, n - 1);
        for (int step = 0; step < 20; ++step) {
            idx.insertEdge(step % n, step % n);          // the self-loop
            checkAllPairs(g, idx, n, checked, mismatches);
            int u = vertexDist(rng), v = vertexDist(rng); // keep the graph moving,
            if (u != v) idx.insertEdge(u, v);             // so a stale label shows up
            checkAllPairs(g, idx, n, checked, mismatches);
        }
        std::cout << "self-loop insertion test  pairs_checked=" << checked
                  << "  mismatches=" << mismatches << "\n";
        totalChecked += checked;
        totalMismatches += mismatches;
    }

    // Test 6: the k / k' clamp -- the one place a caller's argument reaches
    // arithmetic directly. Below 1, kp = 0 makes hash()'s `leafId % kp_` divide
    // by zero; above 128, a bucket index runs past the bitset. landmarks()
    // exposes the clamped k, so the bound is asserted rather than inferred from
    // "it did not crash".
    //
    // Two graph sizes on purpose. selectLandmarks() takes min(k_, n), so on a
    // 30-vertex graph an unclamped k = 200 and a clamped k = 128 both come out
    // as 30 landmarks and the assertion cannot separate them -- the ceiling is
    // only observable above it. The n = 200 cases carry no all-pairs sweep
    // because they need none: an unclamped kp makes build() itself throw, long
    // before any query runs.
    //
    // Sensitivity, measured by mutating the constructor:
    //   clamp dropped entirely             -> SIGFPE, kp = 0 divides by zero
    //   ceiling dropped, lower bound kept  -> std::out_of_range from bitset::set(128)
    {
        std::mt19937 rng(5042);
        struct ClampCase { int n, k, kp, expectLandmarks; bool allPairs; const char* what; };
        std::vector<ClampCase> clampCases = {
            { 30,   0,   0,   1, true,  "k=0, kp=0 -> clamped up to 1" },
            { 30,  -5,  -5,   1, true,  "negative -> clamped up to 1"  },
            { 30, 200, 200,  30, true,  "ceiling still capped by n"    },
            {200, 200, 200, 128, false, "above the ceiling -> 128"     },
            {200, 128, 128, 128, false, "at the ceiling"               },
        };

        int checked = 0, mismatches = 0;
        for (auto& c : clampCases) {
            Graph g = makeRandomGraph(c.n, c.n * 2, rng);
            DBLIndex idx(g, c.k, c.kp);
            idx.build();
            int got = static_cast<int>(idx.landmarks().size());
            if (got != c.expectLandmarks) {
                ++mismatches;
                std::cerr << "CLAMP MISMATCH (" << c.what << "): landmarks="
                          << got << " expected=" << c.expectLandmarks << "\n";
            }
            if (c.allPairs) checkAllPairs(g, idx, c.n, checked, mismatches);
        }
        std::cout << "k/k-prime clamp test  pairs_checked=" << checked
                  << "  mismatches=" << mismatches << "\n";
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
