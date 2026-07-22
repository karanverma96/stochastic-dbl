#pragma once
#include "graph.hpp"
#include <bitset>
#include <vector>
#include <queue>
#include <algorithm>
#include <numeric>
#include <functional>

// DBL: DAG-free dynamic reachability index.
//
// Independent clean-room implementation of the indexing scheme described in
// Lyu et al., "DBL: Efficient Reachability Queries on Dynamic Graphs"
// (arXiv:2101.09441). Combines two complementary labels:
//
//   - DL (Dynamic Landmark) label: for a chosen landmark set L of size k,
//     DL_in(v)  = landmarks that can reach v
//     DL_out(v) = landmarks that v can reach
//     A shared landmark in DL_out(u) ∩ DL_in(v) proves u reaches v.
//
//   - BL (Bidirectional Leaf) label: leaf nodes are vertices with zero
//     in-degree ("sources") or zero out-degree ("sinks"), hashed into a
//     fixed-size bit vector of size k'.
//     BL_in(v)  = (hashed) source leaves that can reach v
//     BL_out(v) = (hashed) sink leaves that v can reach
//     If BL_in(u) ⊄ BL_in(v) or BL_out(v) ⊄ BL_out(u), then u cannot reach v.
//
// DL gives fast *positive* answers; BL gives fast *negative* answers.
// Anything neither label can resolve falls back to a pruned BFS.
//
// This first pass wires up construction (DL/BL labeling). Query and
// dynamic-update logic land in the next commit.

constexpr int kMaxLabelBits = 128; // upper bound on landmark / leaf-hash label size
using Bitset = std::bitset<kMaxLabelBits>;

inline bool isSubset(const Bitset& a, const Bitset& b) {
    // true iff every bit set in `a` is also set in `b`  (a ⊆ b)
    return (a & ~b).none();
}

class DBLIndex {
public:
    // k  = number of landmark nodes for the DL label
    // kp = number of hash buckets for the BL label
    DBLIndex(Graph& g, int k, int kp)
        : g_(g), k_(std::min(k, kMaxLabelBits)), kp_(std::min(kp, kMaxLabelBits)) {}

    // Batch-construct DL and BL labels from scratch (Section 4.1 / Algorithm 1).
    void build() {
        int n = g_.numVertices();
        DLin_.assign(n, Bitset{});
        DLout_.assign(n, Bitset{});
        BLin_.assign(n, Bitset{});
        BLout_.assign(n, Bitset{});

        selectLandmarks();
        selectLeaves();

        // DL construction: one forward BFS + one backward BFS per landmark.
        for (int i = 0; i < static_cast<int>(landmarks_.size()); ++i) {
            bfsMark(landmarks_[i], /*forward=*/true, i, DLin_);
            bfsMark(landmarks_[i], /*forward=*/false, i, DLout_);
        }

        // BL construction: source leaves populate BL_in via forward BFS,
        // sink leaves populate BL_out via backward BFS. Multiple leaves can
        // share a hash bucket, which is why BL only gives negative answers.
        for (int s : sourceLeaves_) {
            bfsMark(s, /*forward=*/true, hash(s), BLin_);
        }
        for (int t : sinkLeaves_) {
            bfsMark(t, /*forward=*/false, hash(t), BLout_);
        }
    }

    // TODO: implement in the next commit (query processing framework).
    bool query(int u, int v) const { return false; }

    // TODO: implement in the next commit (dynamic edge-insertion update).
    void insertEdge(int u, int v) { /* not yet implemented */ }

    // --- accessors, mainly for testing/inspection ---
    const Bitset& dlIn(int v) const { return DLin_[v]; }
    const Bitset& dlOut(int v) const { return DLout_[v]; }
    const std::vector<int>& landmarks() const { return landmarks_; }

private:
    Graph& g_;
    int k_, kp_;
    std::vector<int> landmarks_;
    std::vector<int> sourceLeaves_, sinkLeaves_;
    std::vector<Bitset> DLin_, DLout_, BLin_, BLout_;

    // Approximate centrality M(u) = |Pre(u)| * |Suc(u)|, take top-k (Section 4.1).
    void selectLandmarks() {
        int n = g_.numVertices();
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            long long ma = 1LL * g_.inDegree(a) * g_.outDegree(a);
            long long mb = 1LL * g_.inDegree(b) * g_.outDegree(b);
            return ma > mb;
        });
        int count = std::min(k_, n);
        landmarks_.assign(order.begin(), order.begin() + count);
    }

    // Leaf nodes = zero in-degree (sources) or zero out-degree (sinks).
    void selectLeaves() {
        sourceLeaves_.clear();
        sinkLeaves_.clear();
        for (int v = 0; v < g_.numVertices(); ++v) {
            if (g_.inDegree(v) == 0) sourceLeaves_.push_back(v);
            if (g_.outDegree(v) == 0) sinkLeaves_.push_back(v);
        }
    }

    int hash(int leafId) const {
        // Simple modulo hash into the k'-sized BL bucket space; collisions
        // are expected and are exactly why BL only proves negatives.
        return kp_ > 0 ? (leafId % kp_) : 0;
    }

    // BFS that sets `bit` in label[x] for every vertex x reachable from
    // `start` (forward = via Suc) or every vertex that can reach `start`
    // (forward = false, traversing via Pre).
    void bfsMark(int start, bool forward, int bit, std::vector<Bitset>& label) {
        std::queue<int> q;
        std::vector<char> visited(g_.numVertices(), 0);
        q.push(start);
        visited[start] = 1;
        while (!q.empty()) {
            int p = q.front(); q.pop();
            label[p].set(bit);
            const auto& next = forward ? g_.Suc(p) : g_.Pre(p);
            for (int x : next) {
                if (!visited[x]) {
                    visited[x] = 1;
                    q.push(x);
                }
            }
        }
    }
};
