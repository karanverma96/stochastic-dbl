#pragma once
#include "graph.hpp"
#include <bitset>
#include <vector>
#include <queue>
#include <algorithm>
#include <numeric>
#include <functional>

// DBL: DAG-free dynamic reachability index. Clean-room implementation of
// Lyu et al., "DBL: Efficient Reachability Queries on Dynamic Graphs"
// (arXiv:2101.09441).
//
// Two labels, each of which can only answer one way:
//
//   DL (Dynamic Landmark) proves YES. Pick k landmarks; DL_in(v) is the
//   landmarks reaching v, DL_out(v) the ones v reaches. A landmark in both
//   DL_out(u) and DL_in(v) is a witnessed path u -> L -> v.
//
//   BL (Bidirectional Leaf) proves NO. Leaves are vertices with no in-edges
//   ("sources") or no out-edges ("sinks"), hashed into k' bits. If u reaches
//   v then every source reaching u also reaches v, and every sink v reaches u
//   reaches too -- so a broken containment rules the pair out.
//
// Neither label is complete, so whatever both stay silent about falls through
// to a pruned BFS.

constexpr int kMaxLabelBits = 128; // ceiling on both label sizes
using Bitset = std::bitset<kMaxLabelBits>;

inline bool isSubset(const Bitset& a, const Bitset& b) {
    return (a & ~b).none(); // nothing in a that is missing from b
}

class DBLIndex {
public:
    // k landmarks, kp leaf-hash buckets. Clamped to [1, kMaxLabelBits]: the
    // bitset caps the top, and below 1 a kp stops being a bucket count at all
    // -- 0 divides by zero in hash(), and a negative one asks selectLeaves()
    // for a bucket vector of negative size. (It would not make hash() return a
    // negative: `%` takes the sign of its left operand, which is a vertex id.)
    DBLIndex(Graph& g, int k, int kp)
        : g_(g),
          k_(std::clamp(k, 1, kMaxLabelBits)),
          kp_(std::clamp(kp, 1, kMaxLabelBits)) {}

    // Batch-construct DL and BL labels from scratch (Section 4.1 / Algorithm 1).
    void build() {
        int n = g_.numVertices();
        DLin_.assign(n, Bitset{});
        DLout_.assign(n, Bitset{});
        BLin_.assign(n, Bitset{});
        BLout_.assign(n, Bitset{});

        selectLandmarks();
        selectLeaves();

        // Forward from a landmark fills DL_in, backward fills DL_out
        std::vector<int> single(1);
        for (int i = 0; i < static_cast<int>(landmarks_.size()); ++i) {
            single[0] = landmarks_[i];
            bfsMark(single, /*forward=*/true, i, DLin_);
            bfsMark(single, /*forward=*/false, i, DLout_);
        }

        // Every leaf in a bucket sets the same bit, so one traversal per bucket
        // gives identical labels to one per leaf -- at most 2k' walks (one
        // forward per source bucket, one backward per sink bucket) instead of
        // one per leaf: 5,739 -> 128 on wiki-Vote.
        std::vector<std::vector<int>> srcBuckets(kp_), sinkBuckets(kp_);
        for (int s : sourceLeaves_) srcBuckets[hash(s)].push_back(s);
        for (int t : sinkLeaves_) sinkBuckets[hash(t)].push_back(t);

        for (int b = 0; b < kp_; ++b) {
            if (!srcBuckets[b].empty()) {
                bfsMark(srcBuckets[b], /*forward=*/true, b, BLin_);
            }
            if (!sinkBuckets[b].empty()) {
                bfsMark(sinkBuckets[b], /*forward=*/false, b, BLout_);
            }
        }
    }

    // Query processing framework (Algorithm 2). The order is load-bearing --
    // every line below the DL check argues from that check having already come
    // back false, so reordering them silently breaks correctness.
    bool query(int u, int v) const {
        if (u == v) return true;
        if (dlIntersect(u, v)) return true;   // a landmark witnesses u -> v
        if (!blContain(u, v)) return false;   // a leaf rules the pair out
        // Theorems 1 and 2: if v reaches u, or if either endpoint sits on a
        // cycle with a landmark, then a real u -> v path would have tripped
        // the DL check above. It didn't, so there is none.
        if (dlIntersect(v, u)) return false;
        if (dlIntersect(u, u) || dlIntersect(v, v)) return false;
        return prunedBFS(u, v);
    }

    // Edge-insertion update (Algorithm 3 and its BL counterpart). Propagation
    // runs before the edge is added, so the traversals walk the old adjacency;
    // that is exactly the set of vertices whose labels can change.
    void insertEdge(int u, int v) {
        g_.ensureVertex(std::max(u, v));
        int n = g_.numVertices();
        growLabels(n);

        // If DL already proves u reaches v there was a path before this edge,
        // so no pair becomes newly reachable and no label needs touching.
        // Otherwise both labels grow the same way: what could reach u now runs
        // forward from v, and what v could reach now runs backward from u.
        if (!dlIntersect(u, v)) {
            propagateForward(v, DLin_[u], DLin_);
            propagateBackward(u, DLout_[v], DLout_);
            propagateForward(v, BLin_[u], BLin_);
            propagateBackward(u, BLout_[v], BLout_);
        }

        g_.addEdge(u, v);
    }

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

    void growLabels(int n) {
        if (static_cast<int>(DLin_.size()) < n) {
            DLin_.resize(n);
            DLout_.resize(n);
            BLin_.resize(n);
            BLout_.resize(n);
        }
    }

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

    // Two ifs, not if/else: an isolated vertex is both a source and a sink
    void selectLeaves() {
        sourceLeaves_.clear();
        sinkLeaves_.clear();
        for (int v = 0; v < g_.numVertices(); ++v) {
            if (g_.inDegree(v) == 0) sourceLeaves_.push_back(v);
            if (g_.outDegree(v) == 0) sinkLeaves_.push_back(v);
        }
    }

    int hash(int leafId) const {
        // Collisions are expected, and they are why BL can only prove a
        // negative: a set bit means "some leaf in this bucket", not "this one"
        return leafId % kp_;
    }

    // Sets `bit` on everything reachable from `starts` (forward) or everything
    // that reaches them (backward). Multi-source so a whole BL bucket costs one
    // walk; landmarks pass a one-element vector.
    void bfsMark(const std::vector<int>& starts, bool forward, int bit,
                 std::vector<Bitset>& label) {
        std::queue<int> q;
        std::vector<char> visited(g_.numVertices(), 0);
        for (int s : starts) {
            if (!visited[s]) { visited[s] = 1; q.push(s); }
        }
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

    bool dlIntersect(int x, int y) const { return (DLout_[x] & DLin_[y]).any(); }

    bool blContain(int x, int y) const {
        return isSubset(BLin_[x], BLin_[y]) && isSubset(BLout_[y], BLout_[x]);
    }

    // Last resort, once neither label could answer. Both prunes below drop
    // vertices that provably cannot reach v, so the search stays small.
    bool prunedBFS(int u, int v) const {
        std::queue<int> q;
        std::vector<char> visited(g_.numVertices(), 0);
        q.push(u);
        visited[u] = 1;
        while (!q.empty()) {
            int w = q.front(); q.pop();
            for (int x : g_.Suc(w)) {
                if (x == v) return true;
                if (visited[x]) continue;
                // If a landmark already links u to x, then x reaching v would
                // have given u -> L -> x -> v and tripped the DL check in
                // query(). So x is a dead end for this query.
                if (dlIntersect(u, x)) continue;
                if (!blContain(x, v)) continue;       // leaves rule x -> v out
                visited[x] = 1;
                q.push(x);
            }
        }
        return false;
    }

    // OR `incoming` into everything forward of `start`, cutting a branch as
    // soon as a vertex already has those bits -- its descendants got them the
    // same way, back when they first arrived.
    //
    // By value, not by reference: callers pass label[u], an element of the very
    // vector this writes into. A reference happens to be safe today -- the only
    // write that could reach label[u] ORs it with itself, and nothing here
    // resizes the vector -- but the copy is 16 bytes and keeps that from being
    // a precondition of every future edit.
    void propagateForward(int start, Bitset incoming, std::vector<Bitset>& label) {
        if (incoming.none()) return;
        std::queue<int> q;
        q.push(start);
        while (!q.empty()) {
            int p = q.front(); q.pop();
            for (int x : g_.Suc(p)) {
                if (!isSubset(incoming, label[x])) {
                    label[x] |= incoming;
                    q.push(x);
                }
            }
        }
        // `start` itself is never subset-tested above -- the loop only tests
        // neighbours -- so this last write is bookkeeping, not a precondition
        // of the walk: updating it first visits exactly the same vertices.
        if (!isSubset(incoming, label[start])) {
            label[start] |= incoming;
        }
    }

    // Mirror of the above along Pre edges, for the DL_out / BL_out labels.
    // By value for the same reason.
    void propagateBackward(int start, Bitset incoming, std::vector<Bitset>& label) {
        if (incoming.none()) return;
        std::queue<int> q;
        q.push(start);
        while (!q.empty()) {
            int p = q.front(); q.pop();
            for (int x : g_.Pre(p)) {
                if (!isSubset(incoming, label[x])) {
                    label[x] |= incoming;
                    q.push(x);
                }
            }
        }
        if (!isSubset(incoming, label[start])) {
            label[start] |= incoming;
        }
    }
};
