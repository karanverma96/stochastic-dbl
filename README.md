# dbl-uncertain-reachability

An independent reimplementation exploring source-to-target (s-t) reachability
estimation on large-scale, **dynamic, uncertain directed graphs** — extending
a state-of-the-art deterministic reachability index to the probabilistic
setting.

> This project is inspired by academic research I conducted at IIT Jammu on
> uncertain graph reachability. The original codebase, datasets, and specific
> experimental results are confidential to the institution. This repository
> is a standalone, independently written reimplementation for demonstration
> purposes only — it does not contain any original source code, data, or
> design documents from that work.

---

## Problem Statement

Designing an efficient algorithm to answer s-t reachability queries on
**large, dynamic, uncertain directed graphs** — graphs where edges can be
inserted/deleted over time, and each edge carries a probability of existing.

In such graphs, the answer to an s-t reachability query is itself a
**probability**: the likelihood that a path exists from `s` to `t` across
all possible worlds implied by the edge probabilities. Computing this
exactly requires enumerating an exponential number of possible graphs,
making it intractable at scale — so approximate, sampling-based methods
are used instead.

---

## Approach

### 1. Base Index: DBL (Dynamic reachability via Bidirectional Leaf-DL labeling)
Built on top of **DBL**, a published state-of-the-art dynamic reachability
index for large directed graphs (Lyu et al., *"DBL: Efficient Reachability
Queries on Dynamic Graphs"*, arXiv:2101.09441). DBL avoids maintaining a full
DAG/topological structure and instead uses **Dynamic Landmark (DL)** labels
combined with **Bidirectional Leaf (BL)** labels to answer reachability
queries efficiently, while supporting fast edge insertions/deletions.

### 2. Extension to Uncertain Graphs
The deterministic DBL index answers "is `t` reachable from `s`?" with a
yes/no. To extend this to **uncertain graphs**, each edge `(u, v)` is
assigned an existence probability `p(u, v)`, and reachability is estimated
via **Monte Carlo sampling**:

- Sample the graph `N` times — each edge is included independently with
  probability `p(u, v)`.
- Run the DBL reachability check on each sampled ("possible world") graph.
- Estimate `P(s ⤳ t) ≈ (# samples where t is reachable) / N`.

### 3. Sample Size via Hoeffding's Inequality
Rather than choosing `N` arbitrarily, the required number of samples is
derived from **Hoeffding's Inequality**, based on a target error tolerance
`ε` and confidence level `1 - δ`:

```
N ≥ (1 / 2ε²) · ln(2 / δ)
```

This ensures the sampling-based probability estimate stays within `ε` of
the true value with probability at least `1 - δ`, avoiding both
under-sampling (poor accuracy) and over-sampling (wasted computation).

### 4. Edge Probability Assignment Schemes
Real-world datasets (e.g. from [SNAP](https://snap.stanford.edu/data/)) don't
come with edge probabilities, so probabilities are synthetically assigned
using one of three standard schemes from influence/uncertain-graph literature:

| Scheme | Rule |
|---|---|
| **Uniform** | All edges get the same fixed probability (e.g. 0.1 or 0.2) |
| **Trivalency** | Each edge's probability is drawn uniformly at random from `{0.1, 0.01, 0.001}` |
| **Weighted Cascade** | For edge `(u, v)`, probability = `1 / (deg(u) + deg(v))` |

### 5. Dynamic Query Workload
Queries are issued at runtime against the current graph snapshot; after a
batch of queries, edges are inserted/deleted, and further queries are run
— testing the index's ability to handle reachability estimation under
graph evolution without full recomputation.

---

## Project Structure

```
├── LICENSE
├── README.md
├── graph.hpp                  # Deterministic directed graph representation
├── uncertain_graph.hpp        # Uncertain graph + probability assignment schemes
├── dataset_loader.hpp         # SNAP dataset loader
├── dbl_index.hpp              # DBL reachability index (DL + BL labels)
├── hoeffding.hpp              # Hoeffding's-inequality-based sample size N
├── monte_carlo_sampler.hpp    # Monte Carlo sampling engine
├── test_sampler.cpp           # Sampler unit tests
├── test_dbl.cpp               # DBL correctness tests
├── benchmark.cpp              # Query throughput & accuracy measurement
└── data/
    ├── README.md               # Dataset attribution
    └── wiki-Vote.txt.gz        # SNAP Wiki-Vote dataset
```

---

## Benchmarks

Benchmarked in a Cowork-managed Linux VM (reported: AMD Ryzen 5 5625U with
Radeon Graphics, 2 cores, 2.8Gi RAM, Ubuntu 6.8.0-124-generic kernel),
compiler: g++ 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04.3), flags:
`-std=c++17 -O2`. Note: this reflects the VM environment's reported specs,
not necessarily raw host hardware — `lscpu` shows `Hypervisor vendor:
Microsoft`, confirming this is a virtualized 2-core allocation, not a
dedicated benchmarking machine.

Dataset: `data/wiki-Vote.txt` (7,115 nodes, 103,689 edges).

Benchmarked in a personal/local environment, not an isolated lab setup —
numbers are indicative of relative performance (DBL query speed vs. Monte
Carlo sampling cost), not absolute guarantees. The DBL index build/query/
insertion numbers below are the median of 3 full runs of
`./benchmark data/wiki-Vote.txt`. The Monte Carlo numbers come from the
same unmodified `MonteCarloSampler::estimate()` / `UncertainGraph` code,
run in per-query batches so each measurement fit inside the shell's
per-command time budget (a full run's Monte Carlo section takes several
minutes end-to-end); timings below are medians across 3 repeats of the
5-query, 3-scheme sweep the benchmark performs (k=64 landmarks, k'=64 leaf
hash buckets, epsilon=0.02, delta=0.05, so 4,612 samples per estimate via
Hoeffding's inequality).

| Metric | Median | Range (3 runs) |
|---|---|---|
| Dataset load time | 24.1 ms | 23.3 – 24.7 ms |
| DBL index build time | 444.5 ms | 443.6 – 447.1 ms |
| DBL query throughput | 118.6M queries/sec | 114.9M – 122.9M queries/sec |
| DBL query avg latency | 0.0084 µs | 0.0081 – 0.0087 µs |
| Edge insertion avg latency | 0.52 µs | 0.47 – 0.62 µs |
| Monte Carlo, Uniform(p=0.1) | 8,961 ms/query | 8,886 – 9,073 ms/query |
| Monte Carlo, Trivalency | 7,296 ms/query | 7,191 – 7,405 ms/query |
| Monte Carlo, WeightedCascade | 6,078 ms/query | 6,023 – 6,219 ms/query |

**Why the gap is this large.** A DBL query is a bitset intersection: it
checks whether a landmark in `DL_out(u)` also appears in `DL_in(v)` (a
handful of 128-bit AND/OR ops), with BL leaf-hash bits available to
short-circuit negative answers — no graph traversal at all in the common
case. That's why query throughput lands north of 100M queries/sec with
sub-hundredth-of-a-microsecond latency: the cost is a fixed, tiny number of
word-level bit operations regardless of graph size.

Monte Carlo reachability estimation is solving a fundamentally different
problem: instead of one deterministic graph, `UncertainGraph` assigns each
edge an existence probability, and there's no closed-form way to get
P(s reaches t) without sampling. Each `estimate()` call draws
`N = ceil(ln(2/delta) / (2*epsilon^2))` independent "possible worlds"
(4,612 of them at epsilon=0.02, delta=0.05), and for *each* sampled world
it re-samples every one of the 103,689 edges and runs a fresh BFS from `s`
looking for `t`. That's why one query estimate costs seconds instead of
nanoseconds: it's ~4,600 full graph traversals, not a handful of bit ops,
and the count is dictated by the Hoeffding accuracy bound, not by anything
tunable away without giving up statistical guarantees. The three
probability schemes differ in cost mainly through how dense the sampled
worlds end up (WeightedCascade's degree-based probabilities produce
sparser worlds than the flat p=0.1 Uniform scheme, which is reflected in
the ~35% spread between the fastest and slowest scheme above).

---

## Key Learnings

- How a deterministic reachability index (DBL) can be extended to a
  probabilistic setting via sampling, without redesigning the core index
- Using Hoeffding's Inequality to rigorously bound sample size against
  accuracy/confidence requirements, instead of guessing `N`
- Trade-offs across the three edge-probability assignment schemes and how
  they affect reachability estimates on real-world graphs
- Handling reachability queries under graph evolution (inserts/deletes)
  without full index recomputation

---

## References

- Lyu, B. et al. *"DBL: Efficient Reachability Queries on Dynamic Graphs."*
  [arXiv:2101.09441](https://arxiv.org/abs/2101.09441)
- Hoeffding, W. (1963). *"Probability Inequalities for Sums of Bounded
  Random Variables."*
- SNAP Datasets — https://snap.stanford.edu/data/

---

## License

MIT License — see [LICENSE](./LICENSE).
