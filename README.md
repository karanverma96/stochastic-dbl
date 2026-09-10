# stochastic-dbl

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
**large, dynamic, uncertain directed graphs** — graphs that grow by edge
insertion over time, and where each edge carries a probability of existing.

In such graphs, the answer to an s-t reachability query is itself a
**probability**: the likelihood that a path exists from `s` to `t` across
all possible worlds implied by the edge probabilities. Computing this
exactly requires enumerating an exponential number of possible graphs,
making it intractable at scale — so approximate, sampling-based methods
are used instead.

---

## Approach

### 1. Base Index: DBL (Dynamic Landmark + Bidirectional Leaf)
Built on top of **DBL**, a published state-of-the-art dynamic reachability
index for large directed graphs (Lyu et al., *"DBL: Efficient Reachability
Queries on Dynamic Graphs"*, arXiv:2101.09441). DBL avoids maintaining a full
DAG/topological structure and instead uses **Dynamic Landmark (DL)** labels
combined with **Bidirectional Leaf (BL)** labels to answer reachability
queries efficiently, while supporting fast edge insertions. (DBL targets
*"insertion-only updates"*; the paper lists deletion as future work, and this
reimplementation keeps the same scope.)

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
| **Weighted Cascade** | For edge `(u, v)`, probability = `1 / indeg(v)` |

### 5. Dynamic Query Workload
Queries are issued at runtime against the current graph snapshot; after a
batch of queries, edges are inserted, and the *same* query batch is run
again against the updated index — testing the index's ability to keep
answering under graph evolution without full recomputation. `benchmark.cpp`
reports both rounds, so the shift in `reachable_fraction` shows what the
insertions actually changed.

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

## Build & Run

Everything except the three `.cpp` entry points is header-only, so there is no
build system — each program is a single `g++` invocation.

```bash
# Correctness tests (both exit non-zero on failure)
g++ -std=c++17 -O2 -o test_dbl     test_dbl.cpp     && ./test_dbl
g++ -std=c++17 -O2 -o test_sampler test_sampler.cpp && ./test_sampler

# Benchmark — unzip the dataset first
gunzip -k data/wiki-Vote.txt.gz
g++ -std=c++17 -O2 -o benchmark benchmark.cpp && ./benchmark data/wiki-Vote.txt
```

`test_dbl` runs in seconds. The benchmark's Monte Carlo section takes a few
minutes: at `epsilon=0.02, delta=0.05` Hoeffding's bound calls for 4,612
possible worlds per query estimate.

---

## Benchmarks

Machine: Intel Core i7-8750H (6 cores / 12 threads), Ubuntu 26.04.1, g++
15.2.0, flags `-std=c++17 -O2`. Dataset: `data/wiki-Vote.txt` (7,115 nodes,
103,689 edges, no duplicate rows).

This is a personal machine, not an isolated benchmarking rig, and its CPU
governor is `powersave` with a 800–4100 MHz range — a single run means very
little. Every figure below is the **median of 5 runs pinned to one core**
(`taskset -c 2`), after a discarded warm-up run, with the range given
alongside. The estimates and reachability fractions, by contrast, are
seed-fixed and came out byte-identical in all 5 runs.

Parameters: k=64 landmarks, k'=64 leaf-hash buckets, epsilon=0.02,
delta=0.05 — so 4,612 samples per estimate via Hoeffding's inequality.

| Metric | Median | Range (5 runs) |
|---|---|---|
| Dataset load time | 63.1 ms | 58.2 – 75.5 ms |
| DBL index build time | 47.1 ms | 40.2 – 49.1 ms |
| DBL query throughput | 54.6M queries/sec | 45.7M – 65.6M queries/sec |
| DBL query avg latency | 0.0183 µs | 0.0153 – 0.0219 µs |
| Edge insertion avg latency | 0.51 µs | 0.46 – 0.67 µs |
| Query throughput after 2,000 insertions | 57.2M queries/sec | 38.8M – 60.7M queries/sec |
| Monte Carlo, Uniform(p=0.1) | 53.8 ms/query | 52.4 – 61.0 ms/query |
| Monte Carlo, Trivalency | 14.4 ms/query | 12.8 – 17.2 ms/query |
| Monte Carlo, WeightedCascade | 0.33 ms/query | 0.31 – 0.46 ms/query |

Reachable fraction over 1M random pairs: **0.2353** before the insertions,
**0.3238** after — the same in every run.

**Why queries are this cheap.** A DBL query is a bitset intersection: it
checks whether a landmark in `DL_out(u)` also appears in `DL_in(v)` (a
handful of 128-bit AND/OR ops), with BL leaf-hash bits available to
short-circuit negative answers — no graph traversal at all in the common
case. Instrumenting the 1M-query run shows how rarely the fallback is
needed: 23.5% of queries are answered by DL, 75.9% by BL, and only 0.63%
reach the pruned BFS. (The two theorem-based early exits never fire on this
dataset — everything they could catch, BL has already caught.)

**Why Monte Carlo still costs milliseconds.** There is no closed form for
P(s reaches t), so each `estimate()` call draws 4,612 independent possible
worlds. It does not materialise them: `lazyReachable()` walks outward from
`s` and draws an edge's coin only when the search actually reaches that
edge. Edges the search never reaches cannot lie on an s-t path, so the
outcome is distributed exactly as if the whole world had been sampled first
— but on this dataset a full sample would waste 99% of its coin flips.

That is also why the three schemes differ so sharply in cost: lower edge
probabilities make the search die out sooner, so WeightedCascade's worlds
are explored barely at all compared with the flat p=0.1 Uniform scheme.

### What the estimates mean (and their resolution limit)

| Scheme | avg estimate over the 5 query pairs |
|---|---|
| Uniform (p=0.1) | 0.0573 |
| Trivalency | 0.0076 |
| WeightedCascade | 0.0000 |

The benchmark picks query pairs that are reachable when *every* edge is
present — without that filter roughly 76% of random pairs on this graph have
probability exactly zero, and the column would say nothing about the schemes.

`WeightedCascade` reporting exactly `0.0000` is the estimator's resolution
limit, not a defect. With N = 4,612 samples the smallest non-zero value
expressible is `1/4612 ≈ 2.2e-4`. The benchmark's pairs are 3 to 6 hops
apart, and at WeightedCascade's typical `p ≈ 0.023` a three-hop path carries
roughly `1.2e-5` — about 18x below that floor, so no sample ever succeeds.
Hoeffding's bound is on **absolute** error, so an estimator that always
returned 0 would satisfy it here too; resolving probabilities this small
needs a relative-error method instead.

Trivalency is the instructive contrast: its average `p` is comparable to
WeightedCascade's, yet it estimates `0.0076`. What matters is not the
average but whether paths exist whose edges are *all* high-probability —
Trivalency draws from {0.1, 0.01, 0.001}, so a three-hop path with all three
edges at 0.1 has probability `1e-3`, comfortably above the floor.

### Scale, relative to the paper

`wiki-Vote` is small. The DBL paper's smallest dataset is Email (265,214
vertices / 420,045 edges) and its largest is LiveJournal (4.8M / 69M); this
graph is roughly 40x smaller than the smallest of them. The numbers above
describe this implementation on this graph and should not be read as
reproducing or comparing against the paper's results.

---

## Key Learnings

- How a deterministic reachability index (DBL) can be extended to a
  probabilistic setting via sampling, without redesigning the core index
- Using Hoeffding's Inequality to rigorously bound sample size against
  accuracy/confidence requirements, instead of guessing `N`
- Trade-offs across the three edge-probability assignment schemes and how
  they affect reachability estimates on real-world graphs
- Handling reachability queries under graph evolution (edge insertions)
  without full index recomputation

---

## References

- Lyu, Q., Li, Y., He, B., Gong, B. *"DBL: Efficient Reachability Queries on
  Dynamic Graphs."* [arXiv:2101.09441](https://arxiv.org/abs/2101.09441)
- Hoeffding, W. (1963). *"Probability Inequalities for Sums of Bounded
  Random Variables."*
- SNAP Datasets — https://snap.stanford.edu/data/

---

## License

MIT License — see [LICENSE](./LICENSE).
