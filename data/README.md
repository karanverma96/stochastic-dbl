# Dataset

`wiki-Vote.txt.gz` — Wikipedia administrator-election voting network from the
[SNAP dataset collection](https://snap.stanford.edu/data/wiki-Vote.html)
(Leskovec, Huttenlocher, Kleinberg — CHI 2010 / WWW 2010).

7,115 nodes, 103,689 directed edges. A directed edge A -> B means user A
voted on user B's request for adminship.

Unzip before running the benchmark:
```
gunzip -k wiki-Vote.txt.gz
../benchmark/benchmark wiki-Vote.txt
```
