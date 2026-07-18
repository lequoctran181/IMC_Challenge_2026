# Source guide

## Exact accepted program

The authoritative submission is:

[`submission/submission_20082703.cpp`](../submission/submission_20082703.cpp)

It is 130,973 bytes and has SHA-256 `9195d42a73a6b85c8ae30d731f532175bdcd7c2982d421143d631b4c64b1a92c`. It was fetched back after Kattis reported Accepted 7/7. The dense macro layer and embedded payloads are source-budget engineering, not the recommended entry point for studying the algorithm.

## Readable research core

[`research/compact_qem_lab.cpp`](research/compact_qem_lab.cpp) contains the readable guarded-QEM laboratory implementation, including:

- 10-coefficient quadrics and multiple contraction candidates;
- versioned priority queue and local adjacency maintenance;
- manifold/link, degeneracy, orientation, and radius guards;
- projected-area normal cost;
- additive `clusterNormalSum` propagation and cluster-normal target cost;
- environment-controlled parameter sweeps.

The full renderer-aware replay system evolved through the archived orchestration branches. See [`docs/ALGORITHM.md`](../docs/ALGORITHM.md) for the clean integrated description and [`docs/REPO_MAP.md`](../docs/REPO_MAP.md) for provenance.

Build with:

```bash
make research
```
