# Certified Perceptual Mesh Simplification

> **2nd place · 93.830074 points · Accepted 7/7** in the 2026 IMC Challenge, Problem B (`simplifygeometry`).

This repository is the research artifact behind **NEU.AddictedTribes**' solution to *Perception-Aware Lossless Simplification of Million-Vertex 3D Meshes for Mobile Platforms*. The final program reduces 1,500,780 input vertices to 34,134 output vertices while satisfying the official topology, Hausdorff-distance, and six-view perceptual constraints.

[Read the paper (PDF)](paper/IMC_Challenge_Round2_NEU_AddictedTribes.pdf) · [Editable article (DOCX)](paper/IMC_Challenge_Round2_NEU_AddictedTribes.docx) · [Final Kattis source](submission/submission_20082703.cpp) · [Reproduce the score](docs/REPRODUCIBILITY.md)

![The end-to-end simplification and certification pipeline](paper/figures/pipeline.png)

## Result at a glance

| Case | Input vertices | Output vertices | Retained | Compression |
|---|---:|---:|---:|---:|
| Sphere-like sample | 4,098 | 25 | 0.6101% | 99.3899% |
| Armadillo | 23,201 | 4,340 | 18.7061% | 81.2939% |
| Bunny-like | 35,292 | 2,839 | 8.0443% | 91.9557% |
| Lucy | 49,987 | 3,030 | 6.0616% | 93.9384% |
| Slender | 377,084 | 7,400 | 1.9624% | 98.0376% |
| Nefertiti | 1,009,118 | 16,500 | 1.6351% | 98.3649% |
| **Total / official score** | **1,500,780** | **34,134** | **2.2744%** | **93.830074** |

The exact count-derived value is

$$
100\left(1-\frac{1}{6}\sum_{i=1}^{6}\frac{V'_i}{V_i}\right)
=93.83007422510956,
$$

which matches Kattis submission **20082703** after official rounding.

## What made the solution work

The solution is not a single decimator. It is a fail-closed hybrid pipeline:

1. **Guarded quadric-error contraction.** A compact 10-coefficient QEM core evaluates endpoint, midpoint, weighted, analytic, and segment candidates while enforcing link-condition, orientation, duplicate-face, and degeneracy guards.
2. **Cluster-normal memory.** Every surviving vertex carries the additive area-weighted normal sum of all original faces absorbed into its cluster. Collapse cost is therefore measured against original surface evidence, not only the progressively degraded mesh.
3. **Renderer-aware optimization.** The official six cameras, flat face normals, perspective depth, visibility, and SSIM windowing are reproduced locally. Offline topology flips, fan retriangulation, and coordinate fitting target the pixels that matter.
4. **Conservative geometric certification.** A cluster-radius recurrence provides a fast Hausdorff upper bound; independent topology, exact distance, normal-map, depth-map, and combined-SSIM checks reject unsafe candidates.
5. **Deterministic replay under hard limits.** Expensive search happens offline. The accepted 130,973-byte C++17 submission replays canonicalized, bit-packed edits with checkpoints and reference caches inside the 21-second and 131,072-byte limits.

The main research contribution is the bridge between local geometric simplification and the actual rendered objective: QEM proposes, cluster-normal memory preserves appearance history, the renderer identifies high-value structural edits, and validators decide what may be shipped.

## Artifact map

| Path | Purpose |
|---|---|
| [`paper/`](paper/) | Final Round 2 PDF/DOCX and publication figures |
| [`docs/ALGORITHM.md`](docs/ALGORITHM.md) | Mathematical model and complete method |
| [`docs/RESULTS.md`](docs/RESULTS.md) | Official result, milestones, and ablations |
| [`docs/HIDDEN_CONSTRAINT_WORKFLOW.md`](docs/HIDDEN_CONSTRAINT_WORKFLOW.md) | Controlled score decoding, invariant fingerprints, and frontier reconstruction |
| [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) | One-command checks and evaluator usage |
| [`release/final/`](release/final/) | Machine-readable release record and checksums |
| [`src/research/`](src/research/) | Readable experimental QEM/cluster-normal core |
| [`tools/`](tools/) | Score verifier and local perceptual evaluators |
| [`submission/`](submission/) | Byte-exact final source fetched back from Kattis and its result record |

## Quick verification

Requirements: Python 3 and a C++17 compiler. No third-party Python package is required for the release checks.

```bash
python3 tools/score_from_counts.py
python3 tools/verify_release.py
make verify
```

The verification checks the official score formula, output counts, source byte size, SHA-256 digest, article artifacts, and compilation of the exact accepted source.

## Read before reusing

No official test mesh, hidden input, or third-party model is redistributed. Public proxy meshes and contest inputs used during research remain excluded because their licensing may differ. See [NOTICE.md](NOTICE.md) for provenance and scope.

The code is released under the [MIT License](LICENSE). If this artifact helps your work, please use the metadata in [`CITATION.cff`](CITATION.cff) and cite the accompanying article.

## Further reading

- [Algorithm and mathematical model](docs/ALGORITHM.md)
- [Results and experiment lineage](docs/RESULTS.md)
- [Hidden-constraint discovery and evidence protocol](docs/HIDDEN_CONSTRAINT_WORKFLOW.md)
- [Reproducibility and integrity](docs/REPRODUCIBILITY.md)
- [Repository map and archival policy](docs/REPO_MAP.md)
- [Official problem statement](https://imc2.kattis.com/contests/imc2-2/problems/simplifygeometry)
