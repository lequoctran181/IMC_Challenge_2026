# IMC Challenge 2026 B `simplifygeometry` candidate

## What changed

This is a complete C++17 candidate focused on the requested branch-selection change: instead of selecting only by vertex count or geometric thresholds, it generates several structurally valid simplifications and chooses between them with an internal six-view normal/depth proxy.

Implemented in `solution.cpp`:

- Flexible parser/output for both prefixed OBJ-like (`v`/`f`) and plain numeric formats.
- Fail-closed hard guards before any candidate can be selected:
  - valid indices;
  - no degenerate triangles;
  - closed undirected edge incidence, exactly two faces per edge;
  - output vertices <= input vertices;
  - vertex-wise symmetric nearest-vertex Hausdorff check against `0.04965 * input AABB diagonal`.
- Candidate families:
  - topology-preserving QEM edge-collapse snapshots at multiple target ratios;
  - sphere and ellipsoid remeshers;
  - cylinder/cone-like surface-of-revolution remeshers;
  - torus remesher;
  - row-major UV/sphere-grid subsampler.
- New six-axis visual proxy:
  - renders original and candidate from `+X`, `-X`, `+Y`, `-Y`, `+Z`, `-Z` orthographic cameras using AABB normalization;
  - z-buffers depth and face normals;
  - computes SSIM over depth and three normal channels on the union foreground mask;
  - includes silhouette IoU and a silhouette penalty;
  - ranks safe candidates by proxy threshold first and compression second.

## Expected target cases

This candidate is aimed at hidden cases where the previous router has multiple valid simplification choices but the old proxy/thresholds select the visually wrong one. It should be most useful for smooth closed surfaces, row-major UV spheres, ellipsoids/spheres, torus-like meshes, cones/cylinders, simple surfaces of revolution, and general meshes where QEM snapshots provide several valid compression levels for the proxy to choose among.

The intent is a fail-closed structural jump: parametric remeshers are only accepted after exact closed-mesh and vertex-wise Hausdorff checks; otherwise the code falls back to QEM candidates or the original mesh.

## Risks

- I could not run the official Kattis hidden tests or official renderer here, so the `91.80+` target is not proven.
- This is not a byte-level merge into the exact 131 KB plateau source; it is a smaller standalone candidate with the requested six-view proxy and branch selection.
- SSIM correlation can still differ if Kattis uses a different camera normalization, depth encoding, normal convention, or masking policy.
- The Hausdorff check is intentionally strict and vertex-wise, so some visually good parametric candidates may be skipped.

## Local checks performed

Compile check:

```bash
g++ -std=c++17 -O2 -pipe -static -s solution.cpp -o solution
```

Current source size:

```text
37384 bytes
```

Synthetic smoke tests run in this sandbox:

```text
tetra:   4 / 4       -> 4 / 4       valid closed, Hausdorff ratio 0.000
sphere:  762 / 1520  -> 174 / 344   valid closed, Hausdorff ratio 0.999
sphere2: 3446 / 6888 -> 345 / 686   valid closed, Hausdorff ratio 0.893
torus:   1200 / 2400 -> 240 / 480   valid closed, Hausdorff ratio 0.920
```

`Hausdorff ratio` is measured against the internal tolerance `0.04965 * input AABB diagonal` on vertices only.

## Compile/sample commands

```bash
g++ -std=c++17 -O2 -pipe -static -s solution.cpp -o solution
./solution < sample.in > sample.out
```

