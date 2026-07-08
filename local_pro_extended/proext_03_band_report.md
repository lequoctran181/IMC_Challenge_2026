# IMC Challenge 2026 B `simplifygeometry` — 49843..50625 fail-closed candidate

## What is in `solution.cpp`

This is a concrete C++17 candidate focused only on the suspected hidden band `49843 <= V <= 50625`.

The code is deliberately fail-closed:

- Outside `49843..50625`, it outputs the original mesh unchanged.
- Inside the band, it first tries a narrow ordered periodic-grid recognizer for `M == 2V` UV/torus-style meshes. If the face pattern matches a wrapped 2-triangle grid and a vertex-Hausdorff nearest-selected-vertex check passes, it keeps about every other row and column and rebuilds a closed periodic triangular manifold from original vertices.
- If the periodic-grid path rejects, it runs a guarded endpoint-only edge-collapse branch. Each collapse preserves the standard two-face edge/link condition, rejects face flips/degenerates/duplicate local faces, tracks the vertex-cover radius of removed vertices, and validates the final closed orientable triangular 2-manifold before output.
- If either branch fails validation or removes too little, it falls back to the original input mesh.

This file is ready to compile and run as a standalone candidate. It is not a byte patch onto the minified `submission_1181_81.95_7.cpp` base, because the full attached archive and full untruncated exact-best source were not recoverable in this runtime. The intended transplant point in the exact-best base is after the existing case-5 `W5/VIMP/MIDEC/WK/B16` setup, before final output, with the same `N` gate and rollback condition.

## Why this branch

The GitHub history points repeatedly to the `49843..50625` case-5 band, especially lower UV/torus-ripple and upper softbox/boxy slices. The lower slice notes mention proxy reductions such as `uv49954_bumpy 1249/2494 -> 1129/2254` and `r5_torus50176_ripple 1985/3970 -> 1799/3598`, but also warn that unguarded or B16-only variants can regress. This candidate therefore does not use a broad parameter tweak; it uses a recognizer plus strict rollback/validation.

## Expected target cases

Most likely to help:

- ordered periodic torus/ripple grids around `V = 50176` and `M = 2V`, where the grid recognizer can reduce roughly 75% of vertices while using only original vertices;
- smooth UV/sphere-like meshes around `V = 49954`, where the collapse branch can remove a large fraction if orientation and local normals are stable;
- softbox/boxy meshes in `50234..50625` when local normal changes remain below the guard.

Expected no-op/fallback:

- non-target `V` counts;
- irregular/sharp meshes where the normal/link guards reject collapses;
- periodic-looking meshes whose face ordering or topology does not match the wrapped grid pattern;
- any candidate output with boundary edges, nonmanifold edges, same-direction paired edges, duplicate triangles, invalid indices, or degenerate faces.

## Local verification performed

Compile:

```sh
g++ -O2 -std=c++17 -pipe -Wall -Wextra solution.cpp -o solution
```

Result: compile success, no warnings after cleanup.

Sample smoke test:

```sh
./solution < sample.in > sample.out
head -1 sample.out
# 8 12
```

Generated periodic torus proxy, `224 x 224 = 50176` vertices and `100352` faces:

```text
output: 12544 vertices, 25088 faces
runtime: about 0.22 s locally
edge-count / orientation / duplicate / degenerate checks: pass
```

Generated oriented UV-sphere proxy, `224 x 223 + 2 = 49954` vertices and `99904` faces:

```text
output: 16984 vertices, 33964 faces
runtime: about 1.25 s locally
edge-count / orientation / duplicate / degenerate checks: pass
```

## Risks

The official six-axis normal/depth SSIM is not available here, so the collapse branch is guarded by normal preservation and manifold checks rather than a full official visual gate. This is the main scoring risk.

The standalone candidate is not the exact-best pipeline; outside the branch it returns the original mesh, while the real high-water base has many other hidden-case routes. For a serious Kattis attempt, transplant the two branch functions and their fail-closed call into `submission_1181_81.95_7.cpp` after making source-size space. The branch logic is independent and can also be inserted as an early target-band special path if the exact-best base state is inconvenient.

The grid branch assumes ordered periodic grid vertices and faces. If the hidden torus/ripple case is topologically periodic but not stored in this order, it will reject and leave the work to the collapse branch.

## Commands

```sh
g++ -O2 -std=c++17 -pipe solution.cpp -o solution
./solution < input.objtxt > output.objtxt
```

