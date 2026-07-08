# Worker E visual-shell cover candidate

## Artifact

- `worker_E_visual_shell_cover.cpp`
- `worker_E_visual_shell_cover_renamed.cpp`: source-limit-safe version submitted to Kattis.
- Built binary during local verification: `worker_E_visual_shell_cover`
- Sample output: `sample.out`

## Exact idea

This is a standalone C++17 candidate based on the current strong local simplification family, with a new early recognizer named `VSC` inserted before the ordinary edge-collapse pipeline.

The new part explicitly uses the clarification that Hausdorff distance is vertex-only:

1. Build a tiny mesh-archetype shell candidate:
   - axis-aligned bounding-box shell: 8 vertices, 12 triangles
   - axis-aligned octahedron shell: 6 vertices, 8 triangles
2. Verify it as a vertex cover:
   - every original input vertex must be within `0.0492 * bbox_diagonal` of one of the candidate shell triangles
   - every candidate output vertex must have an original input vertex within the same radius
3. Run the existing manifold validator and, for non-tiny instances, the visual proxy gate.
4. If the visual-shell candidate is accepted, output immediately. Otherwise fall back to the inherited high-score simplifier.

This is intentionally different from another local parameter tweak: it tries to jump by recognizing cases where a whole visual shell can be replaced by a tiny archetype under vertex-only Hausdorff, especially box-like closed meshes where all sampled vertices lie on the AABB shell.

## Limitations

- The new shell recognizer is conservative and only targets AABB-box and axis-octahedron archetypes.
- It will not help arbitrary smooth meshes, noisy organic shapes, or shells whose support vertices are not close to the proposed output vertices.
- For large non-archetype cases, it should normally reject and fall back to the inherited 81.94-ish pipeline.
- The AABB shell is intentionally early-exit; if it accepts incorrectly on a hidden small case, the later fallback pipeline is skipped. The vertex-cover test is the main guard.

## Verification

Compile command used:

```sh
g++ -O2 -std=c++17 -pipe -I/Users/TranAnh/Desktop/Competitive_programming/local_include /Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/local_worker_9180_E_breakthrough/worker_E_visual_shell_cover.cpp -o /Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/local_worker_9180_E_breakthrough/worker_E_visual_shell_cover
```

Official sample command used:

```sh
/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/local_worker_9180_E_breakthrough/worker_E_visual_shell_cover < /Users/TranAnh/Desktop/Competitive_programming/sg_work/1.in > /Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/local_worker_9180_E_breakthrough/sample.out
```

Sample gate result:

```text
8 12
```

The original artifact was 133317 bytes, above the Kattis 131072-byte source limit. The submitted version keeps behavior but shortens long identifiers, producing a 129127-byte file.

Kattis submission:

- Root archive: `submission_1237_81.94_7.cpp`
- Kattis id: `19915979`
- Result: `Accepted (81.938904)`, `7/7`, runtime `20.37 s`

This did not beat the exact best `81.945906`, but it validated the source-size workaround and preserved near-plateau behavior with the new VSC fail-closed branch present.
