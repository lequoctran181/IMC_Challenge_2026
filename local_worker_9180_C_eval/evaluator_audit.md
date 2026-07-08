# Worker C evaluator audit

Scope: local-only audit for `simplifygeometry`. No Kattis submissions were made. New files and run artifacts are confined to:

`/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/local_worker_9180_C_eval/`

## What was inspected

- `local_worker_S8_proxy_eval/s8_proxy_harness.py` and `MANIFEST.md`: compile/run proxy harness, OBJ-like parser, closed-triangle checks, baseline deltas.
- `local_worker_BROAD_30/README_BROAD_30.md` and validator outputs: topology plus vertex-Hausdorff proxy reports, with render/SSIM section usually skipped.
- `local_worker_BROAD_31/outputs_highn_combo/*.validator.json`: richer structural fields, AABB, edge/area quantiles, sampled nearest-vertex ratio.
- `local_worker_S6_queue_verifier/queue_verifier.py`: source inventory, compile-status, source-size, route/pass fingerprinting, and Kattis-score filename mining.
- `local_worker_S12_score_mining/S12_score_mining_manifest.md`: local mapping from fetched/submission ids to rounded Kattis scores and route deltas.
- `worker_patches/briefs/central_brief_20260706_1129.md`: problem essentials and Kattis evidence.
- Official sample `/Users/TranAnh/Desktop/Competitive_programming/sg_work/1.in` and known sample outputs.

## Harness produced

`eval_harness.py` provides three subcommands:

```bash
./eval_harness.py run label=/path/to/candidate.cpp
./eval_harness.py compare label=/path/to/output.out [...]
./eval_harness.py analyze-output /path/to/output.out
```

Default sample input is `/Users/TranAnh/Desktop/Competitive_programming/sg_work/1.in`. Each run writes:

- `report.md`
- `summary.csv`
- `results.json`
- `build/*.bin`, `logs/*`, and `outputs/*.out` for compiled candidates

The `run` command rejects `--run-dir` outside this worker directory.

## Diagnostics implemented

Hard structural checks:

- header count consistency
- finite coordinates
- valid face indices
- zero/degenerate triangles
- duplicate raw faces and duplicate sorted-index faces
- exact duplicate vertices
- undirected edge incidence: boundary edges, nonmanifold `>2` edges, and self-loop edges
- same-direction paired edges, which catches orientation conflicts that can pass a pure two-face count
- face connected components
- Euler characteristic and genus estimate when the mesh is closed and orientable

Cheap score-adjacent features:

- output vertex/face counts
- compression percentage versus the input
- AABB diagonal/min/max
- edge-length and face-area quantiles
- deterministic sampled symmetric nearest-vertex ratio against `0.05 * input AABB diagonal`
- output SHA256 for byte-identical comparisons
- compile/runtime/source-size metadata

The ranking is deterministic: valid outputs first, then lower vertex count, then lower nearest-vertex ratio, then runtime.

## Smoke-test results

Output-only comparison:

```bash
./eval_harness.py compare \
  ans=/Users/TranAnh/Desktop/Competitive_programming/sg_work/1.ans \
  codex91=/Users/TranAnh/Desktop/Competitive_programming/sg_work/codex91_sample.out \
  hybrid=/Users/TranAnh/Desktop/Competitive_programming/sg_work/hybrid_sample.out \
  --run-dir local_worker_9180_C_eval/runs/compare_sample_outputs
```

All three outputs validate as `8 12`, hard-ok cube meshes with no bad edges, no orientation conflicts, no duplicate/degenerate faces, and nearest/tolerance ratio `0.1633`.

Compile/run path:

```bash
./eval_harness.py run \
  best1044=submission_1044_81.94_7.cpp \
  fetched19913561=fetched_sources/kattis_19913561_81.93_7.cpp \
  --run-dir local_worker_9180_C_eval/runs/compile_sample_best \
  --timeout 30
```

Both compiled successfully and returned valid sample outputs:

| candidate | sample V/F | hard_ok | runtime |
| --- | ---: | --- | ---: |
| `best1044` | `8 / 12` | yes | `0.408s` |
| `fetched19913561` | `8 / 12` | yes | `0.440s` |

Both produced the same sample output hash, matching the existing `codex91`/`hybrid` sample output.

## What correlates with Kattis

Strong correlation with avoiding rejection:

- Parseability and exact count consistency: Kattis expects the modified OBJ-like shape; malformed output is immediate failure.
- Valid indices, positive-area triangular faces, and closed two-face edge incidence: these mirror the hard validity constraints in the central brief.
- Same-direction edge conflicts: S8 notes showed some rows had `bad_edges=0` but orientation risk. Tracking this separately is useful because a pure manifold edge count can be too forgiving.
- Source size and compile status: current accepted sources are close to the `131072` byte limit, so source growth is a real submission risk.
- Sample first line `8 12`: not enough for score, but a fast sanity gate for cube simplification and output formatting.

Moderate correlation with accepted-score direction:

- Lower valid output vertex count is the objective, so it matters when the output also preserves the hidden constraints.
- Candidate-vs-baseline byte hash and count deltas are useful for local proxy sweeps: if a candidate is byte-identical on all relevant cases, it cannot improve those local cases.
- Nearest-vertex ratio is a useful approximation because the problem clarification says Hausdorff is vertex-based, not full surface-based.
- AABB, edge-length, area, component, genus, and valence-like summaries help identify hidden-family recognizers and detect destructive over-simplification.

Weak or missing correlation:

- Official score is hidden-case aggregate, not sample behavior. The official sample can be perfect while Kattis score ranges from very low to plateau.
- Local generated proxies are only useful if they resemble hidden cases. S8/BROAD reports show local improvements can plateau or regress on Kattis.
- Vertex-count improvement can be harmful if SSIM drops below the fixed six-view threshold.

## What is missing locally

- No official six-axis normal/depth render SSIM implementation in this harness.
- No hidden test distribution, weights, or per-case feedback.
- No exact Kattis score model; filenames and manifests often preserve only rounded scores such as `81.93`.
- No high-quality nearest-neighbor acceleration dependency. The harness uses deterministic brute-force sampling to stay portable; it is fine for sample/proxy smoke tests, not million-vertex exhaustive validation.
- No surface-intersection or self-intersection detection.
- No strict vertex-link disk/sphere validation. Edge incidence plus connected orientability catches many failures but is not a complete combinatorial 2-manifold proof.
- No automated proxy generation. S8 already owns proxy generation; this harness focuses on sample and reusable output diagnostics.
- No Kattis submission or Kattis manager integration by design.

## Recommendation

Use this harness as the first local gate for every candidate:

1. Compile and sample-run.
2. Require `hard_ok=true`, `bad_edges=0`, `same_direction_edges=0`, `degenerate_faces=0`, and `duplicate_faces=0`.
3. Compare candidate outputs against baseline/proxy outputs by SHA, V/F deltas, and nearest/tolerance ratio.
4. Escalate only structurally valid, non-byte-identical candidates to heavier proxy/render validation or central manager submission.

Do not treat this harness as proof of a 91.80+ improvement. Its value is preventing invalid or redundant candidates from entering the expensive Kattis loop.
