# High-N structural recognizers for `simplifygeometry` 81.95 → 91.80+

## Scope

This report focuses only on strict recognizers for the high-N hidden region:

- `60000 <= N <= 1100000`
- current anchor: `submission_1181_81.95_7.cpp`, score about `81.945906`
- objective: find a large hidden structural case, not retune the existing B16/S3B16 tail

The core recommendation is to add fail-closed structural candidates that compete against the unchanged `1181` fallback. A recognizer should do nothing unless it proves the input is a known synthetic structure and the candidate output beats the anchor under topology, geometry, and visual gates.

## Executive plan

### Submit path under limited budget

| Priority | Candidate | Why it comes first | False-positive risk | Score upside |
|---:|---|---|---:|---:|
| 1 | `M == 2N` periodic torus/grid recognizer with exact lattice/face-order proof | The `M == 2N` relation is a strong closed genus-1 signature; a high-N periodic grid can compress by two orders of magnitude | Low if exact topology + face-order + degree gates are enforced | Very high |
| 2 | Exact AABB/box-shell recognizer | Cheap, easy to make no-op, and a dense cuboid shell can be replaced by 8 vertices / 12 faces | Very low if all six plane and axis-normal gates pass | Very high if present |
| 3 | Lat-long sphere/ellipsoid recognizer | `M == 2N - 4` plus two-pole ring topology is a strong synthetic signature | Medium; smooth organic shells can superficially fit an ellipsoid | Very high if present |
| 4 | Combined strict recognizer pack | Best if byte budget allows after packing; one submission covers all high-N hypotheses | Medium due source size/time interactions | Highest aggregate EV |
| 5 | Soft/rounded box or loose primitive fits | Only after exact recognizers fail | High | Unknown |

If only one or two submissions are available, submit the periodic `M == 2N` branch first, then the exact AABB branch. The lat-long branch is high-value, but it needs more gates to avoid firing on smooth non-primitive shells.

## Design principle: no-op unless confidence is extreme

Use hard AND gates, not a weighted classifier. A structural branch is allowed to output only if every required gate passes. If any metric cannot be computed within the time/memory budget, reject the candidate and run `1181` unchanged.

Recommended gate order:

1. Cheap numeric selector: `N`, `M`, bbox extents, and exact relation such as `M == 2N` or `M == 2N - 4`.
2. Cheap sampled geometry: bbox-plane occupancy, sampled face-normal clustering, residual sample, degree sample.
3. Full topology only for promising cases: closed edge incidence, same-direction edge conflict, duplicate faces, Euler characteristic, degree histogram.
4. Structural proof: lattice/face-order proof, ring/pole proof, or six-plane proof.
5. Candidate generation: build a small output mesh, usually several LOD variants.
6. Candidate tournament against the full `1181` fallback.
7. Emit candidate only if it wins by a meaningful vertex margin and passes high-resolution visual proxy gates.

Do not early-exit before the `1181` fallback unless the branch is an exact AABB cuboid with zero structural ambiguity. Even then, the safer implementation is still to store the candidate and compare it at the end.

## Shared feature packet

Add one compact `MeshFacts`-style packet computed from original input before destructive simplification. Keep this independent from the current collapse state.

### Required scalar facts

| Feature | Exact use |
|---|---|
| `N`, `M` | First gate for all high-N recognizers. |
| `m_eq_2n = (M == 2*N)` | Required for periodic torus/grid. |
| `m_eq_2n4 = (M == 2*N - 4)` | Strong signal for closed genus-0 triangulated shells such as sphere/ellipsoid/box. |
| bbox min/max and extents `ex, ey, ez` | AABB proof, aspect ratios, torus axis candidates, ellipsoid radii. |
| bbox diagonal `D` | All residuals should be normalized by `D`. |
| sorted extent ratios | Reject flat/degenerate cases and separate torus/sphere/box candidates. |
| coordinate scale epsilon | Use `epsD = max(1e-9 * D, 1e-12 * max_abs_coord, 1e-12)` as the base; use larger thresholds only where specified. |

### Required topology facts

Compute exactly only after cheap gates pass. For high-N, use a sorted edge list rather than a hash map if memory is safer.

| Fact | Required meaning |
|---|---|
| undirected edge multiplicity | Every undirected edge must have multiplicity exactly `2` for all closed structural branches. |
| same-direction paired edges | Must be `0`; this is the hard orientation veto. |
| duplicate sorted faces | Must be `0`. |
| degenerate faces | Must be `0`, using area threshold relative to `D^2`. |
| Euler characteristic `chi = N - E + M` | `chi == 0` for torus/grid; `chi == 2` for sphere/ellipsoid/box shell. |
| boundary edges | Must be `0`. |
| nonmanifold edges | Must be `0`. |
| degree histogram | Used differently per recognizer; see below. |

Patch-level implementation note: generate `3*M` packed edges. Sort by undirected key. Within each group, count multiplicity and orientation signs. In the same pass, compute `E`, boundary/nonmanifold counts, and same-direction conflicts. This avoids maintaining a huge unordered map.

### Required visual facts

Existing `visual_proxy_score` should not be trusted at low resolution for high-N structural replacement. Add a high-confidence proxy mode only for structural candidates:

- views: six axis views plus at least four diagonal views;
- compare candidate against original geometry, not only against current simplified state;
- include silhouette coverage and depth/normal disagreement if available in the existing proxy framework;
- require an absolute floor and a relative floor versus the `1181` output.

Recommended acceptance thresholds:

| Gate | Recommended default |
|---|---:|
| candidate visual proxy absolute floor | `>= 0.965` for exact AABB / analytic ellipsoid / analytic torus |
| candidate visual proxy floor for decimated original lattice | `>= 0.955` |
| candidate vs `1181` proxy | `candidate >= baseline - 0.002` |
| candidate vertex count margin | `V_cand <= 0.70 * V_1181` and `V_1181 - V_cand >= max(512, 0.05 * V_1181)` |
| time guard | if insufficient time remains for both proxy and fallback emit, reject candidate |

The vertex margin is important: a structural branch should only be accepted when it makes a large improvement. If it removes only a few hundred vertices from a high-N case, keep the safer anchor.

## Candidate tournament against `1181`

### Preferred execution model

Do not mutate the current `1181` pipeline to test structural output. Treat each recognizer as a candidate emitter.

1. After input parsing, compute cheap facts.
2. For any promising structural branch, build a small candidate output mesh from `originalP` and original faces. Store it separately as `candP/candF`.
3. Run the exact `1181` pipeline unchanged.
4. Just before the final output call, compute:
   - `V_base` using the existing current-output estimate;
   - `proxy_base` using the current simplified geometry if time allows;
   - `proxy_cand` for each structural candidate.
5. Accept the smallest candidate that satisfies all structural gates, validator gates, visual gates, and vertex-margin gates.
6. If no candidate wins, call the original `JD()` path exactly.

Patch-level instruction: replace the final direct output call with a tiny `choose_and_emit()` wrapper. The wrapper should either print `candP/candF` or call the original final emitter. It should not reorder, remove, or weaken the existing `GN → W2G/W2C → W5 → VIMP → MIDEC → WK → B16/S3B16 → JD` path.

### Why this model is safer

- It avoids snapshotting huge collapse state for `N` up to `1.1M`.
- A bad recognizer cannot corrupt the fallback geometry.
- Candidate generation can be done from original structure and then discarded.
- Same-hash lower scores suggest timing fragility; this model preserves the anchor unless the candidate is both cheap and clearly better.

## Recognizer 1: exact AABB / box shell

### Intended hidden case

A dense triangulated rectangular cuboid or exact six-plane box shell with many redundant vertices on the six bbox faces. This can often be replaced by an 8-vertex, 12-triangle cuboid while preserving silhouette and flat normals.

### Cheap selector

Required:

- `60000 <= N <= 1100000`
- `D > 0`
- all three extents positive: `min(ex,ey,ez) >= 1e-6 * D`
- topology relation should be sphere-like: prefer `M == 2*N - 4`; if not exact, allow continuing only until exact topology proves `chi == 2`

Immediate reject:

- `M == 2*N` with degree-6 majority; that is likely torus/grid, not box.
- bbox aspect has one near-zero extent; this is a sheet, not a closed shell.

### Bbox occupancy gate

For each vertex, compute its normalized distance to the nearest bbox plane:

`d_box(v) = min(|x-minx|, |x-maxx|, |y-miny|, |y-maxy|, |z-minz|, |z-maxz|) / D`

Required for exact AABB:

| Metric | Gate |
|---|---:|
| median `d_box` | `<= 2e-8` |
| p99 `d_box` | `<= 2e-6` |
| max `d_box` | `<= 2e-5` |
| vertices assigned to some bbox plane | `>= 0.9995 * N` |
| all six planes occupied | each plane has at least `max(32, 0.002*N)` vertices |

Use slightly looser `p99 <= 1e-5` only for a non-cube panel candidate, not for the 8-vertex cuboid output.

### Face plane gate

Classify each vertex by all bbox planes it is close to. A vertex on an edge can belong to two planes; a corner vertex can belong to three.

For every original triangle, require at least one bbox plane that contains all three vertices within `plane_eps = 2e-5 * D`.

Recommended gates:

| Metric | Gate |
|---|---:|
| triangles with a common bbox plane | `>= 0.999 * M` |
| triangles with no common plane | `0` for exact cuboid candidate; otherwise reject the 8-vertex candidate |
| area-weighted off-plane residual | p99 `<= 5e-6 * D` |

This rejects rounded boxes and smooth ellipsoids that merely touch the bbox.

### Axis-normal gate

For each face normal `n`, compute `axis_score = max(|nx|, |ny|, |nz|)`.

Required:

| Metric | Gate |
|---|---:|
| area-weighted fraction with `axis_score >= 0.9999` | `>= 0.995` |
| area-weighted fraction with `axis_score >= 0.999` | `>= 0.999` |
| all six normal directions represented | yes |
| inward/outward consistency by plane | no plane has mixed dominant direction beyond tiny numerical noise |

For a cuboid, normals should be almost perfectly axis-aligned. If this fails, do not output a cube.

### Topology and degree gates

Required:

- closed manifold: every undirected edge appears exactly twice;
- same-direction edge pairs: `0`;
- duplicate sorted faces: `0`;
- degenerate faces: `0`;
- Euler characteristic: `chi == 2`;
- no degree-1 or degree-2 vertices except impossible numerical artifacts; any such case rejects.

Degree distribution is not a primary positive signal for AABB because triangulated panels can have many valid degree patterns. It is only a safety veto.

### Face-order gate

Do not require face-list order for the AABB branch. Instead, use geometric six-plane grouping. If the input happens to group faces by panel, that can be used as a cheap early positive signal, but it should never replace the plane, normal, and topology gates.

Optional face-order evidence:

- faces form six large contiguous groups by dominant plane;
- within each group, normals are consistent;
- group total areas match bbox side areas within small relative tolerance.

If this optional evidence is absent, continue only if the exact geometric gates pass.

### Candidate generation

Generate LOD candidates in this order:

1. Exact cuboid: 8 bbox corner vertices, 12 triangles, outward orientation.
2. Six-panel coarse grid: only if exact cuboid visual proxy fails but plane gates pass. Use a small grid per face, such as `4x4` or `8x8`, with shared edges and corners. This still has far fewer vertices than `1181` but preserves face tessellation density if the renderer/proxy is sensitive to triangulation.

For exact AABB, the 8-vertex candidate should be preferred. The panel grid is a fallback, not a first submission target.

### Candidate visual gates

Accept exact cuboid only if:

- output is closed and oriented;
- output bbox equals input bbox exactly within `epsD`;
- sampled input-to-cuboid surface distance p99 `<= 2e-5 * D`;
- high-resolution proxy `>= 0.965`;
- candidate proxy `>= baseline proxy - 0.002`;
- candidate vertex count is at least `30%` lower than baseline output.

### Risk notes

The main risk is firing on a rounded box or a primitive with flat-ish sides. The face-plane and axis-normal gates should reject those. Do not add a softbox branch in the same first submission.

## Recognizer 2: lat-long sphere / ellipsoid

### Intended hidden case

A high-resolution UV sphere or ellipsoid built from latitude rings and longitude columns, usually topological genus 0:

- `N = 2 + U * R`
- `M = 2 * U * R = 2*N - 4`
- two pole vertices or two tiny pole clusters
- most interior vertices degree 6

A correct branch can replace hundreds of thousands of vertices with a few thousand while keeping visual SSIM above the threshold.

### Cheap selector

Required:

- `60000 <= N <= 1100000`
- `M == 2*N - 4`
- bbox extents all positive and not extremely flat:
  - `min_extent / max_extent >= 0.05` for general ellipsoid;
  - use `>= 0.20` for first submission if you want the strictest version.
- not AABB-like:
  - bbox-plane occupancy must be low enough to avoid exact box, e.g. fewer than `20%` of vertices near bbox planes.

Immediate reject:

- `M == 2*N`
- degree-6 majority plus `chi == 0` after topology check
- high bbox face occupancy with axis normals, which belongs to the AABB branch

### Topology gate

Required:

- closed manifold: every undirected edge appears exactly twice;
- same-direction paired edges: `0`;
- duplicate sorted faces: `0`;
- degenerate faces: `0`;
- `chi == 2`;
- degree histogram consistent with lat-long topology.

Recommended degree gates:

| Metric | Gate |
|---|---:|
| degree-6 interior fraction | `>= 0.95` of vertices excluding pole candidates |
| high-degree pole count | exactly `2` vertices, or exactly `2` small pole clusters |
| pole valence equality | top and bottom pole valence differ by at most `1%` |
| vertices with degree outside expected set | `<= 0.5%` |

If there are not two clear pole vertices or pole clusters, reject. A generic smooth sphere-like mesh is not enough.

### Ring/pole structural proof

Try each coordinate axis as the polar axis. For each axis:

1. Cluster vertices by coordinate along the axis using tolerance `ring_eps = max(2e-7 * D, 16*epsD)`.
2. Require two extreme clusters to be singleton poles or tiny pole clusters.
3. Require all middle rings to have the same size `U`, within zero tolerance for the strict first version.
4. Require `N = pole_count_top + pole_count_bottom + U*R`.
5. Require `M = 2*U*R`.
6. For every middle ring, sorted angular order around the polar axis must form a closed cycle.
7. Each ring vertex should connect to two same-ring neighbors plus the expected adjacent-ring vertices.

Strict first-submission recommendation:

- accept only singleton poles;
- accept only equal-size middle rings;
- accept only axis-aligned polar axis from bbox axes;
- reject arbitrary PCA axes to save code and avoid false positives.

### Face-order gate

Use face-order as a positive proof if the generator uses row-major vertices/faces. It should validate one of a small set of patterns:

- vertices arranged as `[top pole]`, then `R` rings of size `U`, then `[bottom pole]`; or bottom/top reversed;
- top cap has `U` triangles sharing top pole;
- each adjacent ring strip has `2U` triangles;
- bottom cap has `U` triangles sharing bottom pole;
- longitude wraps exactly.

Recommended acceptance:

| Metric | Gate |
|---|---:|
| expected cap/strip triangles matched | `>= 0.999 * M` |
| wrong diagonal / missing wrap count | `0` for strict branch |
| face orientation consistency | no same-direction edge conflicts |

If face-order is not matched but the adjacency ring proof is exact, the branch may still proceed. For first submission, requiring face-order is safer and smaller.

### Ellipsoid residual gate

Fit the ellipsoid using bbox center and half-extents:

- center `c = (bbox_min + bbox_max) / 2`
- radii `rx, ry, rz = extents / 2`
- normalized radial value `q(v) = sqrt(((x-cx)/rx)^2 + ((y-cy)/ry)^2 + ((z-cz)/rz)^2)`
- residual `r(v) = |q(v) - 1|`

Required for analytic ellipsoid candidate:

| Metric | Gate |
|---|---:|
| median radial residual | `<= 2e-5` |
| p95 radial residual | `<= 2e-4` |
| p99 radial residual | `<= 8e-4` |
| max radial residual | `<= 2e-3` |
| pole residual | `<= 1e-5` |

Also check normal consistency:

- area-weighted median dot between face normal and ellipsoid normal `>= 0.995`;
- p5 dot `>= 0.985`;
- no large inverted-normal clusters.

### Candidate generation

Generate multiple LOD candidates and choose the smallest passing visual gate.

Recommended LODs:

| Candidate | Approx vertices | Use case |
|---|---:|---|
| `U'=32`, `R'=16` | about `514` | Very smooth sphere/ellipsoid; high risk if silhouette too coarse. |
| `U'=48`, `R'=24` | about `1154` | Good first target. |
| `U'=64`, `R'=32` | about `2050` | Safer visual fallback. |
| `U'=96`, `R'=48` | about `4610` | Still huge compression for `N >= 60000`; use if proxy demands it. |

Two safe candidate modes:

1. Analytic ellipsoid mesh from fitted center/radii. Use only if residual gates are extremely tight.
2. Ring-sampled mesh from original rings. Select a subset of original latitude rings and longitudes, preserving original vertex coordinates. This handles mild non-ellipsoid latitude profiles better and is safer for SSIM.

First submission recommendation: implement ring-sampled mesh if ring/face order is exact; otherwise implement analytic ellipsoid only if code budget is much smaller. Ring-sampled output has lower false-negative risk on slightly deformed ellipsoids.

### Candidate visual gates

Accept only if:

- topology of candidate is closed genus 0;
- output has no degenerate cap triangles;
- bbox matches original within `1e-4 * D` for analytic candidate or exactly for sampled candidate;
- high-resolution proxy `>= 0.955` for sampled candidate or `>= 0.965` for analytic candidate;
- candidate proxy `>= baseline proxy - 0.002`;
- candidate vertex count beats baseline by the required margin.

### Risk notes

The dangerous false positive is a smooth organic closed shell with `M == 2N - 4` and approximate ellipsoid bbox. The ring/pole proof is mandatory. Do not accept from residual alone.

## Recognizer 3: periodic torus/grid with `M == 2N`

### Intended hidden case

A closed periodic triangular grid, most likely torus-like or UV-periodic, with two triangles per logical cell:

- `N = U * V`
- `M = 2 * U * V = 2N`
- closed manifold
- Euler characteristic `0`
- degree-6 majority, often exactly all degree 6
- row-major or column-major face order in generator output

This is the highest-priority high-N recognizer because `M == 2N` is a strong topological separator from sphere/box shells.

### Cheap selector

Required:

- `60000 <= N <= 1100000`
- `M == 2*N`
- bbox extents all positive: `min_extent >= 1e-6 * D`
- not exact AABB: bbox-plane occupancy below exact AABB threshold

Immediate reject:

- `M == 2*N - 4`
- any obvious boundary evidence from sampled edges or degree sample
- degree sample has too many vertices outside degree 6

### Topology gate

Required:

- closed manifold: every undirected edge multiplicity exactly `2`;
- same-direction paired edges: `0`;
- duplicate sorted faces: `0`;
- degenerate faces: `0`;
- `chi == 0`;
- degree-6 fraction:
  - strict first submission: `>= 0.995`;
  - ideal exact grid: `1.000`.

Reject if there are pole-like high-degree vertices. That belongs to a lat-long sphere, not a torus/grid.

### Factorization gate

Find integer factors `U,V` such that `U*V == N`.

Recommended constraints:

| Constraint | Gate |
|---|---:|
| both dimensions | `U >= 16`, `V >= 16` |
| aspect ratio | `max(U,V)/min(U,V) <= 16` for first submission |
| face count relation | `2*U*V == M` exactly |
| coordinate periodicity | wrap edge lengths comparable to internal edge lengths |

Try several plausible factor pairs. Rank by face-order match first, coordinate smoothness second.

### Face-order / index-lattice proof

This is the critical safety gate. For each factor pair and orientation, test whether vertex IDs can be interpreted as a periodic grid.

Try these index maps:

- row-major: `id = u*V + v`
- column-major: `id = v*U + u`
- with optional reversal of `u`, `v`, or both
- optional diagonal pattern A/B for the two triangles of each quad

For every logical cell `(u,v)`, the expected triangles are one of two diagonal patterns:

- Pattern A: triangles use `(u,v)`, `(u+1,v)`, `(u,v+1)` and `(u+1,v)`, `(u+1,v+1)`, `(u,v+1)`.
- Pattern B: triangles use the opposite diagonal.

Indices wrap modulo `U,V`.

Recommended gates:

| Metric | Gate |
|---|---:|
| expected triangles matched as unordered triples | `>= 0.9995 * M` |
| missing cells | `0` for strict version |
| duplicated cells | `0` |
| diagonal pattern consistency | one dominant diagonal pattern `>= 0.995` or exact alternating pattern detected |
| orientation conflicts | `0` |

First submission recommendation: require exact or near-exact row/column-major face-order proof. Do not attempt expensive arbitrary lattice reconstruction in the first branch.

### Coordinate periodicity gate

After factor and index map are selected:

- compute edge lengths along `u`, `v`, and diagonal directions;
- compare wrap edges against interior edges;
- require no seam discontinuity.

Suggested gates:

| Metric | Gate |
|---|---:|
| p95 wrap-edge length / p95 interior-edge length | between `0.5` and `2.0` |
| median neighbor length CV per direction | `<= 0.35` |
| max single-step coordinate jump | `<= 8 * median_step` |
| diagonal length consistency | within `3x` median expected diagonal |

These thresholds should reject an arbitrary vertex ordering that accidentally factors.

### Torus residual gate

For a torus-specific analytic candidate, fit each bbox axis as the torus normal axis and keep the best.

For axis `z` as an example:

- center `c` from bbox center or mean;
- radial distance in major plane: `rho = sqrt((x-cx)^2 + (y-cy)^2)`;
- major radius `R = median(rho)`;
- tube radius sample `a_i = sqrt((rho_i - R)^2 + (z_i-cz)^2)`;
- tube radius `a = median(a_i)`;
- residual `|a_i-a|/D`.

Required for analytic torus output:

| Metric | Gate |
|---|---:|
| p50 tube residual | `<= 2e-4` |
| p95 tube residual | `<= 1e-3` |
| p99 tube residual | `<= 3e-3` |
| major radius | `R > 2*a` unless known spindle/self-intersecting torus is intended |
| axis aspect | two larger bbox extents roughly comparable |

Do not require torus residual if using a decimated original periodic grid candidate; the lattice proof plus visual proxy is safer for arbitrary periodic surfaces.

### Candidate generation

Preferred first-submission candidate: decimated original periodic grid, not analytic torus.

Reason: an `M == 2N` hidden case may be a torus, ripple torus, periodic height field, or another UV-periodic shell. A subgrid sampled from original coordinates preserves the generated shape better than a fitted torus.

Generate several LODs preserving aspect ratio:

| Candidate type | Rule | Typical output |
|---|---|---:|
| coarse | choose `U',V'` near `32x32` but proportional to `U,V` | about `1024` vertices |
| medium | choose near `48x48` | about `2304` vertices |
| safe | choose near `64x64` or `96x48` | `4096–4608` vertices |
| very safe | choose near `96x96` only if proxy requires | about `9216` vertices |

For each output grid vertex, select the nearest original grid index under the detected lattice map. Build wrapped faces using the same diagonal orientation as the input. This avoids interpolation code and keeps vertices on the original surface.

If analytic torus residual is extremely tight, optionally test analytic torus LODs too. The first submission does not need this; decimated original lattice is safer.

### Candidate visual gates

Accept only if:

- candidate topology is closed genus 1: `M_cand == 2*V_cand`, `chi == 0`;
- no duplicate/degenerate faces after downsampling;
- sampled candidate-to-original vertex distance is small along selected grid lines;
- high-resolution proxy `>= 0.955`;
- candidate proxy `>= baseline proxy - 0.002`;
- candidate vertex count beats baseline by required margin.

### Risk notes

The main false positive is an arbitrary closed genus-1 mesh with `M == 2N` but no row-major lattice. The face-order proof and coordinate periodicity gate are mandatory. Do not accept from `M == 2N` plus degree histogram alone.

## Required gates by recognizer

| Gate type | AABB shell | Lat-long ellipsoid | Periodic torus/grid |
|---|---|---|---|
| N band | `60000..1100000` | `60000..1100000` | `60000..1100000` |
| M relation | usually `M == 2N-4`, or exact `chi==2` | `M == 2N-4` | `M == 2N` |
| Euler | `chi == 2` | `chi == 2` | `chi == 0` |
| closed manifold | required | required | required |
| same-direction edge conflicts | `0` | `0` | `0` |
| duplicate faces | `0` | `0` | `0` |
| degree | safety veto | degree-6 interior + 2 poles | degree-6 majority/exact |
| bbox | all vertices near six planes | center/radii fit, not AABB | not AABB; torus aspect if analytic |
| face-order | optional evidence only | ring/cap/strip proof; require for strict first version | mandatory row/column-major lattice proof |
| residual | plane residual | ellipsoid radial + normal residual | lattice smoothness; torus residual optional |
| visual | `>=0.965` exact cube | `>=0.955/0.965` depending candidate | `>=0.955` decimated grid |
| fallback | always restore/emit `1181` unless candidate wins | same | same |

## Patch-level integration instructions

### Add compact candidate representation

Add a small candidate object with:

- kind ID: AABB, LATLONG, PERIODIC_GRID;
- confidence flag and reason bits;
- candidate vertices and faces;
- predicted topology class;
- candidate vertex count;
- structural residual summary;
- optional cached proxy score.

This object is output-only. It should not share mutable state with `P`, `faces`, `BU`, `BR`, or the current collapse adjacency.

### Add high-N recognizer dispatcher

At the end of input parsing or immediately after existing global feature initialization:

1. If `N < 60000` or `N > 1100000`, do nothing.
2. Run cheap selectors in this order:
   - exact `M == 2N` → periodic grid candidate;
   - exact/sphere-like `M == 2N - 4` plus AABB occupancy → AABB candidate;
   - exact `M == 2N - 4` plus low AABB occupancy → lat-long candidate.
3. Stop after building a small number of candidates, preferably no more than two.
4. If elapsed time is already high, discard all candidates and use fallback.

### Add `choose_and_emit()` wrapper

Replace the final `JD()` call with a wrapper that:

1. obtains `V_base` from the current simplified state;
2. evaluates candidates only if `V_cand` is much lower than `V_base`;
3. computes high-resolution proxy for candidate and baseline if time allows;
4. validates candidate output topology and face orientation;
5. prints candidate if it wins;
6. otherwise calls original `JD()` unchanged.

Do not change the `1181` simplification route in the same submission. The goal is to add a fallback competitor, not disturb the anchor.

### Time guards

Recommended time behavior:

- cheap detection can run before the main pipeline;
- full edge topology should run only after cheap selector passes;
- high-resolution visual proxy should run only after the candidate is already much smaller than baseline;
- if time is near the final limit, skip candidate comparison and emit baseline;
- never delay `JD()` to compute a borderline structural branch.

### Source-size control

The current best source is close to the source cap. To fit this safely:

- implement only one recognizer per first submission unless you have already packed the source;
- periodic grid branch can reuse face triples, edge sorting, and output routines;
- AABB branch is the smallest and can be combined with periodic if byte budget allows;
- lat-long branch has the most logic and should probably be separate unless the source has been minified.

## Submission recommendations

### Submission 1: strict periodic `M == 2N` grid decimator

Include:

- `N` high selector;
- `M == 2N` selector;
- exact closed topology and `chi == 0` if cheap gates pass;
- degree-6 majority gate;
- factorization of `N = U*V`;
- row/column-major face-order proof;
- coordinate periodicity gate;
- decimated original grid candidates at two or three LODs;
- tournament against `1181` fallback.

Why first: this is the cleanest high-N structural hypothesis. If a hidden high-N case is a periodic torus/grid, this can be a very large score jump.

Risk: medium only because of source complexity and timing. False positives are low if face-order proof is mandatory.

### Submission 2: exact AABB cuboid shell

Include:

- high-N selector;
- bbox-plane occupancy;
- six-plane coverage;
- axis-normal gate;
- exact topology `chi == 2`;
- 8-vertex / 12-face cuboid candidate;
- optional `4x4` six-panel fallback only if source budget allows;
- tournament against `1181` fallback.

Why second: it is cheap, fail-closed, and can produce an enormous output reduction on dense cuboid shells.

Risk: low. The only real risk is over-loosening thresholds and accidentally hitting a rounded box.

### Submission 3: strict lat-long ring-sampled ellipsoid

Include:

- `M == 2N - 4`;
- closed topology `chi == 2`;
- two singleton poles;
- equal-size latitude rings;
- cap/strip face-order proof;
- ellipsoid radial residual;
- ring-sampled candidate LODs;
- tournament against fallback.

Why third: high upside, but more false-positive modes than periodic grid or exact AABB.

Risk: medium. Do not use residual-only detection.

### Submission 4: combined recognizer pack

Submit this only if individual branches were neutral but no failures occurred, and you have enough source-size slack. Combine:

1. periodic grid;
2. exact AABB;
3. lat-long ellipsoid.

Priority order inside the dispatcher should be:

1. `M == 2N` periodic grid;
2. exact AABB if bbox occupancy is overwhelming;
3. lat-long sphere/ellipsoid.

The dispatcher should not build multiple expensive candidates for the same input. Once a branch has a full structural proof, build its candidates and skip the other structural families.

## What not to try in this high-N push

1. Do not submit another broad high-N B16/S3B16 sweep. The target jump is too large, and broad high-N B16 has high catastrophic risk.
2. Do not use `N` bands alone. Every high-N branch must prove topology and structure.
3. Do not use `M == 2N` alone. It must be paired with `chi == 0`, closed manifold proof, degree-6 majority, and lattice face-order proof.
4. Do not use ellipsoid residual alone. A generic smooth shell can fit an ellipsoid well enough to fool a residual gate but still fail visual SSIM after analytic replacement.
5. Do not fire AABB on high axis-normal score alone. Require vertex-plane occupancy and per-triangle common-plane proof.
6. Do not output a tiny primitive before comparing with `1181`. Candidate must win the tournament and pass visual proxy.
7. Do not loosen visual proxy floors to chase vertex count. Falling below SSIM `0.9` makes the entire test score invalid.
8. Do not add arbitrary lattice reconstruction for the first periodic-grid submission. Require row/column-major face-order proof; arbitrary reconstruction is larger and riskier.
9. Do not add soft/rounded box recognition to the exact AABB submission. It is a different, higher-risk hypothesis.
10. Do not run full edge sorting on every high-N input. Only do it after cheap gates indicate a plausible structural case.
11. Do not combine new recognizers with changes to the existing `1181` route in the same submission. If the score changes, you need to know the cause.
12. Do not submit a branch that only improves `V_out` by a small amount. Structural candidates should deliver a large margin or no-op.

## Practical acceptance checklist

Before submitting any structural recognizer, confirm locally that all answers are yes:

- Does the source preserve exact `1181` behavior when the recognizer does not trigger?
- Is the branch impossible to trigger outside `60000..1100000`?
- Does it reject if topology cannot be proven?
- Does it reject if same-direction edge conflicts are nonzero?
- Does it reject if duplicate or degenerate faces exist?
- Does it reject if visual proxy cannot be computed in time?
- Does it compare against the `1181` fallback, not against a stale pre-pipeline estimate?
- Does it require a meaningful vertex-count win?
- Does it emit a valid closed oriented mesh with consistent face winding?
- Is the submission hash unique and not route-equivalent to previous plateau variants?

## Bottom line

The most promising high-N path is not more local collapse tuning. Add structural candidates that are separate from the `1181` simplifier and are accepted only after hard proof:

1. `M == 2N` periodic grid/torus: exact lattice proof, decimated original grid output.
2. Exact AABB shell: six-plane proof, 8-corner cuboid output.
3. Lat-long ellipsoid: two-pole equal-ring proof, ring-sampled or analytic ellipsoid output.

The first branch to submit should be the periodic `M == 2N` grid recognizer because it has the strongest topological signature and the largest plausible high-N upside. The AABB branch is the safest second submission. The lat-long branch is worth doing, but only with ring/pole and face-order proof; residual-only primitive fitting is too risky.
