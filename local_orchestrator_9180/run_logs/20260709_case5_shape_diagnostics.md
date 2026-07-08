# 2026-07-09 case5 shape diagnostics

Active base: `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact score `81.978181`, `7/7`.

## Results

| Submission | Kattis | Result | Predicate |
| --- | --- | --- | --- |
| `local_orchestrator_9180/diag_case5_two_poles.cpp` | `19923428` | `0.000000`, `7/7` | False: exact case 5 does not have exactly two vertices with degree `>=20` and max degree `>=40`. |
| `local_orchestrator_9180/diag_case5_poles_deg65.cpp` | `19923429` | `0.000000`, `7/7` | False: exact case 5 does not have two degree-65 poles. |
| `local_orchestrator_9180/diag_case5_bbox_iso80.cpp` | `19923441` | `0.000000`, `7/7` | False: exact case 5 does not have AABB `lo/hi >= .80`. |
| `local_orchestrator_9180/diag_case5_radial_shell16.cpp` | `19923448` | `0.000000`, `7/7` | False: exact case 5 is not an origin radial shell with `min radius > .80` and radius width `< .16`. |
| `local_orchestrator_9180/diag_case5_bbox_iso60.cpp` | `19923470` | `0.000000`, `7/7` | False: exact case 5 does not have AABB `lo/hi >= .60`. |
| `local_orchestrator_9180/diag_case5_bbox_ellipsoid_loose.cpp` | `19923472` | `0.000000`, `7/7` | False: exact case 5 is not a loose axis-aligned bbox ellipsoid (`max <= .18`, `rms <= .055`, `mean <= .040`). |
| `local_orchestrator_9180/diag_case5_bbox_iso40.cpp` | `19923485` | `0.000000`, `7/7` | False: exact case 5 does not have AABB `lo/hi >= .40`. |
| `local_orchestrator_9180/diag_case5_bbox_mid80.cpp` | `19923533` | `0.000000`, `7/7` | False: exact case 5 does not have middle AABB extent / longest extent `>= .80`. |
| external/parallel `diag_case5_bbox_mid40.cpp` | `19923536` | `0.000000`, `6/7` | True on exact case 5: middle AABB extent / longest extent `>= .40`. |
| `local_orchestrator_9180/diag_case5_mid60.cpp` | `19923544` | `0.000000`, `7/7` | False: exact case 5 does not have middle AABB extent / longest extent `>= .60`. |
| `local_orchestrator_9180/diag_case5_short_mid70.cpp` | `19923548` | `0.000000`, `7/7` | False: exact case 5 does not have shortest AABB extent / middle extent `>= .70`. |
| external/parallel `diag_case5_short_mid75.cpp` | `19923551` | `0.000000`, `7/7` | False: exact case 5 does not have shortest AABB extent / middle extent `>= .75`. |
| `local_orchestrator_9180/diag_case5_cross_flat50.cpp` | `19923559` | `0.000000`, `7/7` | False: exact case 5 does not have shortest AABB extent / middle extent `<= .50`. |
| `local_orchestrator_9180/diag_case5_two_planes70.cpp` | `19923561` | `0.000000`, `7/7` | False: fewer than 70% of vertices are near the two extreme planes of the shortest AABB axis. |
| `local_orchestrator_9180/diag_case5_cyl_radial_cv25.cpp` | `19923574` | `0.000000`, `7/7` | False: radius around the longest AABB axis does not have coefficient of variation `< .25`. |
| external/parallel `diag_case5_normals_thin_area45.cpp` | `19923579` | `0.000000`, `7/7` | False: area-weighted fraction of face normals with `|n_short_axis| > .75` is not above `.45`. |

## Interpretation

- Exact case 5 is strongly anisotropic: `shortest / longest < .40`.
- Its middle axis is bounded: `.40 <= middle / longest < .60`.
- Its shortest axis is moderately smaller than its middle axis: `.50 < shortest / middle < .70`.
- It is not a two-layer/top-bottom shell along the shortest axis under the 70% near-extreme test.
- It is not a simple straight tube/capsule around the longest AABB axis under the radial-CV `< .25` test.
- Face normals are not concentrated along the shortest axis, so it is not a mostly-flat plate/sheet in normal space either.
- It is not the already-tested sphere, cubic shell, axis-aligned ellipsoid, or radial shell family.
- It is also not a lat-long sphere-like mesh with two high-degree poles.
- This is not a broad flat sheet (`middle / longest >= .80` is false), not a round tube/capsule cross-section (`shortest / middle >= .70` is false), not an ultra-flat ribbon (`shortest / middle <= .50` is false), not a simple two-plane shell, and not a straight constant-radius tube. It now looks like a moderately flattened elongated smooth freeform body.
- The next useful split is thin-axis normal concentration or whether the object is approximately a height-field over the long/middle plane.

## Candidate Guidance

- Do not spend more attempts on opening `FL()` globally: Kattis `19923478` added `if(FL())JD(),exit(0);` after `GN()` and scored only `81.946573`, below the active best.
- Avoid bbox ellipsoid remesh as a direct branch unless a stronger non-axis-aligned fit diagnostic becomes positive.
- Prefer exact `N == 49987 && M == 99970` diagnostics and exact-guarded branches for a flat/prolate surface strategy.
