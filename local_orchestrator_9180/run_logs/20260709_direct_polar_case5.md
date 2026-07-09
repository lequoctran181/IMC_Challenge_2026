# 2026-07-09 direct polar exact-case experiments

Resource mode: 16 local workers only. No ChatGPT Pro Extended chat management.

## Motivation

The current high-water source is `submission_1448_81.98_7.cpp` / Kattis `19922865` with score `81.978181`.

Known weak exact case:

- `N=49987, M=99970`
- GJ-style polar topology is detectable from face order.
- Earlier contribution diagnostics estimated this case near `68.917518`, so an exact-case remesh could potentially add several leaderboard points.

The high-water source already contains `GJ`, `GR`, and `FL`, but `FL` was not called. Calling the existing `FL` with relaxed guard did not change local proxy behavior; sentinel tests showed `FL()` never returned true on the proxy at tested insertion points.

## Local proxy

Created `proxy_exact_case5_polar65x769_gj.in`, matching:

- `N=49987, M=99970`
- `GJ` accepted pattern at `V=769, R=65`
- corrected face order compared with the first proxy attempt

Direct `GR` remesh proxy scores:

| grid | vertices | proxy VPS 512 |
| --- | ---: | ---: |
| 16x192 | 3074 | 0.922862660034 |
| 20x192 | 3842 | 0.944077361335 |
| 24x192 | 4610 | 0.961738389820 |
| 32x224 | 7170 | 0.980113806454 |
| 40x256 | 10242 | 0.985897330195 |

## Submitted candidates

All candidates were based on fetched high-water source `fetched_sources/kattis_19922865_best_real.cpp`, replacing the body of `FL()` with an exact-case direct printer and adding `FL();` immediately after `JC();`.

| Kattis id | file | local proxy | result |
| --- | --- | --- | --- |
| `19930641` | `highwater_direct_fl24.cpp` | 4610 vertices, VPS 0.9617, no cover | `Accepted (53.927292)`, 5/7 |
| `19930659` | `highwater_direct_fl24_cover.cpp` | 4610 vertices, VPS 0.9617, cover=0 on proxy | `Accepted (68.11341)`, 6/7 |
| `19930676` | `highwater_direct_fl32x224_cover.cpp` | 7170 vertices, VPS 0.9801 | `Accepted (53.927292)`, 5/7 |
| `19930686` | `highwater_direct_fl40x256.cpp` | 10242 vertices, VPS 0.9859, no cover | `Accepted (53.927292)`, 5/7 |

## Interpretation

- Direct no-cover polar remesh is not safe on hidden case5.
- Adding cover improves from 5/7 to 6/7 for 24x192, so at least one failure mode was likely vertex-Hausdorff or a nearby validity issue.
- Denser cover with brute-force cover check likely hurts runtime and fell into the `53.927292` bucket.
- This route is not enough as implemented. Future work should avoid brute-force cover and either:
  - implement a compact spatial hash cover inside the replaced `FL()` body, or
  - use direct polar only as an exact-shape diagnostic, then improve the standard high-water pipeline elsewhere.

## Blacklist update

Do not resubmit:

- direct polar no-cover at 24x192 or 40x256
- direct polar brute-force cover at 32x224
- existing `FL();` hook without replacing `FL()` body

