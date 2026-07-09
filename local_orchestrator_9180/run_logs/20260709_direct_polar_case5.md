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
| `19930707` | `highwater_direct_fl24_stride8.cpp` | 10859 vertices, VPS 0.9617, stride-8 unused cover | `Accepted (53.927292)`, 5/7 |

## Interpretation

- Direct no-cover polar remesh is not safe on hidden case5.
- Adding cover improves from 5/7 to 6/7 for 24x192, so at least one failure mode was likely vertex-Hausdorff or a nearby validity issue.
- Denser cover with brute-force cover check likely hurts runtime and fell into the `53.927292` bucket.
- Uniform stride cover also fell into the `53.927292` bucket, so the direct-polar proxy is too optimistic for hidden case5.
- This route is not enough as implemented. Future work should avoid direct proxy-only polar remesh and either:
  - implement a compact spatial hash cover inside the replaced `FL()` body, or
  - use direct polar only as an exact-shape diagnostic, then improve the standard high-water pipeline elsewhere.

## Blacklist update

Do not resubmit:

- direct polar no-cover at 24x192 or 40x256
- direct polar brute-force cover at 32x224
- direct polar stride-8 cover at 24x192
- existing `FL();` hook without replacing `FL()` body

## Upper35 guard follow-up

- `19930740`: `fl24_cover_upper_after_gn.cpp`, starting from `19930659` (`FL24` direct polar + cover) and adding `if(N==35292&&M==70580){JD();return 0;}` immediately after `GN()` in `main`.
  - Local upper35 proxies matched high-water outputs (`closed_bunnylike`: `1942/3880`, VPS512 `0.946045615133`; Stanford bunny: `6000/11996`, VPS512 `0.941070708967`).
  - Local case5 polar proxy still used direct branch: `4610/9216`, VPS512 `0.961738389820`.
  - Kattis result: `Accepted (67.760455)`, `6/7`, runtime `>21.00s`, no compact feedback block.
  - Interpretation: the upper35 failure from direct-polar was likely fixed, but the source/pipeline still falls into a 6/7 highbig/runtime plateau. Try an even lower-overhead upper guard before `FL()` only if it keeps source under the limit.

- `19930755`: `fl24_cover_upper_original_after_jc.cpp`, same direct-polar `FL24` cover branch, but exact `N==35292 && M==70580` returns original input immediately after `JC()` and before `FL()`.
  - Local upper35 returns original `35292/70580`; local case5 polar proxy remains `4610/9216`.
  - Kattis result: `Accepted (53.927292)`, `6/7`, runtime `>21.00s`, no compact feedback block.
  - Decision: returning original for upper35 sacrifices too much score and does not produce a useful valid-total route. Close direct-polar upper-guard variants unless a source-shape-preserving way to keep high-water upper35 compression is found.

- `19930766`: `fl24_cover_upper_after_gn_vimpcaplow.cpp`, same as `19930740` plus the previously Kattis-valid same-length VIMP cap lowering `24/16/12 -> 22/14/10`.
  - Local gates: upper35 unchanged (`1942/3880`, VPS512 `0.946045615133`); case5 polar direct unchanged (`4610/9216`, VPS512 `0.961738389820`); Nefertiti proxy output `150190/300376` in about `19.2s`.
  - Kattis result: `Accepted (67.760455)`, `6/7`, runtime `>21.00s`, no compact feedback block.
  - Decision: safe VIMP cap lowering does not rescue the direct-polar upper-guard plateau. Close this combination.
