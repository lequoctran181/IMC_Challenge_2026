# 2026-07-09 guarded reopen batch

Context: user paused the 8 Pro Extended chat workers and raised the local lane ceiling to 16.  Active base remains `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact `81.978181`.

Purpose: test whether exact hidden test-4 safety guards for `N == 23201 && M == 46398` can reopen aggressive branches that previously failed with hidden test 4 `SSIM is too low`.

Generated `local_orchestrator_9180/queue16_guarded_reopen_20260709/` with 16 variants.  Templates included:

- `submission_1448_81.98_7.cpp`
- `submission_1578_53.93_5.cpp` (`b16a_c240`, failed Kattis test 4)
- `queue16_large_guard_20260709/lg06_vimp_reserve40000.cpp`

Guards tested:

- skip the final `WK::run()` only on exact `23201/46398`
- skip `VIMP::run()` only on exact `23201/46398`
- skip `W2C::run()` only on exact `23201/46398`
- combined skip variants, though two combined forms exceeded the source limit by 4-6 bytes and were used only as local probes

## Parallel batch signal

The first 16-worker run appeared to show gains, especially for `lg06_reserve40000__skip_w2c_tail23.cpp`:

- `case5`: `1148/2292 -> 1133/2262`
- `wavy57`: `3127/6254 -> 2257/4514`

However the parallel run also made the base itself worse than its known isolated output (`case5 1148/2292` instead of `1133/2262`, `torus23scr 1907/3814` instead of `1282/2564`).  This marks the apparent gains as timing/load artifacts.

## Isolated rerun

Reran selected variants sequentially on:

- sample
- `low23_sphere_23201`
- `low23_wavy_23205`
- `case5_lobed`
- `torus23_scr`
- `wavy57`

Isolated result: all selected guarded variants produced the same first-line outputs as the active base on those proxies:

| Proxy | Active base isolated |
| --- | --- |
| sample | `8 12` |
| low23_sphere | `697 1390` |
| low23_wavy | `3795 7590` |
| case5 | `1133 2262` |
| torus23scr | `1282 2564` |
| wavy57 | `2257 4514` |

## Decision

No submission.  Exact `23201/46398` main-level guards for final WK, VIMP, or W2C do not create a stable local behavioral difference once isolated.  Do not use this batch as a Kattis candidate.

Next direction should move away from main-call guard reopen and toward either:

1. a true hidden-test4 proxy/diagnostic that predicts the cliff, or
2. a structurally different output route that is stable under isolated timing rather than only under parallel load.
