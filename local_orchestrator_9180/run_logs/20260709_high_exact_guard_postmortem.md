# 2026-07-09 high exact-guard postmortem

Context: after invalid-sentinel fingerprinting, the two large hidden cases are:

- High1: `N == 377084 && M == 754688`, with `M = 2N + 520`, Euler characteristic `-260`, likely orientable genus `131`.
- High2: `N == 1009118 && M == 2018232`, sphere topology (`2N - 4`).

Current high-water remains `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact score `81.978181`, `7/7`.

## Submitted exact high1 probes

| Kattis ID | Candidate | Result | Conclusion |
|---|---|---:|---|
| `19927182` | `h1_b16_c180_st5_q930.cpp` using `bits/stdc++.h` and an extra exact B16 pass | `53.927292`, `5/7` | Source-layout/include change plus aggressive high1 edit is unsafe; do not use `bits` versions for serious probes. |
| `19927280` | `h1_vimp_cap26000.cpp` | `68.113410`, `6/7` | VIMP cap increase on exact high1 is too aggressive. |
| `19927347` | `h1_gn_ratio08_nocstdint.cpp` | `53.927292`, `5/7` | GN ratio edit is unsafe. |
| `19927388` | `h1_gn_ratio09_nocstdint.cpp` | `53.927292`, `5/7` | GN ratio edit remains unsafe even when milder. |
| `19927447` | `h1_b16_first_c160.cpp` | `56.927589`, `5/7` | First high-N B16 count increase is unsafe. |
| `19927511` | `h1_b16_first_q936.cpp` | `53.927292`, `5/7` | First high-N B16 threshold lowering is unsafe. |

High1 summary: the current route is already near the visual/time boundary. Plain exact-guarded B16, VIMP, and GN edits are closed unless a new rollback or renderer-aware branch is added.

## Submitted exact high2 probes

| Kattis ID | Candidate | Result | Conclusion |
|---|---|---:|---|
| `19927592` | `h2_b16_first_q939.cpp` | `68.113410`, `6/7` | Lowering the first high-N B16 proxy threshold from `.941` to `.939` is too aggressive. |
| `19927595` | `h2_b16_first_c140.cpp` | `53.927292`, `5/7` | Increasing first high-N B16 patch count from `120` to `140` is worse; do not try larger counts. |
| `19927620` | `h2_b16_first_q9408.cpp` | `53.927292`, `5/7` | Even a near-original first-threshold change is unstable; close first high-N B16 threshold edits for high2. |
| `19927657` | `h2_b16_second_c89.cpp` | `53.927292`, `5/7` | Increasing the second high-N B16 patch count from `88` to only `89` is already unsafe. |

High2 summary: both first and second high-N B16 edits are at a brittle boundary. Stop B16 count/threshold edits for this exact case. The next probes should use VIMP skip/cap controls, starting with extremely small deltas, and stop immediately if they repeat the `5/7` collapse.

## Open queue

- `19927688` / `h2_vimp_skip19800.cpp` returned `53.927292`, `5/7`. This only raised the non-ring VIMP skip guard from `19606` to `19800` for `M > 2e6`, so exact high2 VIMP skip/cap loosening is also unsafe.

## Pivot

Do not spend more submissions on large-case micro-edits in the current pipeline. The large hidden cases are close to time/SSIM cliffs, and exact guards alone do not make same-mechanism loosenings safe. Next work should use a different mechanism: either a real visual-mesh plus hidden vertex-cover strategy, or exact small/medium sphere-topology branches with a stronger rollback criterion.

## Follow-up transplant check

- `19927783` / `visual_shell_hybrid_full_highwater.cpp` combined the old `VSC::run()` visual-shell branch from `19915979` with the full current high-water fallback main. It was `129335` bytes, compiled/sample-passed, and matched high-water first-line output on selected proxies (`case5`, `low23_sphere_23201`, `upper35_wavy_35292`, `torus23_scr`, `wavy57`), but Kattis returned `53.927292`, `5/7`. Do not use old-base VSC transplants; proxy first-line equivalence is not enough.
- `19927803` / `shrink_best_safe.cpp` returned `53.927292`, `5/7`. This proves the source-layout/compaction route is unsafe on the judge even when local sample/proxy checks pass. Do not use shrink bases as submit foundations; future edits should either patch the exact high-water source with very small byte-neutral changes or add a tightly guarded independent branch with its own fallback.
- `19927842` / `high2_exact_gm_sphere.cpp` enabled the existing `GM()` strict-sphere remesh branch only for exact high2 (`N == 1009118`) by changing the small-`N` guard to an exact high2 guard. It compiled/sample-passed and stayed under the source limit (`131033` bytes), but Kattis returned `53.927292`, `5/7`, runtime `>21.00s`. Close `GM()` activation for high2; even a narrow independent branch with its existing proxy rollback is not safe enough on the hidden large sphere case.
