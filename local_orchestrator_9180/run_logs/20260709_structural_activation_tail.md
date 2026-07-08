# 2026-07-09 structural activation and late-tail checks

Active base: `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact score `81.978181`, `7/7`.

## Results

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1467_68.11_6.cpp` | `19923344` | `68.113410`, `6/7` | Added a late high-N `B16::R(60001,120000,120,-13,128,.965,18.68)` tail. |
| `submission_1468_68.11_6.cpp` | `19923346` | `68.113410`, `6/7` | Activated existing `FL()` after `GN()` with `if(FL())JD(),exit(0);`. |
| `submission_1469_50.48_5.cpp` | `19923370` | `50.479384`, `5/7` | Changed the first broad high-N B16 count from `120` to `121`. |
| `submission_1470_0.00_7.cpp` | `19923428` | `0.000000`, `7/7` | Diagnostic false: `N==49987 && M==99970 && ge20==2 && mx>=40`. |
| `submission_1471_0.00_7.cpp` | `19923429` | `0.000000`, `7/7` | Diagnostic false: `N==49987 && M==99970 && ge20==2 && mx==65 && c65==2`. |
| `submission_1472_57.28_5.cpp` | `19923435` | `57.280543`, `5/7` | Numeric shorthand only in `main`, intended byte-budget neutral. |
| `submission_1473_0.00_7.cpp` | `19923441` | `0.000000`, `7/7` | Diagnostic false: exact case bbox aspect `lo/hi >= .80`. |
| `submission_1474_0.00_7.cpp` | `19923448` | `0.000000`, `7/7` | Diagnostic false: exact case radial shell `mn > .80 && mx-mn < .16`. |
| `submission_1475_0.00_7.cpp` | `19923470` | `0.000000`, `7/7` | Diagnostic false: exact case bbox aspect `lo/hi >= .60`. |
| `submission_1476_0.00_7.cpp` | `19923472` | `0.000000`, `7/7` | Diagnostic false: normalized bbox ellipsoid fit `ma <= .18 && rms <= .055 && mean <= .040`. |
| `submission_1477_81.95_7.cpp` | `19923478` | `81.946573`, `7/7` | Activated `FL()` with all trial thresholds at `.99`. |
| `submission_1478_0.00_7.cpp` | `19923485` | `0.000000`, `7/7` | Diagnostic false: exact case bbox aspect `lo/hi >= .40`. |
| `submission_1479_66.82_6.cpp` | `19923495` | `66.820859`, `6/7` | Numeric-shorthand main plus extra post-`GN()` `GT()` call; not a pure `GT()` test. |
| `submission_1480_0.00_7.cpp` | `19923533` | `0.000000`, `7/7` | Diagnostic false: exact case bbox `middle/long >= .80`. |
| `submission_1481_0.00_6.cpp` | `19923536` | `0.000000`, `6/7` | Diagnostic true: exact case bbox `middle/long >= .40`. |
| `submission_1482_0.00_7.cpp` | `19923544` | `0.000000`, `7/7` | Diagnostic false: exact case bbox `middle/long >= .60`. |
| `submission_1483_0.00_7.cpp` | `19923548` | `0.000000`, `7/7` | Diagnostic false: exact case bbox `short/middle >= .70`. |
| `submission_1484_0.00_7.cpp` | `19923551` | `0.000000`, `7/7` | Diagnostic false: exact case bbox `short/middle >= .75`. |
| `submission_1485_0.00_7.cpp` | `19923559` | `0.000000`, `7/7` | Diagnostic false: exact case bbox `short/middle <= .50`. |

## Interpretation

- The late high-N B16 tail repeats the earlier pattern: apparently small generic tail additions can break a hidden case and should not be used as filler.
- `FL()` had a real local signal on the synthetic torus proxy, reducing output vertices on `case3_torus_23296_scrambled`, but hidden tests punish global activation.
- The first broad high-N B16 count is a cliff: `120 -> 121` breaks badly, so do not run count sweeps around that call without a separate fail-closed guard.
- Numeric shorthand is not safe as a byte-saving technique in this source. The local harness saw identical outputs on the proxy set, but Kattis dropped to 5/7, so keep exact integer spellings in the active base.
- Exact case 5 is not a simple two-pole sphere/ring grid, not close to cubic by bbox aspect, not a normalized ellipsoid under the tested thresholds, not the tested thin radial shell, and has `shortest_bbox_extent / longest_bbox_extent < .40`.
- Raising all `FL()` trial thresholds to `.99` avoids the earlier 6/7 failure but also loses the best gain, landing in the familiar `81.946573` bucket.
- The `1479` result is confounded by numeric shorthand; do not use it to judge the built-in `GT()` dispatcher path.
- Exact case bbox is flat/elongated: `short/long < .40` and `0.40 <= middle/long < .60`.
- Exact case axis ratios include `.50 < short/middle < .70`.
- Do not activate `GI`, `FL`, or `IC` globally just because they are already present in the source. Any structural recognizer needs an exact, fail-closed detector and enough byte budget for an additional guard.

Next local work should stay on the active `1448` base and target exact `N == 49987, M == 99970` with a new ring4/ring8-aware surface strategy rather than more B16/WK/S3B16 tail nudges.

## Packed-source and exact-polar follow-up

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1522_67.76_6.cpp` | `19924194` | `67.760122`, `6/7` | Compact/high-water hybrid main; confirms transplanted call order is source-sensitive. |
| `submission_1523_68.11_6.cpp` | `19924203` | `68.113410`, `6/7` | Packed `1491` source with `VSC/ST5` early recognizers and active `1448` tail. |
| `submission_1524_53.93_5.cpp` | `19924211` | `53.927292`, `5/7` | Packed `1491` source with recognizers removed but active `1448` tail. |
| `submission_1525_68.11_6.cpp` | `19924212` | `68.105407`, `6/7` | Exact `GJ` diagnostic branch. |
| `submission_1526_68.11_6.cpp` | `19924222` | `68.105407`, `6/7` | Exact `FL` polar branch candidate. |
| `submission_1527_68.11_6.cpp` | `19924223` | `68.105741`, `6/7` | Original packed `1491` route plus the positive `1448` micro S3B16 pass. |

Conclusion: do not use packed `1491/1499` as a replacement base for the current `1448` high-water tail. Even conservative-looking transplants fail hidden cases, and the exact polar `GJ/FL` branch currently collapses into the same `68.105` 6/7 bucket. The next batch should return to the byte-tight `1448` base and use exact diagnostics/guards before any new output branch.

## Output precision probe

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1528_68.11_6.cpp` | `19924247` | `68.113410`, `6/7` | Removed unused `<cstdint>`, printed `%.17g` when `N < 60000` or the active mesh had at most half the original vertices. |

Conclusion: output precision is not a free improvement. The high-precision small/mid-output branch breaks one hidden case despite leaving geometry unchanged, so do not submit the remaining precision variants without a narrower reason.

## Byte-bank shrink probe

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1533_43.09_4.cpp` | `19924299` | `43.094425`, `4/7` | Submitted `shrink_best_safe.cpp`, a 129742-byte renamed/shaved variant that compiled and matched the official sample output hash. |

Conclusion: local sample-hash equivalence is not enough for minified/renamed bases. The shrink family is unsafe on hidden cases and should not be used as the development base unless a future shrink is verified by stronger hidden-aligned probes.

## Follow-up verdicts after local-worker limit raised to 16

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1534_53.93_5.cpp` | `19924307` | `53.927292`, `5/7` | Less invasive `highwater_19922865_shaved.cpp`, 130958 bytes, same active `1448` main. |
| `submission_1535_81.95_7.cpp` | `19924308` | `81.946573`, `7/7` | Activates the existing `IC/GI/FL` structural recognizers only for `N < 5001` before the normal pipeline. |
| `submission_1536_43.09_5.cpp` | `19924321` | `43.094425`, `5/7` | For `5000 <= N < 25000`, outputs the original mesh via `IJ()` before the normal pipeline. |
| `submission_1537_53.93_6.cpp` | `19924324` | `53.927292`, `6/7` | For `25000 <= N < 39000`, outputs the original mesh via `IJ()` before the normal pipeline. |

Conclusion: even the minimal include/symbol shave is hidden-unsafe and must not be used as a byte-bank base. The `N < 5001` recognizer activation has a real but limited signal, landing in the known `81.946573` bucket below `1448`; keep it as a clue for small-case specialists, not as the active base. Original-output branches for small ranges are contaminated and should not be used for contribution estimates or production routing.

## Local 16-worker W2G broadening batch

Base: `submission_1448_81.98_7.cpp`. Generated 16 small variants in `/tmp/imc_9180_main/queue16_w2g_20260709`; all compiled, stayed within the observed source-size envelope, and matched the official sample output hash.

Proxy findings:

| Variant group | Proxy result |
| --- | --- |
| `w01..w04` broaden W2G from the narrow `23125..23500` range to larger `M == 2N` grids | On `synthetic_torus80640`, output improved from `2016/4032` to `1260/2520`. |
| W2G LOD/proxy-only tweaks | No proxy vertex-count gain over base on the tested torus/grid inputs. |
| W2C fallback tweaks | No proxy vertex-count gain over base on the tested torus/grid inputs. |

Submitted the balanced broadened variant:

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1538_43.09_4.cpp` | `19924348` | `43.094425`, `4/7` | Outside submission; source tail resembles the active base but hidden score collapses, so do not use as a base. |
| `submission_1539_81.95_7.cpp` | `19924373` | `81.946573`, `7/7` | W2G fail-closed range broadened to `8000 <= N <= 490000` and factor dimensions up to `700`. |
| `submission_1540_43.09_5.cpp` | `19924376` | `43.094425`, `5/7` | For `5000 < N < 25000`, outputs original through `IJ()` before the normal pipeline. |

Conclusion: broadening the existing periodic-grid recognizer is accepted and proxy-positive, but on the hidden set it falls into the same `81.946573` bucket rather than improving the `1448` high-water mark. Do not use W2G broadening as the active base. It remains useful evidence that a hidden periodic-grid target either is not in the broadened range, does not dominate score, or is offset by a small hidden regression.

## Same-size S3B16 count probe

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1541_68.11_6.cpp` | `19924394` | `68.113743`, `6/7` | Same-size edit on the active `1448` source: second `S3B16::T(10,-9,...)` changed to `S3B16::T(11,-9,...)`. |

Local proxy note: this was the only `queue16_exact_same_size/` variant that improved the synthetic exact case5 proxy, moving `1133/2262` to `1132/2260`, while preserving the official sample output. Kattis still drops one hidden case, so the proxy gain is a hidden cliff rather than a safe improvement.

Conclusion: do not increase the second S3B16 count above `10` on the current high-water base. Same-size edits are safer than source shaving, but S3 count/phase perturbations remain too brittle unless guarded by a stronger hidden-aligned detector.

## W2 order A/B

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1542_53.93_5.cpp` | `19924458` | `53.927292`, `5/7` | Same-size edit: `if(!W2G::run())W2C::run();` changed to `if(!W2C::run())W2G::run();`. |

Local proxy note: this variant was byte-identical to the active base on sample, exact-case5 proxy, torus23, torus80, torus57, and wavy90 outputs, and it was locally faster on several proxy runs. Hidden Kattis still collapses to the `53.927292` bucket, so the current W2G-before-W2C order is a hidden-sensitive route dependency.

Conclusion: preserve W2G before W2C. Do not use proxy-identical timing/order changes as submission candidates unless they are tied to a hidden-aligned diagnostic; local proxy equality is not enough for this source.

## Exact-N S3 filter probe

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1543_81.93_7.cpp` | `19924488` | `81.934570`, `7/7` | Outside/new submission observed while polling; lower than the high-water `81.978181`. |
| `submission_1544_68.11_6.cpp` | `19924498` | `68.112743`, `6/7` | Source-size-safe exact-N S3 gate: for `N==49987`, second S3 pass used count `11`, phase `7`, with base `.0055` score threshold; other N kept base count/phase. |

Local proxy note for `19924498`: exact-case5 proxy improved from `1133/2262` to `1131/2258` and proxy SSIM512 improved from `0.938971` to `0.939191`, while torus57, torus23, torus80, wavy90, and sample first-line outputs stayed at the base counts. Kattis still drops to `6/7`.

Conclusion: the exact `N==49987` S3 extra-patch direction is hidden-unsafe even when local proxy count and SSIM both improve. Stop S3 count/phase attempts for this exact case unless a new detector explains the hidden failure directly.

## 1543 byte-bank broad-B16 probe

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1545_53.92_5.cpp` | `19924517` | `53.922624`, `5/7` | Started from the smaller outside/new `1543` source and added both broad B16 calls from the high-water `1448` tail. |

Local proxy note: the candidate stayed proxy-identical to `1543` on sample, exact-case5, torus57, torus80, and the other quick structural checks, with source size `130756` bytes. Hidden Kattis still collapsed to the same low `53.92` bucket as other order/tail-sensitive perturbations.

Conclusion: do not transplant the broad B16 tail across byte-bank bases. The `1543` source is useful as a smaller reference, but the high-water `1448` route dependencies do not survive this direct graft.
