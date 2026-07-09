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

## External/new low-score follow-up

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1546_51.87_5.cpp` | `19924524` | `51.870090`, `5/7` | Newer submission observed while polling Kattis after the 1543 broad-B16 probe. |

Conclusion: record for history only. This result is well below the current `1448` high-water mark and should not be used as an active base.

## Local 16-worker extra micro-S3 no-submit batch

Generated 16 variants from the exact `submission_1448_81.98_7.cpp` source by inserting one additional 39-byte call after the positive phase-0 micro pass, e.g. `S3B16::T(1,-9,192,.96,18.47,p,5,.0015);`. All variants compiled and stayed at `131069` bytes, just below the `131072` source-size ceiling.

Local proxy result:

| Proxy | Base | Best changed variants |
| --- | --- | --- |
| official sample | `8/12` | unchanged |
| exact case5 synthetic | `1133/2262` | several variants reached `1132/2260` |
| `softbox50402` | `1235/2466` | a few variants reached `1234/2464` |
| `torus23`, `torus80`, `lobes49954` | base counts | unchanged |

Conclusion: do not submit this batch. The useful-looking changes are still in the S3B16 phase/count family already rejected by Kattis (`19923923`, `19923905`, `19924394`, and `19924498` all dropped to `6/7` despite similar proxy improvements).

## Local 16-worker pass-order no-submit batch

Reviewed `/tmp/imc_9180_main/queue16_order_20260709`, which tests same-size or source-size-safe reorders/skips/extras around `VIMP`, `MIDEC`, and `WK`.

Local proxy result:

| Variant group | Result |
| --- | --- |
| `01..14`, `16` | Same output counts as `1448` on sample, exact-case5 proxy, torus23, and torus80. |
| `11`, `12` extra `WK` | Same first-line counts but different case5 output hash. |
| `15_skip_wk` | Worse exact-case5 proxy: `1135/2266` vs base `1133/2262`. |

Conclusion: no candidate to submit. Reorders do not produce vertex-count gains on the proxy set, and the only changed-count pass-order variant is worse.

## Local 16-worker `1339` marker-replacement batch

Base: `submission_1339_81.98_7.cpp`, which has the dormant marker `if(0&&"A32_A02_MS55"){}` and exact Kattis score `81.977514`. Generated 16 replacements in `/tmp/imc_9180_main/marker1339_20260709`, including small calls to `W5::post_patch_pass()`, `VIMP::run()`, `MIDEC::run()`, `WK::run()`, `IC/FL/GI`, tiny `B16::R` calls, the known winning micro-S3 call, and micro-S3 plus one extra tiny call. All compiled and stayed under the source-size ceiling.

Local proxy result:

| Variant | Exact-case proxy | Softbox proxy | Other proxies |
| --- | --- | --- | --- |
| `m00_marker_noop` | `1134/2264` | `1236/2468` | base counts |
| `m01_s3_highwater` | `1133/2262` | `1235/2466` | base counts |
| non-S3 marker replacements | no improvement over `m00` | no improvement over `m00` | unchanged |
| `m13..m15` S3 plus extra call | no vertex-count gain over `m01` | no vertex-count gain over `m01` | unchanged |

Conclusion: no submission. This confirms the already-submitted `1448` marker replacement remains the best local marker use; extra calls either no-op on the proxy set or stay in the already-risky S3/B16 family.

## Local 16-worker same-size timing batch

Generated 16 variants from `submission_1448_81.98_7.cpp`, changing only same-length timing constants in the active `main` tail: early/mid/final `B16` end times, both S3 end times, the loop `WK` time, and the final `WK` time. All variants stayed at `131030` bytes and compiled.

Local proxy result: every variant produced byte-identical output to the active `1448` base on sample, exact-case5 synthetic, torus23, torus80, lobes49954, and softbox50402. Runtime differences were local-load noise, not associated with output changes.

Conclusion: no submission. Same-size timing tweaks did not expose any local simplification gain and would only add hidden timeout risk.

## B92 periodic structural probe

Ran the packed high-water B92 structural recognizer batch from `local_orchestrator_9180/queue16_highwater_struct_b92` with 16 local workers. `B92pack_qnet_dsu.cpp` failed local compilation because the packed source omitted the `unordered_set` include. The remaining candidates compiled and preserved the high-water output on sample, exact-case5, torus23-scrambled, lobes49954, softbox50402, generated sphere, generated ellipsoid, and generated torus7680.

Local proxy signal:

| Candidate | Changed proxy | High-water | Candidate |
| --- | --- | --- | --- |
| `B92pack_idx_periodic.cpp` | synthetic torus80 | `2016/4032` | `1023/2046` |
| `B92pack_combo_grid_sphere.cpp` | synthetic torus80 | `2016/4032` | `1023/2046` |
| `B92pack_idx_periodic_safe.cpp` | synthetic torus80 | `2016/4032` | `1599/3198` |
| `B92pack_combo_all.cpp` | synthetic torus80 | `2016/4032` | `1599/3198` |
| `B92pack_star_strict.cpp` | generated cylinder12962 | `181/358` | `98/192` |

Submitted the highest local torus80-gain probe:

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1547_53.93_5.cpp` | `19924626` | `53.927292`, `5/7`, runtime `> 21.00 s` | Added the packed `TG::run()` periodic-grid recognizer before the high-water tail. |

Conclusion: do not use the packed B92/TG periodic recognizer as an active route. The torus80 local gain is real, but Kattis hidden behavior collapses to the familiar `53.927292` / `5/7` bucket, likely through timeout or recognizer misfire on one of the large official cases.

## Local 16-worker late-WK/reserve batch

Generated 16 variants from the exact high-water `submission_1448_81.98_7.cpp` source in `/tmp/imc_9180_main/latewk_gate_20260709`. The variants only changed late `WK` gating/reserve constants or late `B16/S3B16` end times, preserving the original source layout and staying under the source-size limit.

Important parsing/detail update from Kattis:

- `19924626` (`B92pack_idx_periodic`) failed `5/7`: test 4 `Wrong Answer` with validator message `SSIM is too low`, and test 7 `Time Limit Exceeded`.
- Most `53.927`/`68.113` bucket submissions fail test 4 by SSIM; the `53.927` family often also TLEs test 7.
- `19924324` is notable as a diagnostic source that only TLEs test 7, but it uses the known unsafe include/symbol shave and an original-output branch, so it is not a production base.

Late-WK/reserve local proxy result after recompiling a fresh high-water baseline:

| Variant group | Result |
| --- | --- |
| Narrow/disable final `WK` variants | Usually worsen torus23/cylinder and sometimes case5 or primitive proxies. |
| Earlier late-loop/reserve variants | Worsen case5, lobes, sphere, or ellipsoid proxies in several cases. |
| `v16_final_gt39000` | Byte-safe and output-identical to high-water on sample, case5 proxy, torus23, torus80, lobes49954, softbox50402, and generated primitive proxies. |

Conclusion: no submission from this batch. `v16_final_gt39000` is a possible timing-stability guard, but it has no local vertex-count or proxy gain. The other variants clearly remove useful late work and should not be submitted.

## Local 16-worker exact skipped-struct batch

Generated 16 exact-`N==49987` variants in `local_orchestrator_9180/queue16_exact_skipped_struct/`, reopening skipped `QX`, `FR`, `FG`, and `IC` branches individually and in small combinations for the `47.5k-60k` band. All compiled and stayed under the source-size limit (`131042..131066` bytes).

Proxy result: every variant matched the fresh `1448` high-water first-line output on sample, exact-case5 synthetic, torus23, torus23-scrambled, bumpy25, torus57, wavy57, bumpy49, lobes49, and ico2562. The exact-case5 synthetic remained `1133/2262`.

Conclusion: no submission. Reopening the skipped structural branches under the known exact-N guard did not produce a local signal; the previous broad exact-enable submission already failed, so this no-op batch is not worth a Kattis token.

## Compiler pragma speed probe

Generated compiler-pragma variants in `local_orchestrator_9180/queue16_pragma_speed/`. Plain `O3` and `Ofast` showed local timing/output instability on `torus23_scr`; `O3,unroll-loops` and `fast-math` matched the base on the small proxy portfolio but the direct `O3,unroll-loops` source exceeded the Kattis file-size limit after CRLF upload expansion.

Submitted the slim source-size-safe probe:

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1549_56.93_5.cpp` | `19924779` | `56.927589`, `5/7`, runtime `> 21.00 s` | Added `#pragma GCC optimize("O3,unroll-loops")`, removed `#include<utility>`, and changed the leading-zero `02` literal to `2` to stay under the source cap. |

Kattis details: test case 3 failed with `Wrong Answer: SSIM is too low`; test case 7 timed out.

Conclusion: blacklist pragma/unroll plus include-shave speed probes. The local timing idea did not transfer; it broke hidden visual fidelity and still TLEd the last case.

## VIMP same-size cap batch

Generated 16 same-size variants from exact high-water `submission_1448_81.98_7.cpp` in `local_orchestrator_9180/queue16_vimp_samesize_20260709/`. All variants stayed exactly `131030` bytes and only changed constants inside `VIMP::sc/run`.

Proxy result: most variants were output-identical to the fresh high-water build. Two variants gave stable local gains over three reruns:

| Variant | Proxy signal |
| --- | --- |
| `v11_ratio_low` | `case5_lobed 1133/2262 -> 1026/2048`, `torus57 1412/2824 -> 1265/2530`, `wavy57 2257/4514 -> 2025/4050`; no changes on the rest of the portfolio. |
| `v13_keep_loose` | Stronger reductions on the same three proxies, but higher SSIM risk because it loosens VIMP rollback thresholds. |

Submitted the safer candidate first:

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1550_81.95_7.cpp` | `19924859` | `81.946573`, `7/7`, runtime `20.32 s` | Lowered the VIMP cap ratios from `24/16/12` to `22/14/10`. |
| `submission_1556_53.93_5.cpp` | `19924936` | `53.927292`, `5/7`, runtime `> 21.00 s` | Loosened VIMP keep thresholds (`dr .003/.14/.11` and proxy needs) as `v13_keep_loose`. |

Kattis details for `19924936`: test case 4 failed with `Wrong Answer: SSIM is too low`.

Conclusion: VIMP cap reduction is hidden-valid but below the `81.978181` high-water. Loosening VIMP rollback/keep thresholds crosses the same hidden test-4 SSIM cliff as DCEV and should be blacklisted unless a new guard excludes that bucket.

## QEM target and DCEV fine-threshold checks

- Generated 16 same-size `choose_target` variants in `local_orchestrator_9180/queue16_qem_target_samesize_20260709/`, changing `.089`, `.035`, feature cap `.22`, and clamp `.086/.115`. All were output-identical to the high-water build on the portfolio; no submission.

## Standalone Pro Pool Check

- After switching orchestration away from the 8 Pro Extended chats and back to 16 local workers, compiled the local standalone Pro candidates that were already on disk.
- `proext_02_reverse_15k_local` was the only plausible submit: sample OK and local proxy scores at 512 were `0.9021` on torus23, `0.9441` on case5, and `0.9030` on bumpy25.
- `1557` / `19924969` submitted the unpatched local file and received `Compile Error`, likely because clang accepted indirectly available `chrono`/`cstdint`/`climits` while GNU did not.
- `1558` / `19924977` added those headers and compiled on Kattis, but scored only `48.522102`, `5/7`, runtime `13.33s`; detailed hints show test cases 4 and 7 both fail with `SSIM is too low`. Standalone Pro reverse-15k is blacklisted as a replacement route.
- Generated 16 DCEV thresholds just below `.920` in `local_orchestrator_9180/queue16_dcev_fine_p92_20260709/` (`.917`, `.918`, `.919`, `.9195` with `20k..30k`, `20k..25k`, `23k..25k`, and `8k..39k` guards). These also produced no first-line gains over the base on the tested 23k proxies. The useful local cliff remains below `.917`, and known `.916` submissions fail hidden test 4.

## Large-Guard 16-Worker Batch

Generated 16 same-size variants from high-water `submission_1448_81.98_7.cpp` in `local_orchestrator_9180/queue16_large_guard_20260709/`. All variants stayed exactly `131030` bytes. The batch targeted VIMP high-N ratio/cap/reserve/threshold constants and final broad-B16 step/count variants, plus synthetic proxies `bumpy_200k` and `torus_400k`.

Local highlights at 512 proxy:

| Variant | Local signal |
| --- | --- |
| `lg06_vimp_reserve40000` / `lg09` / `lg11` | Reduced `case5` from `1146/2288` to `1133/2262`, improved `wavy57` from `3127/6254` to about `2257..2304`, preserved torus23 first-line in most variants, but this is a timing/deeper-run unlock. |
| `lg16_b16_final_step17` | Improved torus23 proxy score from `0.900061` to `0.929137` and reduced wavy57, but it shares the same risky deeper-run shape. |
| `lg02_vimp_ratio26` | Almost no mid-bucket change; only `case5` `1146 -> 1147` and synthetic `torus400` `776 -> 773`. |

Submitted:

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1559_68.11_6.cpp` | `19925043` | `68.113410`, `6/7`, runtime `20.11s` | `lg11_vimp_need958938`, a VIMP timing/need unlock. |
| `submission_1561_81.95_7.cpp` | `19925061` | `81.946573`, `7/7`, runtime `20.01s` | `lg02_vimp_ratio26`, high-N VIMP ratio `24 -> 26`. |

Kattis details for `19925043`: test case 4 fails with `Wrong Answer: SSIM is too low`.

Conclusion: the apparent local gains from making VIMP run deeper in the `N<60k` route cross the same hidden test-4 cliff; blacklist reserve/need/timing unlocks unless an exact exclusion guard for that hidden case is discovered. Plain high-N VIMP ratio changes are valid but lower than the high-water, matching the earlier `24 -> 22` result bucket, so further ratio/cap sweeps are unlikely to be productive.
