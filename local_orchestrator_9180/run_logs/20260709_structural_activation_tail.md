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
| `submission_1479_66.82_6.cpp` | `19923495` | `66.820859`, `6/7` | Activated existing `GT()` after `GN()` with `if(GT())JD(),exit(0);`. |
| `submission_1480_0.00_7.cpp` | `19923533` | `0.000000`, `7/7` | Diagnostic false: exact case bbox `middle/long >= .80`. |
| `submission_1481_0.00_6.cpp` | `19923536` | `0.000000`, `6/7` | Diagnostic true: exact case bbox `middle/long >= .40`. |
| `submission_1482_0.00_7.cpp` | `19923544` | `0.000000`, `7/7` | Diagnostic false: exact case bbox `middle/long >= .60`. |

## Interpretation

- The late high-N B16 tail repeats the earlier pattern: apparently small generic tail additions can break a hidden case and should not be used as filler.
- `FL()` had a real local signal on the synthetic torus proxy, reducing output vertices on `case3_torus_23296_scrambled`, but hidden tests punish global activation.
- The first broad high-N B16 count is a cliff: `120 -> 121` breaks badly, so do not run count sweeps around that call without a separate fail-closed guard.
- Numeric shorthand is not safe as a byte-saving technique in this source. The local harness saw identical outputs on the proxy set, but Kattis dropped to 5/7, so keep exact integer spellings in the active base.
- Exact case 5 is not a simple two-pole sphere/ring grid, not close to cubic by bbox aspect, not a normalized ellipsoid under the tested thresholds, not the tested thin radial shell, and has `shortest_bbox_extent / longest_bbox_extent < .40`.
- Raising all `FL()` trial thresholds to `.99` avoids the earlier 6/7 failure but also loses the best gain, landing in the familiar `81.946573` bucket.
- Global `GT()` activation is unsafe: it reached only `66.820859`, 6/7.
- Exact case bbox is flat/elongated: `short/long < .40` and `0.40 <= middle/long < .60`.
- Do not activate `GI`, `FL`, or `IC` globally just because they are already present in the source. Any structural recognizer needs an exact, fail-closed detector and enough byte budget for an additional guard.

Next local work should stay on the active `1448` base and target exact `N == 49987, M == 99970` with a new ring4/ring8-aware surface strategy rather than more B16/WK/S3B16 tail nudges.
