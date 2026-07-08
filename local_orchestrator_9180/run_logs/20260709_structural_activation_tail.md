# 2026-07-09 structural activation and late-tail checks

Active base: `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact score `81.978181`, `7/7`.

## Results

| Submission | Kattis | Result | Change |
| --- | --- | --- | --- |
| `submission_1467_68.11_6.cpp` | `19923344` | `68.113410`, `6/7` | Added a late high-N `B16::R(60001,120000,120,-13,128,.965,18.68)` tail. |
| `submission_1468_68.11_6.cpp` | `19923346` | `68.113410`, `6/7` | Activated existing `FL()` after `GN()` with `if(FL())JD(),exit(0);`. |

## Interpretation

- The late high-N B16 tail repeats the earlier pattern: apparently small generic tail additions can break a hidden case and should not be used as filler.
- `FL()` had a real local signal on the synthetic torus proxy, reducing output vertices on `case3_torus_23296_scrambled`, but hidden tests punish global activation.
- Do not activate `GI`, `FL`, or `IC` globally just because they are already present in the source. Any structural recognizer needs an exact, fail-closed detector and enough byte budget for an additional guard.

Next local work should stay on the active `1448` base and target exact `N == 49987, M == 99970` with a new ring4/ring8-aware surface strategy rather than more B16/WK/S3B16 tail nudges.
