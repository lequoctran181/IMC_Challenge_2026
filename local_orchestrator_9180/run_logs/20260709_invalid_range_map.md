# 2026-07-09 invalid-sentinel range map

These diagnostics used tiny echo-or-invalid sources generated under
`local_orchestrator_9180/invalid_diags/`. They do not run the high-water
pipeline; if the predicate is true they print an invalid `0 0` mesh, otherwise
they echo the original input. Therefore the accepted test count is a cleaner
predicate-hit signal than earlier full-pipeline fail-closed probes.

| Submission | Predicate | Result | Interpretation |
|---|---:|---:|---|
| `19925091` | `10 < N < 5000` | `Accepted (0)`, `6/7` | Exactly one non-sample small case in this range. |
| `19925787` | `N >= 11 && N < 2500` | `Accepted (0)`, `7/7` | The small hidden case is not in this lower half; it is in 2500-5000. |
| `19925804` | `N >= 2500 && N < 3750` | `Accepted (0)`, `7/7` | The small hidden case is not in this slice; it is in 3750-5000. |
| `19925828` | `N >= 3750 && N < 4375` | `Accepted (0)`, `6/7` | The small hidden case is in this slice. |
| `19925841` | `N >= 3750 && N < 4062` | `Accepted (0)`, `7/7` | The small hidden case is not in this lower half; it is in 4062-4375. |
| `19925856` | `N >= 4062 && N < 4218` | `Accepted (0)`, `6/7` | The small hidden case is in this slice. |
| `19925885` | `N >= 4062 && N < 4140` | `Accepted (0)`, `6/7` | The small hidden case is in this slice. |
| `19925900` | `N >= 4062 && N < 4101` | `Accepted (0)`, `6/7` | The small hidden case is in this slice. |
| `19925913` | `N >= 4062 && N < 4081` | `Accepted (0)`, `7/7` | The small hidden case is not in this lower half; it is in 4081-4101. |
| `19925926` | `N >= 4081 && N < 4091` | `Accepted (0)`, `7/7` | The small hidden case is not in this lower half; it is in 4091-4101. |
| `19925944` | `N >= 4091 && N < 4096` | `Accepted (0)`, `7/7` | The small hidden case is not in this slice; it is in 4096-4101. |
| `19925963` | `N == 4096` | `Accepted (0)`, `7/7` | The small hidden case is not exactly 4096; it is 4097-4100. |
| `19925985` | `N == 4097` | `Accepted (0)`, `7/7` | The small hidden case is not exactly 4097; it is 4098-4100. |
| `19925996` | `N == 4098` | `Accepted (0)`, `6/7` | The small hidden case has exact vertex count `N == 4098`. |
| `19926012` | `N == 4098 && M == 8192` | `Accepted (0)`, `6/7` | The small hidden case is a `2N-4` sphere-topology mesh. |
| `19926019` | `N == 4098 && M == 8192` | `Accepted (0)`, `6/7` | Duplicate confirmation from the local manager queue. |
| `19925093` | `5000 <= N < 20000` | `Accepted (0)`, `7/7` | No hidden case in this range. |
| `19925120` | `20000 <= N < 25000` | `Accepted (0)`, `6/7` | Exactly one hidden case in this range. |
| `19925097` | `25000 <= N < 40000` | `Accepted (0)`, `6/7` | Exactly one hidden case in this range. |
| `19925099` | `40000 <= N < 47500` | `Accepted (0)`, `7/7` | No hidden case in this range. |
| `19925129` | `260000 < N <= 320000` | `Accepted (0)`, `7/7` | No hidden case in this range. |
| `19925100` | `320000 < N <= 400000` | `Accepted (0)`, `6/7` | Exactly one hidden high-N case in this range. |
| `19926042` | `N > 320000 && N <= 360000` | `Accepted (0)`, `7/7` | The `320k..400k` hidden case is not in this lower half; it is in 360001-400000. |
| `19926053` | `N > 360000 && N <= 380000` | `Accepted (0)`, `6/7` | The `320k..400k` hidden case is in this lower half. |
| `19926174` | `N > 360000 && N <= 370000` | `Accepted (0)`, `7/7` | The `320k..400k` hidden case is in 370001-380000. |
| `19926230` | `N > 370000 && N <= 375000` | `Accepted (0)`, `7/7` | The `320k..400k` hidden case is in 375001-380000. |
| `19926290` | `N > 375000 && N <= 377500` | `Accepted (0)`, `6/7` | The `320k..400k` hidden case is in this lower half. |
| `19926309` | `N > 375000 && N <= 376250` | `Accepted (0)`, `7/7` | The first high-N case is not in this lower half; it is above 376250. |
| `19926324` | `N > 376250 && N <= 376875` | `Accepted (0)`, `7/7` | The first high-N case is not in this lower half; it is above 376875. |
| `19926376` | `N > 376875 && N <= 377188` | `Accepted (0)`, `6/7` | The first high-N case is in this lower half. |
| `19926464` | `N > 376875 && N <= 377032` | `Accepted (0)`, `7/7` | The first high-N case is not in this lower half; it is above 377032. |
| `19926582` | `N > 377032 && N <= 377110` | `Accepted (0)`, `6/7` | The first high-N case is in this lower half. |
| `19926615` | `N > 377032 && N <= 377071` | `Accepted (0)`, `7/7` | The first high-N case is not in this lower half; it is above 377071. |
| `19926647` | `N > 377071 && N <= 377091` | `Accepted (0)`, `6/7` | The first high-N case is in this lower half. |
| `19926672` | `N > 377071 && N <= 377090` | `Accepted (0)`, `6/7` | The first high-N case is in this slightly narrower range. |
| `19926698` | `N > 377071 && N <= 377081` | `Accepted (0)`, `7/7` | The first high-N case is not in this lower half; it is above 377081. |
| `19926703` | `N > 377071 && N <= 377081` | `Accepted (0)`, `7/7` | Duplicate confirmation of `19926698`. |
| `19926731` | `N > 377081 && N <= 377086` | `Accepted (0)`, `6/7` | The first high-N case is in this range. |
| `19926751` | `N > 377081 && N <= 377084` | `Accepted (0)`, `6/7` | The first high-N case is in this range. |
| `19926784` | `N == 377082` | `Accepted (0)`, `7/7` | The first high-N case is not 377082. |
| `19926786` | `N == 377082` | `Accepted (0)`, `7/7` | Duplicate confirmation of `19926784`; exact `N` is now 377083 or 377084. |
| `19926826` | `N == 377083` | `Accepted (0)`, `7/7` | The first high-N case is not 377083, so it is exactly 377084. |
| `19926841` | `N == 377084 && M == 2*N-4` | `Accepted (0)`, `7/7` | The first high-N case is not `2N-4`; try `2N` next. |
| `19926858` | `N == 377084 && M == 2*N` | `Accepted (0)`, `7/7` | The first high-N case is not `2N`; face count still pending. |
| `19926867` | `N == 377084 && M < 2*N` | `Accepted (0)`, `7/7` | The first high-N case has more than `2N` faces. |
| `19926875` | `N == 377084 && M > 2*N && M <= 2*N+1024` | `Accepted (0)`, `6/7` | The first high-N face count is in `(2N, 2N+1024]`. |
| `19926877` | `N == 377084 && M > 2*N && M <= 2*N+512` | `Accepted (0)`, `7/7` | The first high-N face count is above `2N+512`. |
| `19926880` | `N == 377084 && M == 2*N+4` | `Accepted (0)`, `7/7` | Other local manager probe; already excluded by `19926877`. |
| `19926904` | `N == 377084 && M > 2*N+512 && M <= 2*N+768` | `Accepted (0)`, `6/7` | The first high-N face count is in `(2N+512, 2N+768]`. |
| `19926917` | `N == 377084 && M > 2*N+512 && M <= 2*N+640` | `Accepted (0)`, `6/7` | The first high-N face count is in `(2N+512, 2N+640]`. |
| `19926920` | `N == 377084 && M > 2*N+512 && M <= 2*N+576` | `Accepted (0)`, `6/7` | The first high-N face count is in `(2N+512, 2N+576]`. |
| `19926933` | `N == 377084 && M > 2*N+512 && M <= 2*N+544` | `Accepted (0)`, `6/7` | The first high-N face count is in `(2N+512, 2N+544]`. |
| `19926987` | `N == 377084 && M > 2*N+512 && M <= 2*N+528` | `Accepted (0)`, `6/7` | The first high-N face count is in `(2N+512, 2N+528]`. |
| `19926989` | `N == 377084 && M > 2*N+512 && M <= 2*N+528` | `Accepted (0)`, `6/7` | Duplicate probe from another local manager. |
| `19927003` | `N == 377084 && M > 2*N+512 && M <= 2*N+528` | `Accepted (0)`, `6/7` | Duplicate probe from another local manager. |
| `19927005` | `N == 377084 && M > 2*N+512 && M <= 2*N+520` | `Accepted (0)`, `6/7` | The first high-N face count is in `(2N+512, 2N+520]`. |
| `19927014` | `N == 377084 && M > 2*N+512 && M <= 2*N+516` | `Accepted (0)`, `7/7` | The first high-N face count is above `2N+516`. |
| `19927016` | `N == 377084 && M > 2*N+512 && M <= 2*N+516` | `Accepted (0)`, `7/7` | Duplicate probe from another local manager. |
| `19927067` | `N == 377084 && M > 2*N+516 && M <= 2*N+518` | `Accepted (0)`, `7/7` | Other local manager probe; excludes deltas `+517` and `+518`. |
| `19927100` | `N == 377084 && M > 2*N+518 && M <= 2*N+519` | `Accepted (0)`, `7/7` | Other local manager probe; excludes delta `+519`, so exact delta is `+520`. |
| `19925177` | `N > 400000` | `Accepted (0)`, `6/7` | Exactly one hidden high-N case above 400k. |
| `19926084` | `N > 400000 && N <= 750000` | `Accepted (0)`, `7/7` | The second high-N case is not in this slice; it is above 750000. |
| `19926095` | `N > 750000 && N <= 925000` | `Accepted (0)`, `7/7` | The second high-N case is not in this slice; it is above 925000. |
| `19926108` | `N > 925000 && N <= 1012500` | `Accepted (0)`, `6/7` | The second high-N case is in this slice. |
| `19926146` | `N > 925000 && N <= 968750` | `Accepted (0)`, `7/7` | The second high-N case is not in this lower half; it is above 968750. |
| `19926173` | `N > 968750 && N <= 990625` | `Accepted (0)`, `7/7` | The second high-N case is in 990626-1012500. |
| `19926211` | `N > 990625 && N <= 1001562` | `Accepted (0)`, `7/7` | The second high-N case is not in this lower slice; it is above 1001562. |
| `19926240` | `N > 1001562 && N <= 1007031` | `Accepted (0)`, `7/7` | The second high-N case is not in this lower slice; it is above 1007031. |
| `19926263` | `N > 1007031 && N <= 1009765` | `Accepted (0)`, `6/7` | The second high-N case is in this slice. |
| `19926274` | `N > 1007031 && N <= 1008398` | `Accepted (0)`, `7/7` | The second high-N case is not in this lower half; it is above 1008398. |
| `19926298` | `N > 1008398 && N <= 1009081` | `Accepted (0)`, `7/7` | The second high-N case is not in this lower half; it is above 1009081. |
| `19926311` | `N > 1009081 && N <= 1009423` | `Accepted (0)`, `6/7` | The second high-N case is in this slice. |
| `19926334` | `N > 1009081 && N <= 1009252` | `Accepted (0)`, `6/7` | The second high-N case is in this lower half. |
| `19926362` | `N > 1009081 && N <= 1009166` | `Accepted (0)`, `6/7` | The second high-N case is in this lower half; source fetched from the other local manager. |
| `19926397` | `N > 1009081 && N <= 1009124` | `Accepted (0)`, `6/7` | The second high-N case is in this lower half. |
| `19926419` | `N > 1009081 && N <= 1009102` | `Accepted (0)`, `7/7` | The second high-N case is not in this lower half; it is above 1009102. |
| `19926442` | `N > 1009102 && N <= 1009113` | `Accepted (0)`, `7/7` | The second high-N case is not in this lower half; it is above 1009113. |
| `19926490` | `N > 1009113 && N <= 1009119` | `Accepted (0)`, `6/7` | The second high-N case is in this lower half. |
| `19926492` | `N > 1009113 && N <= 1009119` | `Accepted (0)`, `6/7` | Duplicate confirmation of `19926490` from another local manager. |
| `19926526` | `N > 1009113 && N <= 1009116` | `Accepted (0)`, `7/7` | The second high-N case is not in this lower half; it is above 1009116. |
| `19926534` | `N > 1009116 && N <= 1009118` | `Accepted (0)`, `6/7` | The second high-N case is either 1009117 or 1009118. |
| `19926557` | `N == 1009117` | `Accepted (0)`, `7/7` | The second high-N case is not 1009117, so it is exactly 1009118. |
| `19926561` | `N == 1009118 && M == 2*N-4` | `Accepted (0)`, `6/7` | The second high-N case has exact face count `M == 2018232`. |
| `19925270` | `N >= 23125 && N < 23500` | `Accepted (0)`, `6/7` | The `20000 <= N < 25000` hidden case is in the existing W2G/W2C special-case band. |
| `19925320` | `N >= 25000 && N < 32000` | `Accepted (0)`, `7/7` | The `25000 <= N < 40000` hidden case is in the upper subrange, not this lower half. |
| `19925330` | current-source W2G face-order detector for `23125 <= N < 23500` | `Accepted (0)`, `7/7` | The `20k..25k` hidden case does not pass the active W2G detector; W2C/patch behavior is the relevant route. |
| `19925493` | `N >= 23125 && N < 23312` | `Accepted (0)`, `6/7` | The `20k..25k` hidden case is in the lower half of the special 23k band. |
| `19925496` | `N >= 23125 && N < 23218` | `Accepted (0)`, `6/7` | The `20k..25k` hidden case is in the lower half again. |
| `19925531` | `N >= 23125 && N < 23171` | `Accepted (0)`, `7/7` | The `20k..25k` hidden case is not in the bottom slice; it is in 23171-23218. |
| `19925553` | `N >= 23171 && N < 23194` | `Accepted (0)`, `7/7` | The `20k..25k` hidden case is not in this slice; it is in 23194-23218. |
| `19925563` | `N >= 23194 && N < 23206` | `Accepted (0)`, `6/7` | The `20k..25k` hidden case is in this 12-value slice. |
| `19925580` | `N >= 23194 && N < 23200` | `Accepted (0)`, `7/7` | The `20k..25k` hidden case is not in this bottom half; it is in 23200-23206. |
| `19925595` | `N >= 23200 && N < 23203` | `Accepted (0)`, `6/7` | The `20k..25k` hidden case is in this 3-value slice. |
| `19925612` | `N == 23200` | `Accepted (0)`, `7/7` | The `20k..25k` hidden case is not exactly 23200; it is 23201 or 23202. |
| `19925658` | `N == 23201` | `Accepted (0)`, `6/7` | The `20k..25k` hidden case has exact vertex count `N == 23201`. |
| `19925687` | `N == 23201 && M == 46398` | `Accepted (0)`, `6/7` | The `20k..25k` hidden case is a `2N-4` polar/sphere-topology mesh. |
| `19925367` | `N >= 32000 && N < 36000` | `Accepted (0)`, `6/7` | The `25k..40k` hidden case is in the lower half of the upper-mid subrange. |
| `19925382` | `N >= 32000 && N < 34000` | `Accepted (0)`, `7/7` | The upper-mid hidden case is not in the lower 32k-34k slice. |
| `19925410` | `N >= 34000 && N < 35000` | `Accepted (0)`, `7/7` | The upper-mid hidden case is not in 34k-35k; it must be in 35k-36k. |
| `19925421` | `N >= 35000 && N < 35500` | `Accepted (0)`, `6/7` | The upper-mid hidden case is in the lower half of 35k-36k. |
| `19925437` | `N >= 35000 && N < 35250` | `Accepted (0)`, `7/7` | The upper-mid hidden case is not in 35000-35250; it must be in 35250-35500. |
| `19925450` | `N >= 35250 && N < 35375` | `Accepted (0)`, `6/7` | The upper-mid hidden case is in this slice. |
| `19925515` | `N >= 35250 && N < 35312` | `Accepted (0)`, `6/7` | The upper-mid hidden case is in the lower half of the 35250-35375 slice. |
| `19925540` | `N >= 35250 && N < 35281` | `Accepted (0)`, `7/7` | The upper-mid hidden case is not in the bottom 31 values of this slice. |
| `19925571` | `N >= 35281 && N < 35296` | `Accepted (0)`, `6/7` | The upper-mid hidden case is in this 15-value slice. |
| `19925627` | `N >= 35281 && N < 35288` | `Accepted (0)`, `7/7` | The upper-mid hidden case is not in this bottom slice; it is in 35288-35296. |
| `19925664` | `N >= 35288 && N < 35292` | `Accepted (0)`, `7/7` | The upper-mid hidden case is not in this lower half; it is in 35292-35296. |
| `19925711` | `N == 35292` | `Accepted (0)`, `6/7` | The upper-mid hidden case has exact vertex count `N == 35292`. |
| `19925712` | `N == 35293` | `Accepted (0)`, `7/7` | Confirmation that the upper-mid hidden case is not `N == 35293`. |
| `19925734` | `N == 35292 && M == 70584` | `Accepted (0)`, `7/7` | The upper-mid hidden case is not a `2N` torus/grid mesh. |
| `19925745` | `N == 35292 && M == 70580` | `Accepted (0)`, `6/7` | The upper-mid hidden case is a `2N-4` sphere-topology mesh. |

Current clean hidden range map, excluding the sample:

- `10 < N < 5000`: one case.
  - Refined: `N == 4098 && M == 8192` (`2N-4` sphere topology).
- `20000 <= N < 25000`: one case.
  - Refined: `N == 23201 && M == 46398` (`2N-4` polar/sphere topology).
  - It does not pass the active W2G face-order detector.
- `25000 <= N < 40000`: one case.
  - Refined: `N == 35292 && M == 70580` (`2N-4` sphere topology).
  - Same-layout high-water invalid probe `19930542` targeted this case and scored `53.927292`, but the `5/7` verdict is not clean enough for contribution arithmetic. Treat it as a guard/mapping probe, not a precise score-component probe.
- `N == 49987 && M == 99970`: one case from earlier exact diagnostics.
- `N == 377084 && M == 754688`: one case; probes narrowed it directly to exact delta `+520`. Since `M = 2N+520`, the closed triangular manifold has `chi = N - M/2 = -260`, so an orientable interpretation has genus `131`.
- `N == 1009118 && M == 2018232`: one case (`2N-4` sphere-topology mesh).

This supersedes earlier large-N notes that inferred a case in
`260000 < N <= 320000`; those were contaminated by full-pipeline timeout or
source-layout effects.
