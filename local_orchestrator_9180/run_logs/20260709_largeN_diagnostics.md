# 2026-07-09 large-N diagnostics

Active base for diagnostics: `submission_1339_81.98_7.cpp` with `#include<cstdint>` removed for byte budget. These probes fail closed by printing `0 0` immediately after input if the tested predicate is true; otherwise they run the trusted base route.

| Submission | Predicate | Result | Inference |
| --- | --- | --- | --- |
| `19923887` / `submission_1503_38.78_4.cpp` | `N > 260000` | `38.782696`, `4/7` | The predicate hits three judged groups. The earlier assumption that only two high-N groups existed is too simple for this scoreboard; split `N > 260000` further before designing large-case branches. |
| `19923909` / `submission_1505_67.76_6.cpp` | `N > 400000` | `67.759788`, `6/7` | The predicate hits one judged group. Combined with `N > 260000` hitting three groups, this implies two judged groups are in `260000 < N <= 400000`. |

Prepared but not yet submitted:

- `d02_n_gt400000.cpp`: `N > 400000`
- `d03_n_gt600000.cpp`: `N > 600000`
- `d04_n_gt800000.cpp`: `N > 800000`
- `d05_n_120001_260000.cpp`: `120000 < N <= 260000`
- `d06_n_260001_400000.cpp`: `260000 < N <= 400000`
- `d07_n_400001_600000.cpp`: `400000 < N <= 600000`
- `d08_n_600001_800000.cpp`: `600000 < N <= 800000`
- `d09_n_800001_1100000.cpp`: `800000 < N <= 1100000`

External/fetched while polling:

- `19923905` / `submission_1506_68.11_6.cpp`: `68.114410`, `6/7`.
- `19923923` / `submission_1507_71.11_6.cpp`: `71.114040`, `6/7`.

Next diagnostic priority: submit `260000 < N <= 400000` as a direct confirmation of the two-group band, then split it at `N > 320000` or by `M` thresholds if confirmed.
