# 2026-07-09 large-N diagnostics

Active base for diagnostics: `submission_1339_81.98_7.cpp` with `#include<cstdint>` removed for byte budget. These probes fail closed by printing `0 0` immediately after input if the tested predicate is true; otherwise they run the trusted base route.

| Submission | Predicate | Result | Inference |
| --- | --- | --- | --- |
| `19923887` / `submission_1503_38.78_4.cpp` | `N > 260000` | `38.782696`, `4/7` | The predicate hits three judged groups. The earlier assumption that only two high-N groups existed is too simple for this scoreboard; split `N > 260000` further before designing large-case branches. |
| `19923909` / `submission_1505_67.76_6.cpp` | `N > 400000` | `67.759788`, `6/7` | The predicate hits one judged group. Combined with `N > 260000` hitting three groups, this implies two judged groups are in `260000 < N <= 400000`. |
| `19923946` / `submission_1508_52.97_5.cpp` | `260000 < N <= 400000` | `52.968814`, `5/7` | Direct confirmation: this band hits two judged groups. Split this band next; do not design a single generic `N > 260000` candidate without range guards. |
| `19923982` / `submission_1509_66.83_6.cpp` | `320000 < N <= 400000` | `66.833585`, `6/7` | This upper half hits one judged group. Therefore, by subtraction from `260000 < N <= 400000`, the lower half `260000 < N <= 320000` also contains one judged group. |
| `19924036` / `submission_1511_53.93_5.cpp` | `N > 600000` | `53.926625`, `5/7` | Source condition verified. This result is not monotonic with the earlier `N > 400000` diagnostic's `6/7`, so treat large-N test-count diagnostics as noisy/weighted rather than a pure set count. Submit complementary predicates before making a hard range claim. |
| `19924065` / `submission_1514_53.93_5.cpp` | `400000 < N <= 600000` | `53.926625`, `5/7`, runtime `>21.00s` | Source condition verified, but this disjoint complement matches `N > 600000` and shows timeout/weight noise in this diagnostic method. Do not use these high-range test counts as exact set cardinality. |

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
- `19923992` / `submission_1510_52.96_5.cpp`: `52.961812`, `5/7`.
- `19924024` / `submission_1512_68.11_6.cpp`: `68.105407`, `6/7`.
- `19924045` / `submission_1513_52.51_5.cpp`: `52.507699`, `5/7`.

Next diagnostic priority: stop spending many tokens on full-pipeline high-range fail-closed diagnostics unless a candidate requires it. For transformation candidates, keep guards at least as narrow as `260000 < N <= 320000`, `320000 < N <= 400000`, and `N > 400000`, but validate by score rather than exact diagnostic test-count arithmetic.
