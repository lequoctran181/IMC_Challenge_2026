# 2026-07-09 range-guarded B16 candidates

Active base: `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact score `81.978181`, `7/7`.

The large-N diagnostics suggest using range guards instead of global high-N changes, but the high-range fail-closed tests also showed runtime/test-count noise. Candidate validation must be by final score and accepted test count, not by diagnostic arithmetic alone.

## Submitted candidates

| Submission | Candidate | Change | Result | Conclusion |
| --- | --- | --- | --- | --- |
| `19924092` / `submission_1515_53.93_5.cpp` | `c01_first121_r1` | In the first broad `B16::R`, use count `121` only for `260000 < N <= 320000`; otherwise keep count `120`. | `53.927292`, `5/7`, runtime `>21.00s` | Unsafe. The old global `120 -> 121` cliff is not rescued by this `r1` guard. Do not submit equivalent first-count increases for this range without a stronger rollback/validator. |

## External/fetched while polling

- `19924086` / `submission_1516_52.96_5.cpp`: `52.961478`, `5/7`.
- `19924100` / `submission_1517_52.97_5.cpp`: `52.969480`, `5/7`.

## Prepared local candidates

Compiled and sample-checked under `/tmp/imc_9180_main/queue16_20260709_range_candidates/`:

- `c01_first121_r2.cpp`
- `c01_first121_r3.cpp`
- `c04_last89_r1.cpp`
- `c04_last89_r2.cpp`
- `c04_last89_r3.cpp`
- `c07_first119_r1.cpp`
- `c07_first119_r2.cpp`
- `c09_last87_r1.cpp`
- `c09_last87_r2.cpp`
- `c11_vp96_mid260_400.cpp`
- `c12_vp128_low260_320.cpp`
- `c13_vp96_hi400.cpp`
- `c14_first_q945_r1.cpp`
- `c15_last_q950_r1.cpp`
- `c16_first_q939_r1.cpp`

Next priority: prefer proxy-resolution or threshold variants over more count increases, because count increases repeatedly land in the `53.xx`/5 bucket.
