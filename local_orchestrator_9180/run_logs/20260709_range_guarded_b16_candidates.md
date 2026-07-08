# 2026-07-09 range-guarded B16 candidates

Active base: `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact score `81.978181`, `7/7`.

The large-N diagnostics suggest using range guards instead of global high-N changes, but the high-range fail-closed tests also showed runtime/test-count noise. Candidate validation must be by final score and accepted test count, not by diagnostic arithmetic alone.

## Submitted candidates

| Submission | Candidate | Change | Result | Conclusion |
| --- | --- | --- | --- | --- |
| `19924092` / `submission_1515_53.93_5.cpp` | `c01_first121_r1` | In the first broad `B16::R`, use count `121` only for `260000 < N <= 320000`; otherwise keep count `120`. | `53.927292`, `5/7`, runtime `>21.00s` | Unsafe. The old global `120 -> 121` cliff is not rescued by this `r1` guard. Do not submit equivalent first-count increases for this range without a stronger rollback/validator. |
| `19924118` / `submission_1518_68.11_6.cpp` | `c14_first_q945_r1` | In the first broad `B16::R`, raise proxy threshold from `.941` to `.945` only for `260000 < N <= 320000`. | `68.113410`, `6/7` | Still unsafe. Tightening this first-call threshold in `r1` loses one judged group, so the first broad call is extremely brittle in this range. |
| `19924148` / `submission_1519_53.93_5.cpp` | `c11_vp96_mid260_400` | Use proxy resolution `96` instead of `64` for `260000 < N <= 400000` in both broad high-N B16 calls; keep `N > 400000` at `64`. | `53.927292`, `5/7`, runtime `>21.00s` | Unsafe. Raising proxy resolution in the mid-large ranges does not rescue the B16 route and lands in the same bad bucket as count increases. |

## External/fetched while polling

- `19924086` / `submission_1516_52.96_5.cpp`: `52.961478`, `5/7`.
- `19924100` / `submission_1517_52.97_5.cpp`: `52.969480`, `5/7`.
- `19924157` / `submission_1520_53.93_5.cpp`: `53.927292`, `5/7`.

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

Next priority: avoid first-call `r1` edits and mid-large proxy-resolution increases. If continuing this batch, try last-call variants only; otherwise switch back to exact case-5 strategy work.
