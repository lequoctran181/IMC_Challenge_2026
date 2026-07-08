# Round 1320-1326 Notes

Target remains `91.80+`; exact best remains `81.945906`.

## Submitted Results

| root file | kattis id | result | note |
|---|---:|---:|---|
| `submission_1319_81.95_7.cpp` | `19920035` | `81.945906`, 7/7 | exact-best high-N q relaxed slightly; tied best |
| `submission_1320_81.95_7.cpp` | `19920036` | `81.945906`, 7/7 | exact-best second broad count 72 -> 88; tied best |
| `submission_1321_57.28_5.cpp` | `19920108` | `57.279876`, 5/7 | exact-base B92 periodic branch; unsafe |
| `submission_1322_67.76_6.cpp` | `19920109` | `67.759788`, 6/7 | exact-base B92 combo-all branch; unsafe/timeout |
| `submission_1323_53.96_5.cpp` | `19920177` | `53.958233`, 5/7 | external concurrent source; low |
| `submission_1324_57.28_5.cpp` | `19920248` | `57.279876`, 5/7 | macro pack-only exact-base; behavior not preserved |
| `submission_1325_56.31_5.cpp` | `19920311` | `56.305810`, 5/7 | duplicated early `VIMP::run()`; unsafe |
| `submission_1326_53.93_5.cpp` | `19920383` | `53.926625`, 5/7 | broad high-N counts 120/72 -> 160/96; unsafe |

## Current Lessons

- Keep `submission_1181_81.95_7.cpp` / `submission_1255_min1181_PENDING.cpp` as the only trusted base.
- Small high-N B16 q/count variants can tie exact-best but have not improved it.
- Do not use the current `macro_pack` output as a behavior-preserving shrink; it changes hidden behavior despite sample passing.
- B92 structural branches rebased onto exact-best still false-positive or timeout.
- Adding extra early `VIMP::run()` is not a harmless no-op and should be blacklisted.
- Increasing both broad high-N B16 counts is unsafe.
