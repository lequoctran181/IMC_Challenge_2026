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
| `submission_1327_53.93_5.cpp` | `19920478` | `53.926625`, 5/7 | external concurrent source; same low bucket |
| `submission_1328_53.93_5.cpp` | `19920551` | `53.926625`, 5/7 | third B16 76 -> 88 plus last high-N B16 72 -> 88; unsafe |
| `submission_1329_53.93_5.cpp` | `19920693` | `53.926625`, 5/7 | external concurrent source; 131026-byte patch, low |
| `submission_1330_43.09_4.cpp` | `19920733` | `43.093759`, 4/7 | first high-N broad count 120 -> 160 alone; very unsafe |
| `submission_1331_53.93_5.cpp` | `19920738` | `53.926625`, 5/7 | second high-N broad count 72 -> 96 alone; unsafe |
| `submission_1332_53.96_5.cpp` | `19920774` | `53.958233`, 5/7 | high-N q relaxed further than 1319; unsafe |
| `submission_1333_53.93_5.cpp` | `19920778` | `53.926625`, 5/7 | external concurrent source, hash `90beacb7`; same low bucket |
| `submission_1334_53.93_5.cpp` | `19920806` | `53.927292`, 5/7 | exact-base `S3B16::T` count 10 -> 12; unsafe |
| `submission_1335_56.93_5.cpp` | `19920833` | `56.926922`, 5/7 | external/pro minified 102 KB source; unsafe |
| `submission_1336_57.28_5.cpp` | `19920870` | `57.279876`, 5/7 | external/pro branch; unsafe |
| `submission_1337_30.67_3.cpp` | `19920842` | `30.666004`, 3/7 | Pro Extended standalone proxy-router; very unsafe |
| `submission_1338_0.00_6.cpp` | `19920866` | `0.000000`, 6/7 | Pro Extended 49843..50625 standalone branch; accepted but below visual threshold |

## Current Lessons

- Keep `submission_1181_81.95_7.cpp` / `submission_1255_min1181_PENDING.cpp` as the only trusted base.
- Small high-N B16 q/count variants can tie exact-best but have not improved it.
- Do not use the current `macro_pack` output as a behavior-preserving shrink; it changes hidden behavior despite sample passing.
- B92 structural branches rebased onto exact-best still false-positive or timeout.
- Adding extra early `VIMP::run()` is not a harmless no-op and should be blacklisted.
- Increasing both broad high-N B16 counts is unsafe.
- Increasing multiple B16 passes in one source is unsafe unless each pass was proven alone first; Kattis 19920551 repeated the 53.926625 failure mode.
- The first high-N broad B16 count is especially brittle: `120 -> 160` alone fell to 43.093759/4.
- The second high-N broad B16 has a narrow safe window: `72 -> 88` tied exact best, but `72 -> 96` fell to 53.926625/5.
- High-N B16 q relaxation also has a cliff: the first tiny relaxation tied exact best, but the next relaxation fell to 53.958233/5.
- External source hash `90beacb7` also lands in the 53.926625/5 bucket; do not reuse it as a trusted base.
- `S3B16::T` count must stay conservative: changing 10 -> 12 fell to 53.927292/5.
- Pro Extended standalone/weak-fallback rewrites can pass sample but collapse badly on hidden tests; prioritize exact-best-preserving patches.
- A standalone hidden-band branch can pass 6/7 validity yet score 0; visual threshold risk dominates if exact-best fallback is absent.
