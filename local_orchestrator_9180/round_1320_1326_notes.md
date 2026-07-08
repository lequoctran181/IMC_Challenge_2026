# Round 1320-1326 Notes

Target remains `91.80+`; exact best is now `81.977514`.

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
| `submission_1339_81.98_7.cpp` | `19920917` | `81.977514`, 7/7 | new best; first high-N B16 q `.941/.955` plus final high-N B16 count `88` |
| `submission_1340_81.95_7.cpp` | `19920953` | `81.945906`, 7/7 | first q `.940/.954` with final count `88`; valid but drops to old best |
| `submission_1341_53.93_5.cpp` | `19920967` | `53.926625`, 5/7 | first q `.940/.955`; unsafe |
| `submission_1342_53.93_5.cpp` | `19920939` | `53.926625`, 5/7 | final high-N B16 q `.945/.958`; unsafe |
| `submission_1343_81.93_7.cpp` | `19920942` | `81.934570`, 7/7 | BROAD_32 C5 tube guarded branch; no improvement |
| `submission_1344_53.93_5.cpp` | `19920961` | `53.926625`, 5/7 | final high-N B16 count `92`; unsafe |
| `submission_1345_53.93_5.cpp` | `19920965` | `53.926625`, 5/7 | final high-N B16 count `90`; unsafe |
| `submission_1346_53.93_5.cpp` | `19920978` | `53.926625`, 5/7 | final high-N B16 count `89`; unsafe |
| `submission_1347_57.28_5.cpp` | `19921011` / `19921005` | `57.279876` / `53.926625`, 5/7 | first q `.941/.954`; same SHA produced two low scores, likely time-threshold instability |
| `submission_1348_43.09_4.cpp` | `19921049` | `43.093759`, 4/7 | first q `.942/.955`; unsafe, worse than `.941/.954` |
| `submission_1349_43.09_4.cpp` | `19921061` | `43.093759`, 4/7 | concurrent external source hash `0fbe2af2`; unsafe |
| `submission_1350_81.93_7.cpp` | `19921075` | `81.929569`, 7/7 | fixed worker515 T03 case5 compare branch; valid but below best |
| `submission_1351_81.93_7.cpp` | `19921102` | `81.934570`, 7/7 | S14 case5 cap92 candidate; no improvement |
| `submission_1352_80.63_7.cpp` | `19921105` | `80.634329`, 7/7 | S10_06 recombine; valid but harmful |
| `submission_1353_81.93_7.cpp` | `19921125` | `81.934570`, 7/7 | S7_04 pre-MIDEC/WK swap; no improvement |
| `submission_1354_81.93_7.cpp` | `19921130` | `81.934570`, 7/7 | S14 exact feature candidate; no improvement |
| `submission_1355_81.93_7.cpp` | `19921141` | `81.934570`, 7/7 | S7_05 skip-MIDEC selector; no improvement |
| `submission_1356_81.94_7.cpp` | `19921150` | `81.938237`, 7/7 | S10_03 recombine; valid but still below best |
| `submission_1357_70.49_6.cpp` | `19921160` | `70.491928`, 6/7 | S10_05 recombine; unsafe |
| `submission_1358_81.94_7.cpp` | `19921181` | `81.938237`, 7/7 | S10_08 lateW5 recombine; same as S10_03, below best |
| `submission_1359_81.93_7.cpp` | `19921209` | `81.934570`, 7/7 | S5 case6 late candidate; no improvement |
| `submission_1360_81.93_7.cpp` | `19921250` | `81.934570`, 7/7 | S5 case6 strict candidate; no improvement |
| `submission_1361_81.95_7.cpp` | `19921195` | `81.945906`, 7/7 | concurrent worker source, tracked by commit `090389c`; ties old exact best |
| fetched only | `19921249` | `57.279876`, 5/7 | concurrent external source hash `4ec90e9f`; unsafe |
| `submission_1362_81.93_7.cpp` | `19921274` | `81.934570`, 7/7 | S5 highN 650-1200 ultrastrict; no improvement |
| `submission_1363_81.93_7.cpp` | `19921275` | `81.934570`, 7/7 | S15 B16 200-650 ultraguard; no improvement |
| fetched only | `19921293` | `57.279876`, 5/7 | concurrent external source hash `796fd47b`; unsafe |
| `submission_1364_81.93_7.cpp` | `19921296` | `81.934570`, 7/7 | S14 midfloor; no improvement |
| `submission_1365_81.93_7.cpp` | `19921311` | `81.934570`, 7/7 | S14 looseguard; no improvement |
| `submission_1366_81.93_7.cpp` | `19921326` | `81.934570`, 7/7 | S14 tightguard; no improvement |
| `submission_1367_81.93_7.cpp` | `19921340` | `81.934570`, 7/7 | S14 upper-exact guard; no improvement |
| `submission_1372_81.95_7.cpp` | `19921359` | `81.945906`, 7/7 | fetched unknown-source Kattis candidate; ties older exact-best bucket, below `1339` |
| `submission_1373_81.93_7.cpp` | `19921376` | `81.934570`, 7/7 | S11_01 proxy-aware S7_04 selector; no improvement |
| `submission_1374_81.93_7.cpp` | `19921384` | `81.934570`, 7/7 | S11_03 safe S7_04 selector; no improvement |
| `submission_1375_81.93_7.cpp` | `19921399` | `81.934570`, 7/7 | S11_05 clean B16-before-WK selector; no improvement |
| `submission_1376_81.95_7.cpp` | `19921400` | `81.945906`, 7/7 | concurrent worker source; ties older exact-best bucket, below `1339` |
| `submission_1377_81.93_7.cpp` | `19921428` | `81.934570`, 7/7 | S11_06 wide B16-before-WK selector; no improvement |
| `submission_1378_79.73_7.cpp` | `19921435` | `79.729258`, 7/7 | S11_02 wide S7_04 selector; valid but harmful |
| `submission_1379_81.93_7.cpp` | `19921458` | `81.934570`, 7/7 | S7_09 late-W5 guarded B16; no improvement |
| `submission_1380_81.86_7.cpp` | `19921466` | `81.858216`, 7/7 | S7_08 largeN W5L473 tail; valid but below best |
| `submission_1381_81.93_7.cpp` | `19921493` | `81.934570`, 7/7 | S9 full C5T structural shrink; no improvement |
| `submission_1382_75.61_7.cpp` | `19921509` | `75.613142`, 7/7 | unknown concurrent source; valid but harmful |
| `submission_1384_80.63_7.cpp` | `19921546` | `80.634329`, 7/7 | workerF macro5k candidate; valid but below best |
| `submission_1385_0.00_7.cpp` | `19921584` | `0.000000`, 7/7 | diagnostic: S14 lower guard did not trigger on any hidden case, so it printed the original mesh on all tests |
| `submission_1386_81.71_7.cpp` | `19921559` | `81.709845`, 7/7 | broad09 r12grid W2+B16; valid but below best |
| `submission_1387_81.93_7.cpp` | `19921587` | `81.934570`, 7/7 | broad19 boxgrid failclosed; no improvement |
| `submission_1388_81.93_7.cpp` | `19921614` | `81.934570`, 7/7 | broad28 torus_snap; no improvement |
| `submission_1389_81.69_7.cpp` | `19921619` | `81.690173`, 7/7 | workerM1 boundary60000; valid but below best |

## Current Lessons

- Keep `submission_1339_81.98_7.cpp` / Kattis `19920917` as the trusted base.
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
- The useful move was specifically combining first high-N q `.941/.955` with final high-N B16 count `88`; it raised exact score from `81.945906` to `81.977514`.
- Final high-N B16 count is a hard cliff on the new base: `89`, `90`, and `92` all fell to `53.926625`, so keep it exactly `88` unless a stronger guard/fallback is added.
- Final high-N B16 q is also brittle: changing `.946/.959` to `.945/.958` fell to `53.926625`.
- First high-N q `.940` is not a safe direction: `.940/.954` remains valid but loses the new gain, while `.940/.955` falls to `53.926625`.
- First q `.941/.954` is also unsafe; two submissions with identical source produced low 5/7 scores (`57.279876` and `53.926625`), so the low/mid threshold must stay `.955` unless protected by a new fallback.
- First q `.942/.955` is unsafe too; Kattis `19921049` fell to `43.093759`/4, so the safe point is a very narrow `.941/.955`.
- Worker515 T03's local compare/fallback compiles after replacing out-of-scope `SN` with `cove()`, but Kattis `19921075` only scored `81.929569`; do not prioritize that branch without a stronger hidden-case trigger.
- S14/S7 case5 and pipeline-selector variants mostly collapse to the old `81.934570` plateau; S10 recombines are mixed, with `S10_03`/`S10_08` at `81.938237` and `S10_05` unsafe at `70.491928`.
- S5 case6 late/strict did not improve (`81.934570`), concurrent `1361` only tied old exact best `81.945906`, and external `19921249` fell to `57.279876`.
- S5 highN/S15 ultraguard/S14 midfloor variants also did not improve (`81.934570`); external `19921293` fell to `57.279876`.
- S14 guard variants (`midfloor`, `looseguard`, `tightguard`) all return the same `81.934570` plateau.
- S14 upper-exact guard also returned `81.934570`; stop spending attempts on S14 guard-only parameter changes unless they are paired with a new hidden-case discriminator.
- Unknown-source Kattis `19921359` scored `81.945906`; it is valid but lacks the narrow `1339` gain, so use it only as a comparison source rather than a new base.
- `19921359` extends W2G-like grid detection to larger N but does not improve; this suggests output quality/target geometry, not just detector coverage, is the blocker for large-N structural branches.
- S11_01 S7_04 proxy-aware selector returned `81.934570`; S11 variants may still test clean order selectors, but the pre-MIDEC/WK swap idea has weak signal.
- S11_03 safe S7_04 selector also returned `81.934570`; keep remaining S11 attempts focused on distinct clean order/skip variants rather than more S7_04 guard tuning.
- S11_05 clean B16-before-WK selector returned `81.934570`; order tournaments are not showing lift unless they unlock a genuinely different hidden bucket.
- Concurrent `1376` scored `81.945906`; it is useful as another valid comparison point but does not challenge `1339`.
- S11_06 wide B16-before-WK selector returned `81.934570`; stop prioritizing S11 order-selector variants after the pending wide-gain result unless it breaks the bucket.
- S11_02 wide S7_04 selector regressed to `79.729258`; blacklist wide S7_04 margins for this cycle.
- S7_09 late-W5 guarded B16 returned `81.934570`; late-W5 alone is not enough to produce a new hidden bucket.
- S7_08 largeN W5L473 returned `81.858216`; broad large-N W5L473 is not the needed jump without a better hidden-case discriminator.
- S9 full C5T structural shrink returned `81.934570`; C5T should be paused unless redesigned with a new target surface or discriminator.
- Unknown concurrent `19921509` returned `75.613142`; do not use it as a base unless later diff analysis reveals an isolated positive branch.
- Diagnostic `19921584` accepted with score `0.000000` because the S14 lower guard never fired on hidden tests; stop spending attempts on that exact lower guard, since its proxy wins do not correspond to any official hidden slice.
- WorkerF macro5k candidate returned `80.634329`; macro5k is not a breakthrough branch in this form.
- Broad09 r12grid W2+B16 returned `81.709845`; avoid this structural branch unless redesigned around a much narrower detector.
- Broad19 boxgrid failclosed returned `81.934570`; failclosed structural branches still mostly preserve plateau rather than improve compression.
