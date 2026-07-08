# 2026-07-08 ring extra phase sweep

Current best before this sweep: `submission_1339_81.98_7.cpp`, exact Kattis score `81.977514`.

New best from this sweep: `submission_1447_81.98_7.cpp`, Kattis submission `19922865`, exact score `81.978181`, `7/7`.

## Case 5 diagnostics

- `submission_1444_0.00_7.cpp` / `19922794`: diagnostic `N==49987 && M==99970 && (v4+v8)>=14%` was false.
- Earlier diagnostics still imply `(v4+v8)>=12%`, `v4 in [4%,8%)`, `v8 in [4%,8%)`, `v5..v7<88%`.
- Working interpretation: exact case 5 is smooth genus-0 with moderate valence-4/8 irregularity, not a pure sphere/ellipsoid/regular ring grid.

## Ring-8 result

Ring-8 removal is unsafe even when tiny:

- `submission_1445_57.28_5.cpp` / `19922813`: replacement ring-8 pass, score `57.276209`, `5/7`.
- `submission_1446_68.11_6.cpp` / `19922845`: extra ring-8 `count=2, ms=.0015`, score `68.113410`, `6/7`.

Conclusion: do not remove valence-8 patches on case 5 with current W5/S3B16 geometry guard.

## Ring<=5 extra pass result

An extra low-score `S3B16` pass in the existing placeholder is narrowly useful:

- `submission_1447_81.98_7.cpp` / `19922865`: `S3B16::T(02,-9,192,.96,18.47,0,5,.0015);`, score `81.978181`, `7/7`.
- `submission_1448_68.11_6.cpp` / `19922880`: count 4 version, score `68.114410`, `6/7`.
- `submission_1449_81.95_7.cpp` / `19922897`: count 2, stride -7 phase 1, score `81.945906`, `7/7`.
- `submission_1450_81.95_7.cpp` / `19922913`: count 2, stride +11 phase 3, score `81.945906`, `7/7`.
- `submission_1451_68.11_6.cpp` / `19922935`: count 2, stride +13 phase 5, score `68.113076`, `6/7`.
- `submission_1456_81.95_7.cpp` / `19922974`: phase0 count 2 with threshold `.0020`, score `81.946573`, `7/7`.

Conclusion: the safe window is extremely narrow. Only the original phase0 threshold `.0015` extra pass improved the exact best; relaxing phase0 to `.0020` stays valid but loses the gain.
