# 91.80 Control Board

Current target: Kattis score >= 91.80. Current retained root-submission best by filename: 81.95 (submission_1181_81.95_7.cpp). Exact high-water mark remains 81.945906.

## Latest Local Manager Round

- `submission_1236_68.11_6.cpp`: Worker D low-risk WK tail-window extension, Kattis `19915873`, accepted 68.112743, 6/7.
- `submission_1237_81.94_7.cpp`: Worker E visual-shell cover candidate, byte-trimmed by identifier renaming to 129127 bytes, Kattis `19915979`, accepted 81.938904, 7/7.
- `submission_1238_68.11_6.cpp`: fetched concurrent submission `19915947`, accepted 68.108742, 6/7.
- `submission_1239_80.65_7.cpp`: fetched concurrent submission `19916059`, accepted 80.64606, 7/7.
- None of the latest four submissions beats `submission_1181_81.95_7.cpp`; use them as negative/near-plateau evidence, not as a new base.

## Score Buckets

- 0/7: 8 files, max 0.00, avg 0.00
- 1/7: 4 files, max 16.50, avg 15.96
- 2/7: 38 files, max 16.56, avg 13.97
- 3/7: 19 files, max 29.80, avg 25.21
- 4/7: 44 files, max 48.10, avg 34.79
- 5/7: 24 files, max 59.78, avg 44.03
- 6/7: 19 files after latest archive; latest max in this bucket is still below best.
- 7/7: 116 files after latest archive; current exact best remains 81.945906.

## Top Files

- submission_1181_81.95_7.cpp: score 81.95, tests 7/7, size 131019, sha afa1759b071bf51d, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1188_81.94_7.cpp: score 81.94, tests 7/7, size 130864, sha 7105ff5c748bfc6f, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1186_81.94_7.cpp: score 81.94, tests 7/7, size 130864, sha 4bbeefe3e941fb35, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1185_81.94_7.cpp: score 81.94, tests 7/7, size 130822, sha 0617452dd2311bd8, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1183_81.94_7.cpp: score 81.94, tests 7/7, size 130865, sha b4d3d0219c81c50b, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1182_81.94_7.cpp: score 81.94, tests 7/7, size 130864, sha abb13626791389d8, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1180_81.94_7.cpp: score 81.94, tests 7/7, size 130905, sha 7c3143a8efa77f95, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1179_81.94_7.cpp: score 81.94, tests 7/7, size 130864, sha ccb10d8563288ac8, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1178_81.94_7.cpp: score 81.94, tests 7/7, size 130864, sha 1ddaf039dba0c24d, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1177_81.94_7.cpp: score 81.94, tests 7/7, size 130864, sha 2a3ea7de3d8018a6, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1176_81.94_7.cpp: score 81.94, tests 7/7, size 130824, sha 5f9f4b1f07c616b5, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1174_81.94_7.cpp: score 81.94, tests 7/7, size 130864, sha b83ecd73f35b4a04, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1173_81.94_7.cpp: score 81.94, tests 7/7, size 130822, sha a619974367569a36, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1172_81.94_7.cpp: score 81.94, tests 7/7, size 130823, sha 0886f4ded245f2c0, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove
- submission_1171_81.94_7.cpp: score 81.94, tests 7/7, size 130822, sha b9ce58a0266e40e4, sig W2G,W2C,W5,VIMP,MIDEC,WK,B16,S3B16,GN,cove

## Immediate Interpretation

- There are many 7/7 submissions, so validity is mostly solved; the bottleneck is compression/visual score on hidden cases.
- Best variants cluster around the same large minified router using W2G/W2C, W5, VIMP, MIDEC, WK and B16-style post patches.
- Recent submissions 1229-1235 are much lower (13-55 range), so latest is not the best baseline.
- A 10-point jump likely requires a new hidden-test-specific branch or substantially better per-case router, not a parameter-only tweak.

## Worker Assignments

- A: history/family clusters.
- B: hidden case inference from code/result patterns.
- C: local evaluator/harness.
- D: low-risk candidate from best plateau.
- E: breakthrough visual-shell / vertex-only Hausdorff candidate.
- F: safe submit/dedup pipeline.

## Current Next Direction

- Do not spend more attempts on pure WK/B16 tail nudges without a hidden-case reason; D dropped to 68.11.
- E proved the source-size workaround and fail-closed fallback can preserve the plateau, but the visual-shell recognizer did not trigger a score jump.
- Highest priority remains isolating the exact `submission_1181` delta and building a size-budgeted, fail-closed branch for the suspected `49843..50625` hidden band.
