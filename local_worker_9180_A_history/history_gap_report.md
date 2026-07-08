# Worker A History / Gap Report

Scope: current root `submission_*.cpp` files plus git/log/fetched-source clues available in `IMC_Challenge_2026_remote`. No Kattis submission was made.

## Executive Findings

- Parsed 272 current root submissions, sequence range 967..1239.
- Current root best filename score is 81.95; exact git/Kattis metadata identifies the high-water mark as `kattis_19912924_81.945906_7.cpp`, mapped to `submission_1181_81.95_7.cpp` by the round-84 commit.
- The dominant plateau is 81.93-81.95: 86 current root files, 91 files at >=81.90.
- Accepted-count is not enough evidence: 7/7 appears in high, plateau, and failure families. The root has 114 files with 7 accepted tests, but many are far below the high score.
- The best evidenced code family is the central B16/WK/W5/VIMP/MIDEC route with `S3B16::T`; the only exact >81.945 result adds broad B16 high-N calls around that central route. Later round 90-93 plateau commits return to exact `81.93457`.
- The gap to 91.80 is not a small constant-tuning gap: exact best `81.945906` is still about 9.854 points short, and many new structural families either tie the plateau or drop into 0-56 score bands.

## Score And Test Distribution

### Accepted-count buckets

| accepted tests | count |
|---:|---:|
| 0 | 8 |
| 1 | 4 |
| 2 | 38 |
| 3 | 19 |
| 4 | 44 |
| 5 | 24 |
| 6 | 19 |
| 7 | 116 |

### Score bands

| score band | count |
|---|---:|
| 0 | 16 |
| 0-20 | 75 |
| 20-40 | 29 |
| 40-60 | 41 |
| 60-80 | 12 |
| 80-81.90 | 8 |
| 81.90-81.93 | 5 |
| 81.93-81.95 | 85 |
| 81.95+ | 1 |

### Most common rounded filename scores

| score token | count |
|---:|---:|
| 81.93 | 58 |
| 81.94 | 27 |
| 0.00 | 16 |
| 16.30 | 7 |
| 29.33 | 6 |
| 53.93 | 6 |
| 81.91 | 5 |
| 16.47 | 5 |
| 17.00 | 4 |
| 70.49 | 4 |
| 16.54 | 4 |
| 16.53 | 4 |
| 17.01 | 3 |
| 15.58 | 3 |
| 81.69 | 3 |
| 16.50 | 3 |
| 16.56 | 3 |
| 48.10 | 3 |

## Top Submissions And Exact Metadata

| root file | filename score | exact Kattis score | Kattis id | tests | bytes | family | route clue |
|---|---:|---:|---:|---:|---:|---|---|
| `submission_1181_81.95_7.cpp` | 81.95 | 81.945906 | 19912924 | 7 | 131019 | central_s3b16_plus_broad_b16 | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1180_81.94_7.cpp` | 81.94 | 81.941905 | 19912852 | 7 | 130905 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1177_81.94_7.cpp` | 81.94 | 81.941571 | 19912823 | 7 | 130864 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1188_81.94_7.cpp` | 81.94 | 81.941571 | 19913102 | 7 | 130864 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1179_81.94_7.cpp` | 81.94 | 81.940905 | 19912850 | 7 | 130864 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1182_81.94_7.cpp` | 81.94 | 81.940905 | 19912927 | 7 | 130864 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1183_81.94_7.cpp` | 81.94 | 81.940905 | 19912946 | 7 | 130865 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1186_81.94_7.cpp` | 81.94 | 81.940905 | 19912973 | 7 | 130864 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1178_81.94_7.cpp` | 81.94 | 81.940571 | 19912824 | 7 | 130864 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1042_81.94_7.cpp` | 81.94 |  |  | 7 | 128611 | central_s3b16_bestline | GN() > W2G::run())W2C::run() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.05) > B16::R(39000,60000,76,... |
| `submission_1043_81.94_7.cpp` | 81.94 |  |  | 7 | 127692 | central_s3b16_bestline | GN() > W2G::run())W2C::run() > VIMP::run() > MIDEC::run() > JD() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.05) > B16::R(39000,60... |
| `submission_1044_81.94_7.cpp` | 81.94 |  |  | 7 | 122675 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1074_81.94_7.cpp` | 81.94 |  |  | 7 | 127298 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1077_81.94_7.cpp` | 81.94 |  |  | 7 | 127918 | central_s3b16_bestline | GN() > JD() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,... |
| `submission_1078_81.94_7.cpp` | 81.94 |  |  | 7 | 127477 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1079_81.94_7.cpp` | 81.94 |  |  | 7 | 127469 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1089_81.94_7.cpp` | 81.94 |  |  | 7 | 125820 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > JD() > WK::run() > B16::R(39000,60000,220,-7,192,... |
| `submission_1095_81.94_7.cpp` | 81.94 |  |  | 7 | 125794 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1096_81.94_7.cpp` | 81.94 |  |  | 7 | 123622 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |
| `submission_1098_81.94_7.cpp` | 81.94 |  |  | 7 | 116347 | central_s3b16_bestline | GN() > if(!W2G::run())W2C::run() > W5::post_patch_pass() > VIMP::run() > MIDEC::run() > WK::run() > B16::R(39000,60000,220,-7,192,.96,18.... |

## Family Clusters

| family label | count | best filename score | test buckets | dominant score tokens | representative high files |
|---|---:|---:|---|---|---|
| direct_primitive_or_structural_probe | 69 | 81.93 | 0:2, 1:2, 2:13, 3:9, 4:16, 5:3, 6:7, 7:17 | 16.30:6, 0.00:6, 81.91:5, 81.93:4, 17.00:4, 17.01:3, 81.69:3, 15.41:2, 27.62:2, 45.18:2, 16.49:2, 16.53:2 | `submission_1005_81.93_7.cpp`, `submission_1006_81.93_7.cpp`, `submission_1010_81.93_7.cpp`, `submission_1103_81.93_7.cpp`, `submission_1062_81.91_7.cpp` |
| b16_wk_plateau | 58 | 81.93 | 6:5, 7:53 | 81.93:49, 70.49:3, 79.88:1, 81.60:1, 10.39:1, 67.73:1, 68.11:1, 80.65:1 | `submission_968_81.93_7.cpp`, `submission_971_81.93_7.cpp`, `submission_973_81.93_7.cpp`, `submission_978_81.93_7.cpp`, `submission_980_81.93_7.cpp` |
| small_standalone_worker_code | 45 | 43.88 | 0:2, 1:2, 2:16, 3:2, 4:10, 5:5, 6:1, 7:7 | 16.47:5, 0.00:4, 16.56:3, 16.52:2, 15.58:2, 1.51:2, 16.50:2, 16.54:2, 16.53:2, 30.64:1, 43.88:1, 16.30:1 | `submission_1000_43.88_6.cpp`, `submission_1166_42.50_4.cpp`, `submission_1168_40.88_5.cpp`, `submission_1165_40.65_4.cpp`, `submission_997_30.64_4.cpp` |
| standalone_qem_or_generic_rewrite | 44 | 64.44 | 0:3, 2:7, 3:8, 4:17, 5:4, 6:3, 7:2 | 29.33:6, 48.10:3, 0.00:3, 46.08:2, 28.83:2, 41.89:2, 40.45:2, 16.54:2, 38.28:1, 15.67:1, 59.78:1, 38.26:1 | `submission_1220_64.44_6.cpp`, `submission_975_59.78_5.cpp`, `submission_1233_55.27_5.cpp`, `submission_1216_55.23_5.cpp`, `submission_1167_50.25_5.cpp` |
| central_s3b16_bestline | 34 | 81.94 | 5:5, 7:29 | 81.94:27, 81.61:2, 56.31:2, 56.66:2, 52.61:1 | `submission_1042_81.94_7.cpp`, `submission_1043_81.94_7.cpp`, `submission_1044_81.94_7.cpp`, `submission_1074_81.94_7.cpp`, `submission_1077_81.94_7.cpp` |
| central_s3b16_plus_broad_b16 | 10 | 81.95 | 2:1, 5:7, 6:1, 7:1 | 53.93:6, 81.95:1, 56.93:1, 15.14:1, 68.11:1 | `submission_1181_81.95_7.cpp`, `submission_1236_68.11_6.cpp`, `submission_1194_56.93_5.cpp`, `submission_1145_53.93_5.cpp`, `submission_1189_53.93_5.cpp` |
| box_gate_plus_b16_wk_plateau | 6 | 81.93 | 6:1, 7:5 | 81.93:5, 67.75:1 | `submission_1035_81.93_7.cpp`, `submission_1040_81.93_7.cpp`, `submission_1197_81.93_7.cpp`, `submission_1206_81.93_7.cpp`, `submission_1219_81.93_7.cpp` |
| other_large_solver | 3 | 67.75 | 0:1, 2:1, 6:1 | 67.75:1, 12.91:1, 0.00:1 | `submission_1105_67.75_6.cpp`, `submission_1190_12.91_2.cpp`, `submission_1225_0.00_0.cpp` |
| standalone_axial_impostor_or_cubemap | 3 | 44.20 | 4:1, 7:2 | 0.00:2, 44.20:1 | `submission_1232_44.20_4.cpp`, `submission_1201_0.00_7.cpp`, `submission_1202_0.00_7.cpp` |

## Concrete Log / Commit Clues

- `688da89 Record round 84 Kattis submissions` added `submission_1177..1186` and fetched exact files including `kattis_19912924_81.945906_7.cpp`; by filename/order, this maps to `submission_1181_81.95_7.cpp`.
- `ce527d7 Record 81.9459 central variant batch` names the same central-variant high-water batch.
- `worker_outputs/round83_targeted_feedback_prompt.md` identifies `submission_1044_81.94_7.cpp` / `kattis_19903544_81.938904.cpp` as the earlier best and explicitly preserves the main route: `GN(); if(!W2G::run()) W2C::run(); W5; VIMP; MIDEC; WK; B16; B16; S3B16; WK tails; JD();`.
- `local_worker_S12_score_mining/S12_score_mining_manifest.md` records the older `19901232 = 81.93457` plateau and shows nearby B16 parameter variants staying rounded-plateau, while broad high-N and S3 tiny routes fell to roughly 67.74/70.49.
- `local_worker_S8_proxy_eval/MANIFEST.md`, `local_worker_S14_case5_next/RESULTS_S14.md`, `local_worker_BROAD_09/BROAD09_NOTES.md`, and `local_worker_BROAD_30/README_BROAD_30.md` repeatedly point at the `49843..50625` UV/torus/ripple band, but with orientation/time/proxy-risk caveats.
- B91/B92/B93/B97/SGX manifests show many attempted structural branches: AABB/box, packed structural recognizers, shadow carrier, direct torus/sphere/tube, cubemap/grid. The submitted/fetched history around those attempts mostly ties `81.93` or drops to low buckets, so they are evidence for careful fail-closed integration, not for replacing the core solver.

## What Likely Blocks 91.80

1. The current solver is overfit to a stable 7/7 validity pipeline but not to the hidden visual/objective optimum. The plateau files preserve manifold/Hausdorff/SSIM well enough but fail to remove enough vertices on at least one high-weight hidden family.
2. Source-size and time budgets are severe. The high-water `submission_1181` is 131019 bytes, only 53 bytes under the 131072-byte source limit, so adding a real new branch requires packing/pruning or an in-place replacement.
3. Broad structural rewrites are too brittle. The root distribution and fetched logs show many sample-pass/generic candidates at 0, 16-17, 29, 40-56, 67-70; they either lose hidden validity/SSIM or do not trigger safely.
4. The 47500-60000 band is still the most evidence-rich suspect, but proxy improvements have not reliably translated to exact Kattis improvement. Orientation conflicts and local-vs-official SSIM mismatch are recurring risk flags.
5. The exact `81.945906` improvement seems to come from a very narrow central variant interaction, not a broadly understood new method. Until the delta from `1044/1177/1181/1206` is isolated, workers risk unknowingly deleting the only positive signal.

## Prioritized Hypotheses / Next Experiments

1. **Diff and isolate submission_1181 / Kattis 19912924 against 1044, 1177, and 1199/1206.** Only exact score above the rounded plateau found in git metadata: 19912924 = 81.945906. Route includes central S3B16 plus two broad B16::R(8000,1100000,...) calls; later round 90-93 plateau submissions are back to 81.93457.
2. **Use a packed or pruned exact-best base before adding any new branch.** High-water root source is 131019 bytes, leaving only 53 bytes under the 131072 source cap; B92-style packing manifests show 17-24 KB possible slack but those packed branches need exact-best rebasing.
3. **Target the 47500-60000 case5/ripple band with a fail-closed local tournament, not a standalone remesher.** S8/S14/BROAD_09/BROAD_30 proxy reports repeatedly point to uv/torus/ripple cases around 49843-50625; standalone structural rewrites usually score 0-56, while plateau-preserving guarded inserts keep 7/7.
4. **Treat accepted-count as a weak signal; require exact score metadata or paired A/B submissions.** Root has 114 submissions with 7 accepted tests, including scores from 0.00 through 81.95; many 7/7 submissions score 16.30, 42.12, 63.67, or 81.93.
5. **Blacklist generic QEM/primitive/cubemap/hidden-carrier replacements unless wrapped as fail-closed exact-best branches.** Root score buckets and manifests show QEM/standalone remesh/primitive attempts mostly in 0-56 or tied 81.93; they do not close the 9.85-point gap to 91.80.

## Machine-Readable Companion

See `history_gap.json` for the full parsed distribution, top submissions, route signatures, exact SHA clusters, fetched Kattis score metadata, and experiment list.
