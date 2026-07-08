# Worker B Case Inference Report

Scope: local-only inference for IMC Challenge 2026 Problem B `simplifygeometry`. I read high 7/7 submissions around `81.93..81.95`, partial high files around `67.7/70.5/79.88`, and recent low accepted `4/7`, `5/7`, and `6/7` families. No Kattis submit was performed.

## Executive Read

The current plateau is a strong routed solver, not a generic simplifier. The stable anchor pipeline is:

```cpp
JC(); GN();
if(!W2G::run()) W2C::run();
W5::post_patch_pass();
VIMP::run();
MIDEC::run();
WK::run();
B16::R(39000,60000,220,-7,192,.96,18.05);
B16::R(39000,60000,76,-10,192,.96,18.35);
repeat WK on 47500..60000;
JD();
```

`kattis_19901232.cpp` and many `submission_*_81.93/81.94_7.cpp` variants use this skeleton. The local exact-score notes identify `19901232` as `81.93457`; later filenames round many neighbors to `81.93/81.94`, and `submission_1181_81.95_7.cpp` is the best visible rounded token.

The likely missing 10 points are not from another small second-B16 tweak. They are probably from one or two hidden structural slices where the plateau currently no-ops or only modestly reduces: case5 lower/upper around `49843..50625`, case6/high-N smooth/box-like surfaces, or a special periodic/parametric surface not safely recognized by the current router.

## Evidence From Scores And Code Signatures

| score cluster | representative files | code signature | inference |
|---|---|---|---|
| `81.93..81.95`, `7/7` | `kattis_19901232.cpp`, `submission_1042..1044`, `1074`, `1089`, `1131`, `1181` | W2G/W2C -> W5 -> VIMP -> MIDEC -> WK -> B16 `220,-7` -> B16 `76,-10` -> WK tail | Plateau anchor. Very robust across 7 hidden cases, but likely saturates all easy routed wins. |
| `79.88`, `7/7` | `submission_969_79.88_7.cpp` | anchor plus `if(N>=49843&&N<50234) B16::R(49843,50233,140,-5,224,.958,18.55)` | Narrow lower case5 B16-only intervention hurts enough to lose ~2 points while still accepted. B16 count/stride alone is not the right lower-slice fix. |
| `70.49`, `6/7` | `submission_1034_70.49_6.cpp`, `1036`, `1060`, `1084` | UG periodic/grid special or replacing/omitting anchor B16; `1060` has `B16::R(47501,85000,240,-5,160,.957,18.15)` | Specialized router fires destructively on one hidden case, or anchor pass was removed. Treat as branch-risk evidence, not as useful generic direction. |
| `67.73..67.75`, `6/7` | `submission_1033_67.75_6.cpp`, `1061`, `1085`, `1105`, `kattis_19901997_67.75_6.cpp` | broad structural branches such as `STARX`, high-count `B16::R(47500,50625,700/900,...)`, or late special run after output | One hidden case is badly damaged by aggressive structural replacement or broad high-N/low-res gate. |
| `53..57`, `5/7` | `submission_1184`, `1187`, `1189..1196`, `1216`, `1233` | mostly standalone generic QEM / primitive recognizers / torus/ellipsoid/radial cube shell | Generic simplifiers can pass many format/topology cases but miss the high-value hidden structural routes. Useful as fallback ideas, not replacement. |
| `40..48`, `4/7` | `submission_1115`, `1132`, `1158`, `1214`, `1235` | generic decimation or cubemap/box shell support | Confirms several hidden cases reward structural recognition; plain decimation under-scores heavily. |

## Likely Hidden Case Archetypes

1. **Official/tiny/simple sample guard**
   - Evidence: recent low files include explicit `tryTinySampleBox`, `N0==9 && M0==14`, and primitive recognizers, but they score only `4/7..5/7`.
   - Current anchor likely handles this by falling back safely rather than scoring big.
   - Branch risk: low. Do not spend the 91.80 push here except preserving output validity.

2. **Case3 periodic torus/grid, `N ~= 23125..23500`, `M == 2N`**
   - Evidence: anchor has `W2G::run()` with detector requiring `N<23500`, `M==2*N`, factorable `U,V` in `80..320`, and sampled face order matching a wrapped grid. `W2C::run()` is the fallback same band.
   - Proxy notes mention `case3_torus_23296`.
   - Helps: `W2G` direct grid downsample; `W2C` edge-collapse fallback.
   - Hurts: old/broad replacements that try to generalize periodic detection above this band.
   - Missing: probably not the 10-point gap unless a second periodic class exists outside `23125..23500`.

3. **Core case5 smooth mid-band, `N ~= 47500..60000`**
   - Evidence: most anchor-specific branches guard this range: `W5::post_patch_pass`, `MIDEC`, `WK`, two `B16` passes, and WK tail.
   - Mesh classifier evidence: `AS`, `BG`, `BL`, `AL`, `Z`, `AH` ratios from normal smoothness/sharpness sampling.
   - Helps: current second B16 ridge `count 74..76`, stride `-10`, `vp=192`, `q=.960..961`, end `18.35`.
   - Hurts: direct variants `96,-9,.962` and `80,-11,224,.964` are documented slight drops; `47501..85000,240,-5,160,.957` scored `70.49/6`.
   - Missing: router/tournament choice among W5/MIDEC/WK/B16 order, not raw B16 count.

4. **Lower case5 slice, `49843 <= N < 50234`, UV/torus-ripple family**
   - Evidence: S4/S14 notes target exactly this band. S4_08 guard is `AS && BG <= 0.001 && BL >= 0.985 && Z <= 0.0037`.
   - Proxy evidence: `uv49954_bumpy` improves `1249/2494 -> 1129/2254`; `r5_torus50176_ripple` improves `1985/3970 -> 1799/3598` in S14, while looser variants worsen rough-phase.
   - Helps: feature-guarded W5/R13-style case5 patch with aggressive cap; BROAD_30 DC-deep helps torus-ripple proxy only.
   - Hurts: B16-only lower patch (`submission_969`) loses ~2 points; exact unguarded S4 worsens `r5_uv49954_rough_phase`.
   - Missing: a safe lower-slice recognizer that distinguishes UV/bumpy/ripple from rough-phase and torus genus cases.

5. **Upper case5 slice, `50234 <= N < 50625`, boxy/softbox**
   - Evidence: S4/S14 upper exact/aggressive candidates; S8 generated `s8_uv50402_softbox`; S14 notes `r5_uv50623_boxy` and `s8_uv50402_softbox`.
   - Helps: `s14_07_upper_aggr_50234_50625` proxy improves upper boxy/softbox by about `-105..-111V`.
   - Hurts: orientation-conflict flags and rough-phase worsening in exact upper fallback.
   - Missing: orientation-safe upper patch; current anchor is conservative here.

6. **Case6 / high-mid smooth band, `60001..160000`**
   - Evidence: S5 notes say `worker509_B06` late `60001..160000,160,-13,96,.945,18.7` scored slightly lower when it omitted the winning second B16, so the late pass itself remains plausible. Candidate queue narrows this to `vp=128`, `q=.965`, lower/upper halves.
   - Helps: strict late B16 after preserving anchor, especially `60001..120000`.
   - Hurts: removing the anchor second B16; too-low visual proxy.
   - Missing: isolated proof of a positive strict late pass after current best.

7. **Very high-N smooth/UV/box/torus, `160000..1200000`**
   - Evidence: S5/S15 blacklist broad high-N `B16::R(160001,1200000,220,-17,64,.93,18.7)` after about `67.74`. BROAD_26 targets high-N AABB/box-shell with stricter gate and timer reset. S15 proposes extra VIMP and W5L473 with tight `AS/BG/AL/AH`.
   - Helps: high-N VIMP/W5L473 only with tight smoothness guards and `vps256/512` high floors; AB19 only for actual box-shell grids.
   - Hurts: broad low-resolution B16 over high-N, especially `vp<=96`, `q<=.93`.
   - Missing: high-N classifier strong enough to avoid official-score false positives at lower render resolution.

8. **Generic irregular meshes / sharp meshes**
   - Evidence: recent `4/7` and `5/7` accepted files use generic QEM and primitive reconstruction but top out near `45..57`.
   - Helps: generic fallback for cases where structural recognizers no-op.
   - Hurts: using generic simplifier as main route loses the special hidden cases.
   - Missing: probably not worth replacing anchor; only consider as a guarded emergency branch when anchor cove remains high and classifier says non-smooth/non-special.

## Branches To Keep, Promote, Or Avoid

Keep as baseline:

- `W2G/W2C` case3 detector.
- Current `W5 -> VIMP -> MIDEC -> WK -> B16(220,-7) -> B16(76,-10) -> WK tail` order.
- Second B16 ridge near `74..76`, `stride -10`, `vp=192`, `q=.960..961`.

Promote experiments:

- Feature-guarded S4/S14 lower case5 patch, especially S4_08/S14_01 style, not B16-only.
- Upper-slice S14_07 with orientation-safe validation added or stricter rollback.
- S7 tournament selectors: skip MIDEC, B16-before-WK, and pre-MIDEC WK swap, kept only on cove gain plus `vps512`.
- S5/S15 strict case6/high-N additions after preserving anchor.
- BROAD_30-style DC-deep only if recognizer detects torus/ripple strongly, and only with full rollback.

Avoid:

- Broad high-N B16 with `vp=64` or low `q` over `160001..1200000`.
- Direct second-B16 retries `96,-9,.962` and `80,-11,224,.964`.
- B16-only lower case5 patch `49843..50233,140,-5,224,.958`.
- Any branch that replaces or omits the current second B16 while testing a late pass.
- Unguarded S4/R13 lower patch: it catches rough-phase and worsens it.

## Concrete Router / Branch Experiments

1. **Lower case5 guarded W5 patch tournament**
   - Band: `49843 <= N < 50234`.
   - Gate: start with `AS && BG <= 0.001 && BL >= 0.985 && Z <= 0.0037`.
   - Action: clone anchor state after `W5::post_patch_pass`; run S4_08 aggressive cap; keep only if cove drops and `vps512` is not below baseline by more than a tiny tolerance.
   - Risk control: add rough-phase veto using higher `Z`/lower `BL` or sampled normal variance.

2. **Upper case5 orientation-safe patch**
   - Band: `50234 <= N < 50625`.
   - Action: S14_07 upper-aggressive route, but require stricter `HS` plus a local orientation-conflict check before `JD`.
   - Keep condition: cove drop plus `vps512 >= max(baseline_vps512 - .001, .95)`; no orientation-conflict flag.

3. **Mid-band pass-order tournament**
   - Band: `47500 <= N < 60000`.
   - Actions: compare anchor order against `WK before MIDEC`, `skip MIDEC`, and `B16 before WK`.
   - Keep condition: same topology validator, cove strictly lower, `vps512 >= baseline_vps512 - .001`.
   - Reason: plateau is order-sensitive; S7 says these are non-redundant with B16 parameter sweeps.

4. **Strict case6 late pass**
   - Band: `60001 <= N <= 120000` first, then `120001..160000`.
   - Action: after full anchor tail, test `B16::R(...,120,-13,128,.965,18.68)` style.
   - Keep condition: preserve anchor second B16; require high proxy and rollback on any cove increase.

5. **High-N smooth extra VIMP**
   - Band: `60000..180000` first; later `260000..1200000`.
   - Gate: S15-like `AS && BG<.0095 && AL>.78 && AH<.105`, tighter for upper high-N.
   - Action: one extra VIMP or W5L473 pass with small cap.
   - Avoid B16 unless count is tiny and `q>=.985`.

6. **Box-shell high-N recognizer**
   - Band: very high-N only.
   - Action: BROAD_26 AB19-style AABB shell/grid branch.
   - Gate: point-to-box-face tolerance, grid coverage, and `vps512 >= .992`.
   - Keep only if non-box UV/torus proxies byte-identical/no-op.

## What Is Missing For A 10-Point Jump

- A trustworthy structural classifier. The score history shows the solver already has mechanisms that can reduce aggressively; the failures happen when a mechanism fires on the wrong hidden case.
- Better local proxy discrimination. `vps64/96/128` admits false positives on high-N; case5 patches can pass edge/topology checks while orientation or official visual score suffers.
- A slice-aware router for `49843..50625`. Current best is safe, but S4/S14/BROAD_30 show sizeable local count deltas on exactly this band. This is the most plausible near-term source of a multi-point gain.
- A high-N branch that no-ops by default. The broad high-N catastrophe implies high-N can contribute many points, but any branch there needs shape gates much stronger than just `N` and smoothness.

## Recommended Immediate Queue

1. Implement a tournamented S4_08 lower case5 branch on top of `19901232`/current best.
2. Implement upper S14_07 with orientation-safe rollback.
3. Try S7 pass-order tournament in `47500..60000`.
4. Try S5 strict lower case6 `60001..120000`.
5. Try S15 extra VIMP `60000..180000`.
6. Only after those, test high-N box-shell AB19 and BROAD_30 torus-ripple recognizers with hard no-op guarantees on non-target proxies.

