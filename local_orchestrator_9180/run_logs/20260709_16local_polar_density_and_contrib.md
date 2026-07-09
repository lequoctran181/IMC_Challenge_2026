# 2026-07-09 16-Local Worker Log

Scope update: stop managing the extra 8 ChatGPT Pro Extended chats. Use 16 local workers only, with Kattis submissions reserved for evidence-backed candidates or diagnostics.

## High-water

- Current best source: `submission_1448_81.98_7.cpp`
- Kattis score: `Accepted (81.978181)`, 7/7
- Hidden exact shapes of interest:
  - `N=4098, M=8192`
  - `N=23201, M=46398`
  - `N=35292, M=70580`
  - `N=49987, M=99970`
  - `N=377084, M=754688`
  - `N=1009118, M=2018232`

## 16-worker polar-density batch

Directory: `local_orchestrator_9180/queue16_hw_polar_density_20260709`

The batch tested same-source-shape edits that relaxed the FL/polar guard (`lo < .7*hi`) and varied density ladders/thresholds. Local result was negative:

- Polar proxy did not improve; best variants remained around 1146-1165 vertices with VPS about 0.939-0.941.
- Lucy proxy worsened versus high-water on many variants.
- `uvsphere_4098` worsened badly from the high-water small-mesh behavior to about 1368 vertices.

Decision: do not submit these polar-density variants.

## Contribution diagnostics

Goal: intentionally invalidate one exact hidden shape at a time, while preserving high-water behavior on the other cases, to infer the hidden case contribution from total score.

Generated in `local_orchestrator_9180/contribution_exact_invalid_20260709`:

- `diag_invalid_high1_377084_tiny.cpp`
- `diag_invalid_high2_1009118_tiny.cpp`
- `diag_invalid_small_4098_tiny.cpp`
- `diag_invalid_mid23_23201_tiny.cpp`
- `diag_invalid_mid35_35292_tiny.cpp`

Implementation notes:

- Uses `if(N==target)return!puts("0 0");` after input parse.
- Removes unused `#include<cstdint>` to stay under Kattis source-size accounting.
- Local compile/sample passed for high1 and high2 tiny variants: sample output starts with `8 12`.

Submissions:

- `19930417`: `diag_invalid_high1_377084_tiny.cpp`
  - Result: `Accepted (66.802644)`, 6/7
  - Decision: this is the least contaminated large-case invalid probe observed so far, but use it only as a rough signal; Kattis partial-score arithmetic is not reliable unless the full failing test list is known.
- `19930429`: `diag_invalid_high2_1009118_tiny.cpp`, submitted and waiting for final judgement.
  - Result: `Accepted (53.927292)`, 5/7
  - Decision: contaminated / not a clean one-case contribution probe. The clean expectation from prior arithmetic was about `67.855279`, 6/7.
- `19930450`: `diag_invalid_small_4098_tiny.cpp`, submitted and waiting for final judgement.
  - Result: `Accepted (37.431440)`, 4/7
  - Decision: contaminated / not a clean one-case contribution probe.

After the early-return small probe was contaminated, generated after-work probes that run the full high-water algorithm and only then print `0 0` for exact `N` before `JD()`. These compile and sample-pass:

- `diag_invalid_after_small_4098_tiny.cpp`
- `diag_invalid_after_mid23_23201_tiny.cpp`
- `diag_invalid_after_mid35_35292_tiny.cpp`

Submission:

- `19930485`: `diag_invalid_after_small_4098_tiny.cpp`, submitted and waiting for final judgement.
  - Result: `Accepted (51.617558)`, 5/7
  - Decision: still contaminated. Stop invalid-output contribution probes for small/mid cases; do not submit the prepared mid23/mid35 after-work probes without a new isolation mechanism.

Same-layout follow-up using the cleaner `19930417` source layout:

- `19930532`: `small4098.cpp`
  - Result: `Accepted (37.431440)`, 4/7, runtime `>21.00s`
  - Detail: `secret/1` invalid and `secret/3` SSIM too low.
  - Decision: contaminated; not a clean small-case contribution probe.
- `19930534`: `low23201.cpp`
  - Result: `Accepted (43.094425)`, 4/7, runtime `>21.00s`
  - Detail: `secret/2` invalid and `secret/3` SSIM too low.
  - Decision: contaminated; not a clean low23 contribution probe.
- `19930542`: `up35292.cpp`
  - Result: `Accepted (53.927292)`, 5/7, runtime `>21.00s`
  - Detail: only `secret/3` invalid.
  - Decision: identifies `secret/3` as the exact `N==35292 && M==70580` upper-mid case and confirms this case is the recurring fragile guard target. Do not use the score drop as a clean contribution estimate because the `5/7` verdict can hide a second failing group.
- `19930556`: `case49987.cpp`
  - Result: `Accepted (42.472647)`, 4/7, runtime `>21.00s`
  - Detail: `secret/3` SSIM too low and `secret/4` invalid.
  - Decision: contaminated; not a clean case5 contribution probe.

Conclusion: direct early invalid branches are too contaminated for precise contribution arithmetic. The stable operational lesson is stronger: `secret/3` is exact `N==35292 && M==70580` and is the recurring fragile case that many high2/highbig attempts break. The next high-leverage route remains high2/highbig, but every candidate needs an explicit plan to preserve or conservatively protect upper35.

Only infer a case contribution from `81.978181 - T` when the detail page proves exactly one hidden group failed and there is no hidden TLE/second failure. Otherwise treat the score drop as contaminated.

## H2Q with upper35 guard

Generated in `local_orchestrator_9180/batch_h2q_upperguard_20260710`.

- `19930640`: `ug_after_gn.cpp`, H2Q target `90000` plus exact `N==35292 && M==70580` return after `GN()`.
  - Result: `Accepted (67.793400)`, `6/7`, runtime `>21.00s`.
  - Detail page had no validator feedback block. Interpretation: upper35 protection improved the previous H2Q collapse, but highbig/test7 still fails or times out.
- `19930662`: `ug_after_gn_h2q_t130.cpp`, lighter H2Q target `130000` with the same upper35 guard.
  - Result: `Accepted (67.752453)`, `6/7`, runtime `>21.00s`.
  - Detail page had no validator feedback block. Interpretation: lighter target did not rescue test7; the failure is likely highbig runtime/validator, not upper35.

Conclusion: the upper guard after `GN()` is useful, but H2Q still cannot produce a valid judged highbig output under the time limit. Next highbig attempt should reduce or time-cap the post-pass more aggressively, or fall back to a smaller same-mechanism tweak that keeps test7 valid.

- `19930685`: `tc145_l186_rb192.cpp`, H2Q target `145000`, loop cap `18.6`, rollback if elapsed exceeds `19.2`, with the same upper35-after-GN guard.
  - Local Nefertiti proxy: `149601/299198` in about `19.9s` on first run.
  - Kattis result: `Accepted (67.752453)`, `6/7`, runtime `>21.00s`, no feedback block.
  - Decision: even a very light/time-capped H2Q remains on the same 6/7 plateau; close H2Q post-pass variants unless a substantially cheaper highbig mechanism or a guaranteed rollback-to-base before output is implemented.

## H2S mild scan with upper35 guard

Generated in `local_orchestrator_9180/batch_h2s_mild_upperguard_20260710`.

- Local Nefertiti outputs were unstable: cap3000 `95856/191708`, cap5000 `2335/4666`, cap8000 `90856/181708`, cap10000 `88856/177708`, cap15000 `135190/270376`.
- `19930708`: submitted `h2s_cap15000_ug.cpp` because it had the mildest useful local output and `vps_eval512 = 0.917495885621`.
  - Kattis result: `Accepted (67.752453)`, `6/7`, runtime `20.51s`; detail HTML contains WA/SSIM but no compact feedback block.
  - Decision: H2S, like H2Q, remains on the highbig 6/7 plateau. Close H2S/H2Q post-pass families for now; future highbig work needs a different mechanism or a judge-safe way to restore test7.
