# 2026-07-09 Armadillo exact-23 repair experiment

Context: local worker ceiling is 16; browser/Pro Extended chat management is paused. Active high-water remains `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact `81.978181`.

## Proxy construction

External literature/search confirmed a common Armadillo mesh size of `23201` vertices and `46398` facets, matching hidden test-4 count evidence. The first public raw Armadillo downloaded from `common-3d-test-models` was `49990/99976`, so it was normalized and decimated with `fast_simplification` to exact `23201/46398` proxies at agg `{3,5,7,10,15,20}`.

## Discovery

On exact Armadillo proxies, `GN()` frequently creates a compressed `11601/23198` active mesh that final `JD()/HS()` rejects, causing the official output to fall back to the original `23201/46398`. Dumping raw post-`GN` state for `agg3` showed:

- no degenerate faces
- one duplicate face-key
- one edge with incidence `4`
- the duplicate pair used an isolated degree-2 vertex

Removing both opposite duplicate cap faces and compacting the isolated degree-2 vertex restored a valid closed manifold while preserving `F = 2V - 4`.

## Local repaired proxy results

Final `RP()` cap-pair repair after the normal pipeline:

| Proxy | Repaired first line | `vps_eval512` |
| --- | ---: | ---: |
| `agg3` | `11600 23196` | `0.954637628481` |
| `agg5` | `11601 23198` | `0.956342353946` |
| `agg7` | `11599 23194` | `0.960040597173` |
| `agg10` | `11597 23190` | `0.961526388699` |
| `agg15` | `11597 23190` | `0.961551649584` |
| `agg20` | `11597 23190` | `0.961551649584` |

This is a real local breakthrough for an Armadillo-like exact-23 case, but the source is highly layout/timing sensitive.

## Source-limit attempts

- Direct repair body before `B16/main`: over limit and changed non-exact behavior badly.
- Removing `W2G` to fit: under limit, but non-exact behavior regressed badly (`Nefertiti` output around `895867` vertices).
- Body after `main` with forward declaration: over limit; local behavior was closer to base (`Lucy` byte-identical, `Nefertiti` `150190` vs base `150070`).
- Macro-shave candidate `submission_1579_arm23_repair_macrofit.cpp`: under limit at `131000` bytes and preserved `Lucy` first line, but still changed Nefertiti timing (`149982/299960` locally).

## Kattis result

Submitted macrofit as `19928275`, recorded as `submission_1579_68.11_6.cpp`:

- Judgement: `Accepted (68.113410)`
- Tests: `6/7`
- Runtime: `20.30s`
- Detail page: `Test case 4/7: SSIM is too low`

This is far below the high-water and not a base candidate.

## Decision

Do not submit more naive `RP()` macrofit/source-shaved variants. The cap-pair repair idea is valuable, but it must be integrated in a way that is byte/layout-neutral against `1448` on high-N timing-sensitive cases, especially Nefertiti/high2. A viable next attempt would need either:

1. an exact-23 early-return branch that avoids touching the high-N hot path, or
2. a preprocessor/macro-only source reduction preserving preprocessed hot tokens and measured high2 behavior before submission.
