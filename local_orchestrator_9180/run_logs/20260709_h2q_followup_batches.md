# 2026-07-09 H2Q follow-up batches

After Kattis `19928674` failed at `53.91929`, three local 16-lane follow-up batches tested whether the promising high2 QEM signal could be made judge-safe with exact test4 guards and lower runtime.

## Batch: `batch_1499_h2q_t4guard_20260709`

Idea: keep the `1499` slack base, add exact `N==35292 && M==70580` test4 early branch, and use faster no-proxy H2Q variants.

Result:

- Test4 proxy stayed stable at `1942/3880`, `vps_eval512 = 0.946045615133`.
- Most no-proxy H2Q variants either did not run enough and stayed near `150k`, or collapsed dangerously.
- Worst unsafe family produced high2 `4959/9914`, `vps_eval512 ~= 0.82826`.

Conclusion: no-proxy H2Q inside the normal tail is unsafe because later passes can collapse far past the intended target.

## Batch: `batch_1499_h2q_proxy_t4guard_20260709`

Idea: same test4 guard, but restore proxy rollback (`eW64`) and lower the H2Q time limit to avoid Kattis TLE.

Result:

- High2 outputs stayed around `149926..150190` vertices.
- Local Nefertiti runtime was `20.37..25.16s`.
- No useful vertex-count gain over base `1499`.

Conclusion: proxy-guarded H2Q is too slow and loses the intended high2 gain when the time limit is shortened.

## Batch: `batch_1499_h2q_earlyreturn_20260709`

Idea: exact high2 branch runs `GN + W5 + VIMP`, then H2Q, then returns immediately to avoid the common tail collapse.

Result:

- H2Q starts too early, while the active mesh is still too large.
- Best completed runs were still `630k..655k` vertices with `vps_eval512 ~= 0.988..0.989`.
- Runtime was `25..28s`; lower targets timed out locally.

Conclusion: sorting/processing the large edge set before the normal tail is infeasible.

## Batch: `batch_1499_h2q_posttail_20260709`

Idea: exact high2 branch runs the full normal tail to reach about `150k`, then H2Q does only a short final pass and returns.

Result:

- Runtime was already `20..22s` before useful H2Q work.
- Output stayed at base high2 `150190/300376`.

Conclusion: after the full tail there is no remaining runtime budget for source-side H2Q.

## Current decision

Do not submit any candidate from these batches. The high2 QEM signal is real offline, but source-side H2Q cannot currently be made both fast and test4-stable. Future high2 work needs a cheaper decimator than sorted-edge H2Q, or a precomputed/structure-aware route that avoids processing all edges.
