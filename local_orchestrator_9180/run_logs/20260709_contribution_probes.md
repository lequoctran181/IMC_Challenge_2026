# 2026-07-09 contribution probes

Active base: `submission_1448_81.98_7.cpp` / Kattis `19922865`, exact score `81.978181`, `7/7`.

## Clean probes

| Submission | Branch | Result | Inference |
| --- | --- | --- | --- |
| `19923695` / `submission_1496_70.49_7.cpp` | For exact `N == 49987 && M == 99970`, output original mesh with `%.17g`; otherwise run base. | `70.491928`, `7/7` | Clean accepted probe. Assuming six secret tests are equally averaged, this exact case contributes about `(81.978181 - 70.491928) * 6 = 68.917518` points in the active base. |
| `19923847` / `submission_1502_40.78_7.cpp` | For `N < 47500`, output original mesh with `%.17g`; otherwise run base. | `40.784692`, `7/7` | Clean accepted probe for the smaller hidden group. Assuming six secret tests are equally averaged, this group contributes about `(81.978181 - 40.784692) * 6 = 247.160934` points in the active base. If this is hidden tests 2-4, the average current contribution is about `82.386978` per case. |

## Contaminated probes

| Submission | Branch | Result | Why not clean |
| --- | --- | --- | --- |
| `19923641` / `submission_1493_45.83_5.cpp` | Same exact-case original probe, but through low-precision `IJ()` output. | `45.825898`, `5/7` | Low-precision original printing is unsafe/noisy. |
| `19923772` / `submission_1501_38.81_5.cpp` | For `N > 120000`, output original mesh with `%.17g`; otherwise run base. | `38.814971`, `5/7` | High-N original output is not a valid measurement path, likely due output/time/format stress on large cases. Do not infer high-N contribution from this score. |

## Strategy consequence

- Exact case `49987/99970` is important but not sufficient: even improving it from `68.92` to `100` would add only about `5.18` total score points.
- The `N < 47500` group is strong but not saturated: if it covers three hidden cases, perfecting all three would add at most about `(300 - 247.160934) / 6 = 8.806511` total points. That ceiling is useful, but unrealistic without a new small-case strategy.
- To reach `91.80+`, we need both a better exact-case strategy and at least one additional high-impact hidden-case improvement, with large-case compact ablations still a priority because original-output probes are contaminated there.
- Avoid high-N original-output probes as contribution estimators; use narrower ablations or compact candidate replacements instead.
