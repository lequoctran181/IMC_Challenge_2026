# 2026-07-09 primitive-threshold same-size batch

Base: `submission_1448_81.98_7.cpp`.

Purpose: test whether existing primitive branches (`GI`, `FS`, `GM`) can be made useful for the exact sphere-topology hidden buckets by same-size numeric edits only, avoiding new source layout risk.

Generated `local_orchestrator_9180/queue16_primitive_threshold_20260709/` with 16 variants, all size `131030` bytes. The variants lowered/changed:

- `GI` threshold and `(lat, lon)` for `N < 40000`.
- `FS` low/upper-mid thresholds.
- `GM` strict-sphere thresholds and trial resolution choices for `N < 50000`.

The parallel 16-worker run initially showed apparent gains such as:

- `p16_combo_all`: `case5 1148 -> 1133`, `torus23scr 1907 -> 1282`, `wavy57 3127 -> 2257`.
- `p12_gm_mid22`: `case5 1148 -> 1133`, `wavy57 3127 -> 2268`.
- `p05_gi_26x52`: `low23_wavy 4003 -> 3795`.

However, an isolated rerun showed the parallel baseline itself was timing-contaminated:

| Variant | Isolated result |
| --- | --- |
| `base` | `low23_wavy 3795`, `case5 1133`, `torus23scr 1282`, `wavy57 2257`. |
| `p16_combo_all` | Identical to isolated base on tested proxies. |
| `p05_gi_26x52` | Identical to isolated base on tested proxies. |
| `p12_gm_mid22` | Same on most proxies but regressed `torus23scr 1282 -> 1907`. |

Decision: no submission. Treat the primitive-threshold same-size batch as a timing false-positive. Do not submit `p16`, `p12`, `p05`, or the related `GI/FS/GM` threshold-lowering variants without a new isolated local gain.
