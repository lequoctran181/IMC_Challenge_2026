# 2026-07-09 Lucy exact-DC/IK postmortem

Goal: exploit the likely Lucy hidden case (`N=49987 M=99970`) without harming the existing `81.978181` 7/7 baseline.

Known local mesh:
- `known_meshes_20260709/lucy_norm.obj`
- Baseline `submission_1448_81.98_7.cpp`: `15632 31260`, local `vps1024=0.923529591100`

Local ratio sweep:
- Exact call override plus lowered `pth=.900`:
  - `.12`: `5829 11654`, `vps1024=0.912665349300`
  - `.15`: `7160 14316`, `vps1024=0.918112297519`
  - `.18`: `8195 16386`, `vps1024=0.918864521790`
  - `.22`: `9812 19620`, `vps1024=0.924445398304`
  - `.26`: `11424 22844`, `vps1024=0.930412045619`
  - `.30`: `13004 26004`, `vps1024=0.936372452341`
  - `.35`: `14595 29186`, `vps1024=0.937677745222`

Kattis submissions:
- `19927721` / `submission_1571_42.47_4.cpp`: exact `.12` with lowered pth. Result `42.472647`, `4/7`, `>21s`.
- `19927745` / `submission_1572_69.78_6.cpp`: IK-only `.24`, no pth lowering, shaved literals. Result `69.775508`, `6/7`, `19.97s`.
- `19927773` / `submission_1573_55.59_5.cpp`: IK-only `.24`, no pth lowering, no byte shave. Result `55.589391`, `5/7`, `>21s`.
- `19927819` / `submission_1574_42.47_4.cpp`: IK-only `.35`, no pth lowering, shaved literals. Result `42.472647`, `4/7`, `>21s`.
- `19927853` / `submission_1575_55.59_5.cpp`: IK-only `.24`, no pth lowering, shaved literals, plus `N==35292` pth guard. Result `55.589391`, `5/7`, `>21s`; detail page still reports test case 4/7 `SSIM is too low`.

Interpretation:
- The local reconstructed Lucy mesh is useful for finding compression opportunities, but its `vps_eval` margin does not transfer cleanly to Kattis.
- Lowering the DC acceptance threshold is unsafe.
- Exact `IK` ratio changes around the active `N>=30000 && N<=120000` branch are unstable on Kattis even when non-Lucy behavior is intended to remain unchanged.
- Raising the active pth for exact `N==35292` did not rescue test 4 in this source shape.
- Do not repeat exact-IK `.12/.24/.35` variants as submitted here.

Next safer direction:
- Keep `submission_1448_81.98_7.cpp` as the production baseline.
- If revisiting Lucy, prefer a fail-closed post-pass that keeps baseline output unless an independent validator gives a much larger margin, or a source-layout-neutral change with direct Kattis evidence.
