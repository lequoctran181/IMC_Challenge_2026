# 2026-07-09 H2Q high2 postmortem

Context: user paused the 8 ChatGPT Pro Extended/browser workers and raised the local lane ceiling to 16. This run tested whether the source-slack `submission_1499_81.94_7.cpp` base could carry a stronger exact high2 post-pass without touching the source-size-limited `1448` high-water source.

## Local evidence

- Base `1499` is Kattis-valid at `81.938571`, `7/7`, and locally matches `1448` on `uvsphere_4098`, `closed_bunnylike_35292`, and high2 Nefertiti.
- Offline `trimesh` QEM decimation of the current high2 output showed real headroom:
  - current high2 output: about `150190/300376`
  - QEM 130k vertices: `vps_eval256 ~= 0.943839`
  - QEM 110k vertices: `vps_eval256 ~= 0.942482`
  - QEM 90002 vertices: `vps_eval256 ~= 0.941124`
  - QEM 70002 vertices: `vps_eval256 ~= 0.939775`
  - QEM 50002 vertices: `vps_eval256 ~= 0.935889`
- Inline H2Q post-pass on exact high2 reached about `98856/197708`, with:
  - `vps_eval256 ~= 0.931953`
  - `vps_eval512 ~= 0.903396`
  - isolated local high2 runtime about `17.82s`
- Smoke tests before submit:
  - sample: `8 12`
  - Lucy proxy: `15656/31306`
  - closed bunny-like test4 proxy: `1942/3880`

## Submitted candidate

- Local file: `submission_1585_h2q_high2_pending.cpp`
- Fetched Kattis source: `fetched_sources/kattis_19928674_h2q_53.91929_5.cpp`
- Size: `107998` bytes
- Kattis submission: `19928674`
- Result: `Accepted (53.91929)`, `5/7`, runtime `>21.00s`
- Detail page first hint: `Test case 4/7: Wrong Answer: SSIM is too low`

## Conclusion

The geometry idea is promising for high2, but the implementation route is not judge-safe. Adding a sizeable H2Q namespace/post-pass to the slack `1499` source perturbs the same source-layout/timing-sensitive region that repeatedly breaks hidden test4. It also leaves too little runtime margin on the largest case.

Blacklist:

- More raw `1499 + H2Q` exact-high2 submits.
- Source-layout-heavy post-passes unless test4 is independently protected.
- H2Q variants with target under about 100k vertices; local runtime is too close to the judge limit.

Keep:

- The offline evidence that high2 has real QEM headroom.
- The idea of a high2 visual/QEM branch, but only if implemented as a compact, runtime-safe, test4-stable route or a standalone early branch that returns before the fragile pipeline without affecting medium cases.
