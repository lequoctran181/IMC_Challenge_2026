# 2026-07-09 Armadillo exact-23 RP-gate postmortem

Context: after `19928275` showed that unconditional `RP()` cap-pair repair was too aggressive on hidden test 4, a gated variant was tested. The branch kept the `RP()` repair but snapshotted the state and restored it when internal `vps(512) < .97`.

Local proxy behavior for `rpgate_97`:

- `agg3`: rollback to `23201 46398`
- `agg5`: retained the base compressed family `11507 23010`
- `agg10`: rollback to `23201 46398`
- `agg20`: rollback to `23201 46398`
- `closed_bunnylike_35292`: unchanged `1942 3880`
- `Lucy`: unchanged first line `15632 31260`

Kattis result:

- Submission: `19928418`
- Local file: `submission_1581_43.09_4.cpp`
- Judgement: `Accepted (43.094425)`
- Tests: `4/7`
- Runtime: `> 21.00 s`
- Feedback: secret/2 and secret/3 failed with `SSIM is too low`

Decision: blacklist `RP()` plus internal `vps` rollback. Even a high `.97` proxy gate changes timing/source behavior enough to fail multiple hidden tests. Future exact-23 attempts must avoid post-final render gates and should either be byte-neutral inside an already judge-stable path or early-return only with an independently valid mesh.
