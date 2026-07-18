# Kattis submission 20079731

- Judgement: Accepted
- Tests: 7/7
- Score: 93.796552
- Output vertex counts: `[25, 4340, 2848, 3040, 7800, 17000]`
- Authoritative source bytes: 130397
- SHA-256: `bfbc9d57715e30bb592af3a1f6d5b2cb07fa03d97739144ca855e7e6d72e0a71`

This isolated improvement lowers the accepted Nefertiti output from 17100 to
17000 vertices while preserving the other five outputs. The count-only score
is `93.796552207`, matching the rounded Kattis result.

Local Nefertiti QA before submission produced 17000 vertices and 33996 faces,
passed manifold/Hausdorff validation, and scored `0.903253296700` with the
exact VPS evaluator at resolution 1024. The 7/7 status and score were verified
with `kattis_manager.py list`; the archived source was fetched back from Kattis.
