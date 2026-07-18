# Kattis submission 20082128

- Judgement: Accepted
- Tests: 7/7
- Score: 93.812395
- Output vertex counts: `[25, 4340, 2839, 3030, 7800, 16500]`
- Authoritative source bytes: 130434
- SHA-256: `9696df574ae41e7df56a9fd975884ca814de72014b43e2c14ad5c5528055c038`

This isolated Nefertiti change replaces the 16,650-vertex output from
submission `20081873` with a 16,500-vertex QEM continuation using five
pre-passes and two late passes. The five non-Nefertiti outputs are
byte-identical to `20081873`. The count-only score is `93.812394698355`,
matching Kattis.

The exact source was fetched back after the 7/7 result. The Nefertiti output
is VALID and its local judge-separating metric exceeds the previously accepted
16,650-vertex mesh.
