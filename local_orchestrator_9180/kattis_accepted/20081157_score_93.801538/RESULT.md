# Kattis submission 20081157

- Judgement: Accepted
- Tests: 7/7
- Score: 93.801538
- Output vertex counts: `[25, 4340, 2848, 3030, 7800, 16900]`
- Authoritative source bytes: 130469
- SHA-256: `c4ae7bfa77d710e6177edd541ad44dce9f44c326f82f911d7370b2fc3731083c`

This isolated improvement rebases the judge-proven Lucy topology replay from
3040 to 3030 vertices onto submission `20079783`. The other five judged
outputs are unchanged. The count-only score is `93.801538014`, matching the
rounded Kattis result.

Local QA produced a VALID 3030-vertex Lucy mesh with exact rotation-00
VPS1024 `0.909291696953` and Hausdorff/tolerance ratio `0.935986134`.
The source passed the official sample, compiled with GCC 16, and was fetched
back from Kattis after the 7/7 result.
