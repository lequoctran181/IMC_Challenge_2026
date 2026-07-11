# Kattis result 20024979

- Problem: `simplifygeometry`
- Judgement: Accepted
- Tests: 7/7
- Score: 92.049968
- Source bytes: 131020
- SHA-256: `b57c32ca511b80c89d42fe208d19b238de806e7beb0b16a959719a021532df9f`
- Official-score-consistent output counts: `[25, 6106, 3750, 3100, 7900, 18669]`
- Improvement over 20024892: Lucy replay/coordinate/flip tail from 3150 to 3100 vertices.

The source is rebased on submission 20024530, whose Arm branch has 6106
vertices.  The score delta from 20024530 is exactly the 50-vertex Lucy gain;
the earlier 6104 proxy count was therefore corrected to 6106.

Exact GCC 14.2 x86 QA for Lucy rotation 00 produced 3100 vertices / 6196
faces, VPS1024 `0.909001323108`, and a VALID mesh with Hausdorff ratio
`0.935986134`.  The GCC output SHA-256 was
`b4eed83be3a0e490209df9d6d065f1606b91d9273418202926a388a13c07832d`.
Lima/QEMU timing was 157.70 s wall / 82.47 s user / 73.43 s sys and is not
representative of native Kattis; the native Lucy regression used about 7.84 s
user under contention.  `kattis_manager.py` did not expose official runtime.

The score and 7/7 status above were rechecked with `kattis_manager.py list` before archival. The stored `solution.cpp` was fetched directly with `kattis_manager.py source 20024979`.
