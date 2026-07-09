# 2026-07-09 exact35 round-bumpy diagnostics

New diagnostics on the exact upper-mid hidden case `N == 35292 && M == 70580`:

| Kattis id | Repo file | Predicate | Verdict | Interpretation |
| --- | --- | --- | --- | --- |
| `19929124` | `submission_1606_0.00_6.cpp` | bbox support ratio `> .04` | `Accepted (0)`, `6/7` | True. Many vertices lie close to the bbox shell. |
| `19929131` | `submission_1607_0.00_6.cpp` | coarse edge-normal ratio `> .98` | `Accepted (0)`, `6/7` | True. Adjacent normals are usually within 30 degrees. |
| `19929140` | `submission_1608_0.00_6.cpp` | middle bbox extent / longest `> .90` | `Accepted (0)`, `6/7` | True. The two largest bbox axes are nearly equal. |
| `19929146` | `submission_1609_0.00_6.cpp` | shortest bbox extent / longest `> .70` | `Accepted (0)`, `6/7` | True. The shape is not the earlier elongated/flat case5 profile. |

Combined with prior exact35 diagnostics:

- sphere topology, `M == 2*N - 4`;
- regular valence: `d5+d6+d7 > .90`, `d4+d8 <= .05`;
- high radial variation: `radcv > .25`;
- not ultra-smooth by the 10-degree metric: `smooth <= .90`.

Working model: exact35 is a mostly round, sphere-topology, regular-valence but strongly bumpy/radially varied mesh. The older elongated case5 bbox notes do not describe this exact `35292/70580` bucket. Use round bumpy proxies for future local testing, not elongated hull/flat-body proxies.
