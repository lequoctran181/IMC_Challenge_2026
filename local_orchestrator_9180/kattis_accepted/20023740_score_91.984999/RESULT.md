# Kattis submission 20023740

- Judgement: Accepted, 7/7
- Score: 91.984999
- Source SHA-256: `e9cc1208576f0534f4c5d07eb0f851f5a289ae0f11f3d8792c35bba748745d3d`
- Hidden output vertex counts: `[25, 6116, 3800, 3150, 8296, 18669]`

This submission replaces the proven Bunny3850 tail with a renderer-aware
two-stage topology search, 82 independent edge flips, and a compact
three-bit coordinate correction.  The isolated Bunny reduction from 3850 to
3800 vertices passed all seven Kattis tests.
