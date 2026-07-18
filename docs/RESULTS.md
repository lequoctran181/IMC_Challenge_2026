# Results, milestones, and evidence

## Official final result

- Team: **NEU.AddictedTribes**
- Rank: **2nd**
- Submission: **20082703**
- Judgement: **Accepted, 7/7**
- Official score: **93.830074**
- Exact reconstructed score: **93.83007422510956**
- Source size: **130,973 bytes**
- Source SHA-256: `9195d42a73a6b85c8ae30d731f532175bdcd7c2982d421143d631b4c64b1a92c`

The byte-exact source and Kattis fetch-back record are archived in [`submission/`](../submission/).

## Per-case score reconstruction

| Case | $V$ | $V'$ | $V'/V$ | Compression | Score contribution |
|---|---:|---:|---:|---:|---:|
| Sphere-like sample | 4,098 | 25 | 0.610054% | 99.389946% | 16.564991 |
| Armadillo | 23,201 | 4,340 | 18.706090% | 81.293910% | 13.548985 |
| Bunny-like | 35,292 | 2,839 | 8.044316% | 91.955684% | 15.325947 |
| Lucy | 49,987 | 3,030 | 6.061576% | 93.938424% | 15.656404 |
| Slender | 377,084 | 7,400 | 1.962427% | 98.037573% | 16.339595 |
| Nefertiti | 1,009,118 | 16,500 | 1.635091% | 98.364909% | 16.394151 |
| **Sum** | **1,500,780** | **34,134** | — | — | **93.830074** |

The sample-like 4,098-vertex case is included because it is one of the six score-bearing counts exposed by the final run; the separate format-debugging sample is not part of the official mean.

![Final retained-vertex ratios and score contributions](../paper/figures/final_results.png)

## Selected accepted milestones

| Submission | Score | Main change |
|---:|---:|---|
| 19928674 | 81.945906 | First robust seven-test QEM baseline |
| 19930642 | 86.998654 | Case-aware targets and stronger validity guards |
| 19932621 | 87.913148 | Cluster-normal accumulation; Bunny-like 15.5% |
| 20021691 | 91.801173 | Compressed replay and large-case structural branch |
| 20025898 | 92.932731 | Armadillo/Bunny/Lucy/Nefertiti checkpoint advances |
| 20039214 | 93.600980 | Renderer-aware topology replay |
| 20051927 | 93.769981 | Armadillo and Bunny refinements |
| 20082128 | 93.812395 | Nefertiti reduced to 16,500 vertices |
| 20082703 | **93.830074** | Slender reduced from 7,800 to 7,400 with flips and cache |

![Selected accepted submissions over the optimization campaign](../paper/figures/score_progression.png)

## Clean ablations and negative results

The competition process contained many intertwined edits, so only interventions with an otherwise stable comparison are treated as clean ablations.

1. **Cluster-normal memory.** Carrying additive original-face normal evidence enabled the hidden Bunny-like case at a smaller retained ratio than the preceding local-normal-only branch. This exposed accumulated-reference drift as a real failure mode.
2. **Reference cache on Slender.** The 7,400-vertex geometry could fail the runtime test without the cache; the otherwise identical cached release passed 7/7. The final improvement therefore required systems engineering as well as geometry.
3. **Structural replay.** The final Slender step reduced 400 vertices while preserving the other five output counts byte-for-byte relative to the selected lineage, isolating a score gain of approximately 0.01696.

Important approaches that did **not** survive the validation funnel:

- isotropic remeshing: good depth silhouettes, poor flat-normal SSIM;
- shared-vertex multi-layer triangulations: valid topology but unstable z-buffer ordering;
- third-party generic decimators: weaker under the exact renderer than the tuned guarded QEM;
- lowering target counts without structural search: frequent hidden perceptual failures;
- proxy-only parameter wins: often failed a different rotation or the official hidden mesh.

Negative results are preserved in the worker and orchestration directories instead of being silently removed. They document why the final method is hybrid.

## Evidence policy

Claims in the article and this repository are assigned one of three evidence levels:

- **Official:** Kattis judgement, tests passed, score, and fetched-back source.
- **Reconstructed:** deterministic arithmetic from official output counts or checksum/byte-size measurements.
- **Experimental:** local evaluator results on legally obtained public proxies; these are not presented as hidden-test measurements.

No hidden mesh is included or reverse-engineered into this repository.
