# Results, milestones, and evidence

## Official submission result and provisional standings context

<!-- BEGIN GENERATED: RELEASE_SUMMARY -->
Kattis submission **20082703** was **Accepted (7/7)** with displayed score **93.830074**; the count-derived value is **93.83007422510956**. It reduces 1,498,780 input vertices to 34,134 (2.277452% retained; 97.722548% global compression). A Kattis snapshot captured at 2026-07-19T12:35:16.490Z displayed the team at rank 2, and the same page explicitly labelled the standings *Unfinalized*. Both fields are recorded verbatim and kept separate from the submission result.

| Case | Input vertices | Output vertices | Retained | Compression |
|---|---:|---:|---:|---:|
| Sphere-like sample | 4,098 | 25 | 0.610054% | 99.389946% |
| Armadillo | 23,201 | 4,340 | 18.706090% | 81.293910% |
| Bunny-like | 35,292 | 2,839 | 8.044316% | 91.955684% |
| Lucy | 49,987 | 3,030 | 6.061576% | 93.938424% |
| Slender | 377,084 | 7,400 | 1.962427% | 98.037573% |
| Nefertiti | 1,009,118 | 16,500 | 1.635091% | 98.364909% |
| **Global count aggregate** | **1,498,780** | **34,134** | **2.277452%** | **97.722548%** |
<!-- END GENERATED: RELEASE_SUMMARY -->

- Source size: **130,973 bytes**
- Source SHA-256: `9195d42a73a6b85c8ae30d731f532175bdcd7c2982d421143d631b4c64b1a92c`

The byte-exact source and Kattis fetch-back record are archived in [`submission/`](../submission/).

## Per-case score reconstruction

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

## Controlled comparisons, ablation evidence, and negative results

The competition process contained many intertwined edits, so only interventions with an otherwise stable comparison are treated as clean ablations.

1. **Cluster-normal memory.** In the same-parent 5,471-vertex Bunny public-proxy experiment, carrying additive original vertex-face-incidence normal evidence increased sixteen-rotation mean combined SSIM from 0.90625420 to 0.90899736 and reduced final support drift from 5.8627 to 5.0502 degrees. This is strong Experimental evidence for that proxy; hidden transfer remains Inference.
2. **Reference cache on Slender.** Submission 20082666 in the 7,400 lineage completed 6/7, whereas cache-equipped submission 20082703 passed 7/7. Runtime causality remains Inference because the curated artifact does not contain a canonical geometry-equivalence audit for the failed source.
3. **Structural replay.** The final Slender step reduced 400 vertices while preserving the other five output counts byte-for-byte relative to the selected lineage, isolating a score gain of exactly 0.017679526754 points under the published formula.

Important approaches that did **not** survive the validation funnel:

- isotropic remeshing: good depth silhouettes, poor flat-normal SSIM;
- shared-vertex multi-layer triangulations: valid topology but unstable z-buffer ordering;
- third-party generic decimators: weaker under the specification-matching local renderer than the tuned guarded QEM;
- lowering target counts without structural search: frequent hidden perceptual failures;
- proxy-only parameter wins: often failed a different rotation or the official hidden mesh.

Negative results are preserved in the worker and orchestration directories instead of being silently removed. They document why the final method is hybrid.

## Evidence policy

Claims in the article and this repository are assigned one or more of four evidence levels:

- **Official:** a Kattis judgement, displayed score, test count, or fetched-back source;
- **Reconstructed:** deterministic arithmetic, byte count, checksum, or count-derived score;
- **Experimental:** a controlled measurement on a legally obtained public proxy;
- **Inference:** an interpretation not directly exposed by Kattis.

No hidden mesh is included or reverse-engineered into this repository.
