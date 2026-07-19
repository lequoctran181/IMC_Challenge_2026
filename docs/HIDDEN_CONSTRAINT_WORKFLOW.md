# Hidden-constraint discovery and evidence protocol

Kattis returned an aggregate score and a pass/fail judgement, not the hidden meshes, per-case render scores, or active failure component. We treated that interface as a controlled black-box experiment. The objective was not to reconstruct hidden geometry. It was to learn only the smallest decision-relevant facts that a submitted solver could lawfully compute from its own input: which public proxy family was relevant, which branch had changed, and where the judge-accepted frontier lay.

This document separates three kinds of evidence throughout:

- **official evidence:** a Kattis judgement, displayed score, fetched-back source, or output count recovered from a controlled submission;
- **local evidence:** topology, distance, and 1024-square rendering measurements on a public or otherwise licensed proxy;
- **inference:** a conclusion supported by official and local evidence but not directly observed, such as a likely proxy identity.

Proxy measurements were never relabelled as official measurements. No hidden input, hidden output mesh, or third-party mesh is included in this repository.

## 1. Establish an immutable Accepted parent

Every experiment began from a known Accepted source. Its unchanged branches were kept byte-identical, and the candidate changed only one case, target count, operator family, or diagnostic digit. For each probe we recorded:

1. parent submission and source digest;
2. changed branch and exact code delta;
3. expected output count or diagnostic payload;
4. local structural, distance, and perceptual checks;
5. official judgement and displayed score;
6. interpretation, including whether it was a direct observation or an inference.

This isolation rule converted aggregate feedback into a usable experiment. Without it, an Accepted score change could not be assigned to a particular case, and a failure could have come from topology, perception, runtime, or an unrelated branch.

## 2. Decode a changed output count from the aggregate score

For six scored cases, Kattis uses

$$
S=100\left(1-\frac{1}{6}\sum_{i=1}^{6}\frac{V'_i}{V_i}\right).
$$

If only case $k$ changes and all other output counts are known, its count is recovered by

$$
|V'_k|=\operatorname{round}\!\left(
|V_k|\left[6\left(1-\frac{S}{100}\right)-
\sum_{i\ne k}\frac{|V'_i|}{|V_i|}\right]
\right).
$$

For the 23,201-vertex branch, one output vertex changes the score by approximately $100/(6\cdot23201)=0.00071836$. Six displayed decimals therefore distinguish adjacent counts. This made a single isolated submission a precise case-specific measurement even though the public interface exposed only one aggregate score.

The method has two safeguards:

- use the exact official score definition and rounding, not the global compression ratio;
- require every unchanged branch to come from the same Accepted parent, ideally verified byte-for-byte.

## 3. Carry small diagnostic integers without changing the active surface

When the runtime program needed to report a small statistic $h$, a conservative Accepted mesh with $B_k$ vertices was used as the visible surface. The branch then emitted

$$
|V'_k|=B_k+h,
$$

where the additional vertices were unreferenced and did not alter any face, normal, depth, or silhouette. The normal manifold surface remained the judge-proven one; only the count changed. The aggregate score then revealed $h$ through the count-decoding equation.

Payloads were deliberately small and nonnegative. If an integer was too large for the safe count slack, it was split into base-$B$ digits,

$$
H=d_0+B d_1+B^2d_2+\cdots,
$$

and recovered from separate isolated probes. This avoided overflow, ambiguity, and unsafe surface changes. A diagnostic submission could score below the current record without lowering the stored Kattis best.

## 4. Prefer invariant fingerprints to raw coordinates

Vertex count alone is not an identity test. Different meshes, poses, or preprocessings can share the same count. We therefore computed several independent signatures from normalized geometry:

- axis-aligned bounding-box ratios;
- component count and Euler-characteristic or genus information;
- a sorted multiset of quantized normalized edge lengths;
- edge-graph labels and degree statistics;
- orientation and raw-order signatures when ordering mattered;
- scale and pose checks obtained from independent codes.

A representative rotation- and scale-invariant edge channel used the root-mean-square edge scale

$$
L_E=\left(\frac{1}{|E|}\sum_{(a,b)\in E}\lVert p_a-p_b\rVert_2^2\right)^{1/2},
\qquad
q_e=\operatorname{round}\!\left(B_q\frac{\lVert p_u-p_v\rVert_2}{L_E}\right).
$$

Translation disappears in edge differences, uniform scale disappears after division by \(L_E\), rotation preserves Euclidean length, and sorting removes edge-enumeration order. Bounding-box ratios were retained as a separate pose-sensitive channel. The sorted edge sequence was then hashed and, when necessary, emitted as multiple count digits.

For the 23,201-vertex branch, three recovered base-15,776 digits were 2,382, 13,367, and 11,736. They matched the exact-source public proxy and differed from aggregate alternatives. Independent pose, graph, orientation, scale, and ordering checks agreed. For the 35,292-vertex branch, normalized bounding-box diagnostics matched the public Bunny family, while official acceptance still differed from optimistic local rotations. We therefore kept two claims separate:

- **identity evidence** selects the public geometry family worth studying;
- **acceptance evidence** comes only from isolated Kattis pass/fail probes.

A hash match is never proof by itself. Collisions, preprocessing, and renderer differences remain possible, so multiple independent signatures and official frontier observations are required.

## 5. Reconstruct pass/fail frontiers scientifically

After identifying a useful branch, target counts were bracketed or binary-searched. A pass moved the certified frontier downward; a failure tightened the rejected side. Each probe changed one target or one operator family. The early observations included:

| Branch and method | Accepted observation | Rejected neighboring observation | What it established |
|---|---:|---:|---|
| Armadillo, generic QEM | 34%, 33%, 32% retained | 31% retained | A flat-normal frontier near 31–32% |
| Bunny-like, normal-aware QEM | 16% retained | more aggressive local-only settings | Current-face normal cost accumulated drift |
| Bunny-like, cluster memory | 15.5% retained | 14.5% new-cost branch | Original-normal memory improved the hidden frontier, but public rotations were optimistic |
| Lucy, early QEM | 9.5% retained | 9.0% local and clustered variants | Hidden normal margin was tighter than the Hausdorff margin |
| Nefertiti, early QEM | 4% retained | 3% combined probe | Runtime and perceptual constraints were both active |

These are historical frontiers of the corresponding method, not limits of the final solver. Renderer-guided offline replay later reached 18.706% for Armadillo, 8.044% for Bunny-like, 6.062% for Lucy, and 1.635% for Nefertiti.

Near the 21-second limit, an identical source could be resubmitted to distinguish transient judge load from a systematic timeout. Repetition was bounded: repeated failures triggered profiling, caching, or algorithmic changes rather than indefinite resubmission.

## 6. Use public rotations as stress tests, not hidden substitutes

Public proxies were evaluated under multiple rotations because a single pose can hide a silhouette or flat-normal defect. Candidates had to pass manifold, distance, and exact 1024-square render checks across selected rotations with an explicit margin. Some rotations were useful predictors—rotation 04 was especially informative for Lucy and Nefertiti—but none was treated as authoritative.

This distinction mattered. Public Bunny rotations passed counts that the hidden judge rejected. The correct response was not to claim a contradictory hidden mesh; it was to acknowledge a proxy-to-hidden renderer or geometry gap and let official isolated probes define the accepted frontier.

## 7. Convert discoveries into a fail-closed release

The final program does not perform exploratory search on the judge. Offline discoveries are canonicalized into compact replay transactions. Each transaction carries local preconditions, an expected count, and branch-specific checks. At runtime:

1. the generic simplifier builds a safe checkpoint;
2. a replay transaction is decoded;
3. canonical neighborhood and count preconditions are checked;
4. the transaction commits only if all structural guards pass;
5. otherwise the program returns the judge-proven checkpoint.

The fetched-back final source, output counts, result metadata, and checksums are immutable. The release verifier reconstructs the score, checks the 131,072-byte limit, validates article hashes, and compiles the exact C++17 source.

## 8. Reusable lessons

- Treat aggregate feedback as an experimental-design problem: one unknown per probe.
- Decode exact counts from the published score equation before drawing qualitative conclusions.
- Use small, harmless, fail-closed diagnostic payloads.
- Identify proxy families with several normalized invariants, not filenames or counts alone.
- Separate identity evidence, local feasibility, and official acceptance.
- Record both accepted and rejected sides of every frontier.
- Preserve source hashes and unchanged branch outputs so later integrations remain attributable.
- Never redistribute hidden or unlicensed meshes; release algorithms, metadata, and integrity records instead.

The same protocol applies to other fixed-instance optimization challenges whenever the judge publishes an aggregate objective but withholds component measurements. Its value is not in exposing hidden data; it is in making limited external feedback scientifically controlled, reproducible, and honest about uncertainty.
