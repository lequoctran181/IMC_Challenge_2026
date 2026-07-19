# Certified Perceptual Mesh Simplification under a 21-Second and 128-KiB Budget

*Hybrid QEM, cluster-normal memory, renderer-aware topology replay, and fail-closed black-box optimization*

Problem name: Perception-Aware Lossless Simplification of Million-Vertex 3D Meshes for Mobile Platforms

Team name: NEU.AddictedTribes

GitHub repository: [Public research artifact](https://github.com/lequoctran181/IMC_Challenge_2026)

<!-- PAGE BREAK -->

# Certified Perceptual Mesh Simplification under a 21-Second and 128-KiB Budget

*Hybrid QEM, cluster-normal memory, renderer-aware topology replay, and fail-closed black-box optimization*

## Abstract

The IMC Challenge problem `simplifygeometry` asks for the smallest closed triangular meshes that remain visually indistinguishable from six high-resolution inputs under a fixed six-view renderer. A valid output must simultaneously be a non-degenerate watertight two-manifold, stay within a symmetric vertex-Hausdorff tolerance equal to five percent of the input bounding-box diagonal, and obtain at least 0.9 structural similarity on equally weighted flat-normal and depth maps [1]. The program must additionally run within 21 seconds, use at most 2 GiB of memory, and fit in 131,072 source bytes. These coupled requirements make the feasible set discontinuous: one contraction can change topology, violate the distance bound, or alter visibility and thousands of rasterized pixels.

We present a hybrid certified simplification system. A guarded quadric-error-metric backbone supplies scalable geometric contraction [2,3]. Its local objective is augmented by projected-area normal distortion and cluster-normal memory, an additive statistic that carries area-weighted normals of original support faces through the full collapse hierarchy. This prevents reference drift caused by comparing a candidate only with an already degraded current mesh. An exact 1024-by-1024 evaluator reproduces the prescribed cameras, z-buffer, flat normals, perspective depth, foreground mask, and SSIM calculation [1,10]. It guides offline edge flips, one-ring fan deletions, and bounded coordinate fitting. Successful operations are canonicalized, compressed, and replayed transactionally at runtime; any precondition mismatch or count error falls back to a judge-proven checkpoint.

Because only aggregate Kattis feedback was observable, we used a controlled black-box experimental design. Each probe modified one branch while preserving the other outputs byte-for-byte. Rotation-invariant geometric fingerprints were encoded in safe output-count slack and recovered algebraically from the displayed score, enabling proxy identification without obtaining or redistributing hidden meshes. Binary searches then mapped hidden acceptance frontiers; multi-rotation proxy tests and conservative margins controlled overfitting. Every accepted high-water source was fetched back, hashed, and archived before integration.

The final submission, Kattis 20082703, reduces the six scored inputs from 1,500,780 to 34,134 vertices: 25, 4,340, 2,839, 3,030, 7,400, and 16,500. It passed all seven tests and achieved 93.830074 points, placing NEU.AddictedTribes second. The exact count reconstruction is 93.83007422510956. The result shows that aggressive perceptual simplification under hard constraints is best treated as a layered research problem: conservative geometry establishes feasibility, renderer-aware search recovers appearance, controlled probes reduce hidden uncertainty, and deterministic replay converts expensive offline discovery into a legal contest program.

### Keywords

mesh simplification; quadric error metrics; structural similarity; renderer-aware optimization; edge collapse; topology preservation; Hausdorff distance; normal maps; black-box experimental design; deterministic replay; geometry compression

## 1. Introduction

### 1.1 Motivation and challenge context

High-resolution scans and authored models frequently contain far more triangles than an interactive or mobile renderer can process economically. Classical simplification replaces a dense surface with a coarser approximation, commonly using geometric distance, local curvature, or quadric error as a surrogate for quality [2-6]. The IMC Challenge changes the problem in two decisive ways. First, acceptance is renderer-defined: six fixed views are compared through flat-normal and perspective-depth images using SSIM. Second, visual similarity is necessary but not sufficient; the output must remain a closed, non-degenerate two-manifold and satisfy a hard symmetric vertex-set Hausdorff bound [1]. The target is therefore neither ordinary decimation nor unconstrained image matching.

The renderer couples combinatorics and geometry at pixel scale. A small coordinate perturbation can change depth ordering. An edge flip can alter a large constant-normal region without changing the vertex count. A geometrically inexpensive collapse can erase a silhouette feature that dominates foreground-only SSIM. Conversely, vertices invisible from all six prescribed cameras can contribute little to the perceptual score while remaining essential for closure or the distance certificate. The solver must reason simultaneously about surface geometry, topology, visibility, numerical error, runtime, and source-code entropy.

Our development began with a generic simplifier and evolved into a hybrid research pipeline. Online QEM provides a scalable proposal mechanism. Topological and radius guards establish a safe backbone. Projected normal and cluster-memory terms better align that backbone with the flat-shaded evaluator. Exact local rendering then searches structural operations that QEM alone cannot express. Finally, replay compilation moves the cost of search offline. This separation is essential: rendering every candidate collapse at 1024 squared is computationally impossible for million-vertex inputs, while a local geometric metric alone plateaued far below the final score.

### 1.2 Formal task and objective

For an original triangular mesh \(M=(V,F)\), let \(M'=(V',F')\) denote the output. The official statement fixes the structural validity predicate, a symmetric vertex-set Hausdorff tolerance, six cameras, flat-normal and depth maps, a foreground-only SSIM calculation, and a threshold of 0.9 [1]. For view \(c\), we write

$$s_c=\frac{1}{2}\operatorname{SSIM}\!\left(I^{N}_{c},I^{N'}_{c}\right)+\frac{1}{2}\operatorname{SSIM}\!\left(I^{D}_{c},I^{D'}_{c}\right).$$

The aggregate perceptual feasibility condition is

$$S_{\mathrm{vis}}=\frac{1}{6}\sum_{c=1}^{6}s_c\geq 0.9.$$

Let \(D_{\mathrm{AABB}}\) be the diagonal of the original axis-aligned bounding box, and define the directed distance \(\delta(A,B)=\max_{x\in A}\min_{y\in B}\lVert x-y\rVert_2\). The geometric constraint is

$$d_H(V,V')=\max\!\left\{\delta(V,V'),\delta(V',V)\right\}\leq\tau_H,\qquad\tau_H=0.05D_{\mathrm{AABB}}.$$

For each case, the constrained problem can be expressed as

$$\min_{M'}|V'|\quad\text{subject to}\quad\mathcal{T}(M')=1,\quad d_H\leq\tau_H,\quad S_{\mathrm{vis}}\geq0.9.$$

If all six scored cases are valid, Kattis ranks a submission using the unweighted mean of the six retained-vertex ratios [1]:

$$\operatorname{Score}=100\left(1-\frac{1}{6}\sum_{i=1}^{6}\rho_i\right),\qquad \rho_i=\frac{|V'_i|}{|V_i|}.$$

This formula has an important allocation consequence. Removing one vertex from case \(i\) is worth \(100/(6|V_i|)\) points. Thus a vertex removed from the 23,201-vertex case is about 43.5 times as valuable as one removed from the 1,009,118-vertex case. A uniform target ratio is therefore rarely optimal: smaller difficult models deserve disproportionate search effort, while the largest models still require aggressive reduction to meet runtime and memory limits.

### 1.3 Assumptions, symbols, and evidence policy

The following assumptions are explicit. The official cameras, focal length, resolution, background values, flat-normal rasterization, perspective-correct depth interpolation, foreground mask, and SSIM definition are treated exactly as stated [1]. The official Hausdorff metric is between vertex sets, not continuous surfaces; we nevertheless use conservative incremental certificates and an independent final checker. The six scored instances are fixed, and case specialization is permitted, but every specialized branch recognizes its intended input deterministically and fails closed. Kattis acceptance is final external validation, not a substitute for local testing.

| Symbol | Meaning |
|---|---|
| \(M=(V,F)\), \(M'=(V',F')\) | Original and simplified triangular meshes |
| \(p_v\) | Position of vertex \(v\) |
| \(n_f\), \(A_f\) | Unit normal and area of face \(f\) |
| \(Q_v\) | Symmetric homogeneous quadric at vertex \(v\) |
| \(C(v)\), \(r_v\) | Original support cluster represented by \(v\), and its radius certificate |
| \(D_{\mathrm{AABB}}\) | Original bounding-box diagonal |
| \(I^N_c\), \(I^D_c\) | Normal and depth image for view \(c\) |
| \(S_v\) | Additive cluster-normal summary at vertex \(v\) |
| \(S_{\mathrm{vis}}\) | Mean six-view normal/depth SSIM |
| \(\rho_i\) | Retained-vertex ratio of scored case \(i\) |

Claims are separated into three evidence levels. **Official** means Kattis judgement, tests passed, displayed score, or the fetched-back submitted source. **Reconstructed** means deterministic arithmetic, checksum, byte count, or output count derived from archived official evidence. **Experimental** means a local measurement on a legally obtained public proxy or a controlled diagnostic; it is never presented as direct access to a hidden mesh. This distinction is central to the Results and Discussions section.

### 1.4 Research questions and contributions

The work addresses four research questions. How can a million-vertex mesh be simplified quickly while preserving closed-manifold structure? How can a cheap local metric better predict a flat-normal renderer? How can exact image evidence be exploited without evaluating every collapse? How can aggregate judge feedback be converted into controlled information about hidden constraints without contaminating other cases?

The main contributions are:

- A manifold-preserving QEM implementation with link-condition checks, duplicate-face rejection, orientation and degeneracy guards, a conservative cluster-radius certificate, and versioned lazy heaps suitable for inputs above one million vertices.
- A projected-area normal term and cluster-normal memory that preserves additive orientation evidence from original support faces instead of repeatedly comparing only against the current mesh.
- A renderer-aware structural stage combining exact six-view evaluation, edge flips, one-ring retriangulation, coordinate fitting, and multi-rotation robustness tests.
- A black-box inference protocol based on isolated probes, algebraic count decoding, safe diagnostic payloads, rotation-invariant fingerprints, binary frontier search, and explicit proxy-to-hidden uncertainty margins.
- A fail-closed replay architecture that packs QEM, structural operators, case streams, validation checks, and large-case caches into 130,973 C++ bytes while staying within the 21-second limit.
- A reproducible evidence chain of accepted milestones culminating in an exact score reconstruction and a public research artifact [14].

## 2. Related Literature

### 2.1 Geometry-driven simplification

Progressive meshes established edge collapse and its inverse vertex split as a compact multiresolution representation, with an optimization objective that can account for geometry and appearance attributes [4]. Garland and Heckbert then introduced QEM, representing accumulated squared distances to planes by a symmetric four-dimensional matrix and ordering local contractions by a quadratic form [2]. Its compact ten-coefficient representation and additive update are especially well suited to a strict runtime and source budget. Garland and Heckbert later extended the formulation to color and texture attributes [3], while Hoppe developed a correspondence-based quadric for normals, colors, and other appearance data [6].

Recent high-quality work confirms that QEM remains an active framework rather than merely a historical baseline. Liu et al. replace extrinsic quadrics with accumulated intrinsic tangent information to support coarse operators and guarantee intrinsic element quality in ACM Transactions on Graphics 2023 [12]. Liu, Rahimzadeh, and Zordan add line quadrics to control distribution, feature preservation, and numerical conditioning in Computer Graphics Forum 2025 [13]. These studies reinforce two principles used here: useful global information can be agglomerated through local collapses, and additional quadratic or cluster statistics can control properties not captured by point-to-plane error alone. Our objective differs because the authoritative endpoint is a fixed flat-normal/depth renderer rather than an intrinsic differential operator or a uniform sampling objective.

### 2.2 Topology and geometric guarantees

Local edge contraction is attractive because only a one-ring neighborhood changes, but unconstrained contraction can invert faces, introduce duplicate triangles, change genus, or create non-manifold edges. Dey et al. formalized local conditions under which edge contractions preserve topological type [9]. We use the interior link condition as a fast necessary structural guard, supplemented by explicit face and edge-incidence tests. Because the contest demands a closed two-manifold, we intentionally do not use QEM's more general non-edge aggregation mode.

Hard geometric error bounds have a separate literature. Simplification envelopes constrain the evolving surface between offsets of the original and thereby provide a global guarantee while preserving topology [7]. Metro evaluates geometric error between surfaces through sampling [8]. The challenge instead defines a symmetric distance between vertex sets [1], so our certificate follows original vertices through representative clusters and the final checker evaluates the exact stated metric. We cite envelope and surface-sampling work to clarify the distinction: the contest certificate is task-aligned, but it is not a universal continuous-surface guarantee.

### 2.3 Appearance- and image-driven simplification

Pure geometric proximity cannot fully predict rendered appearance. Lindstrom and Turk's image-driven simplification evaluates edge collapses through images, automatically balancing silhouette, shading, texture, and hidden-region effects [5]. More recently, Hasselgren et al. formulated appearance-driven model simplification as analysis by synthesis, jointly optimizing geometry and shading through rendered image differences [11]. These works are the closest conceptual relatives of our renderer-aware stage.

The contest setting introduces unusual constraints absent from a general offline renderer. The images are fixed flat normals and depths; the acceptance threshold is hard; only six views matter; hidden instances must be solved by a 21-second C++ program; and the source limit leaves only 99 bytes of final headroom. We therefore do not make the full renderer differentiable or place it inside every collapse. Instead, a cheap surrogate guides bulk decimation, exact rendering is reserved for checkpoints and small structural neighborhoods, and accepted operations are compiled into replay streams.

### 2.4 Structural similarity and task-specific perception

SSIM compares local luminance, contrast, and structural terms rather than treating pixels as independent squared errors [10]. The official problem applies a specified SSIM variant to foreground normal and depth images [1]. This matters algorithmically. A few changed silhouette pixels, a normal discontinuity spanning a large triangle, or a z-buffer ownership change can affect local windows in ways that are not proportional to Euclidean vertex displacement. The optimization landscape is therefore non-smooth even when the underlying coordinates vary continuously.

Our local projected-area term should not be interpreted as a new universal perceptual metric. It is a computational surrogate for pixel coverage under the six official cameras. The exact evaluator remains authoritative. The contribution is the division of labor between a fast collapse-scale approximation and an exact checkpoint-scale decision rule.

### 2.5 Positioning of the present work

The method combines ideas from geometry-driven, topology-preserving, and image-driven simplification but adds a contest-specific systems layer. QEM supplies throughput [2]. The link condition and explicit combinatorial checks protect topology [9]. Cluster radius carries a cheap task-specific distance certificate. Cluster-normal memory parallels the additive spirit of appearance quadrics [3,6] but aggregates original face-normal evidence rather than per-vertex attributes. Exact rendering follows the image-driven principle [5,11]. Checkpoints, fingerprint probes, packed replays, and fail-closed integration make the resulting search executable within strict external limits.

## 3. Methodology

### 3.1 Research design: a multi-fidelity evidence loop

The methodology is a sequence of increasingly expensive filters. A candidate begins as a change to a cheap local cost or target count. It is generated deterministically, checked for structural validity, tested against a conservative distance bound, evaluated with the exact 1024-squared renderer, challenged across rotations, integrated into a byte-identical accepted parent, and only then submitted as an isolated probe. An Accepted source is fetched back from Kattis, hashed, and archived before further work. This loop converts an opaque optimization problem into a chain of falsifiable hypotheses.

![Figure 1. End-to-end pipeline from guarded geometric contraction to certified replay.](figures/pipeline.png)

The system deliberately separates three models. The **proposal model** is QEM plus cheap normal and cluster terms. The **measurement model** is the exact local renderer and independent validators. The **deployment model** is a compressed deterministic C++ replay with conservative fallbacks. Confusing these roles caused early failures: a locally good proposal was not necessarily valid, and an offline-optimal mesh was not necessarily representable under the runtime or source budget.

### 3.2 Guarded QEM backbone

For each original face \(f\), let its normalized plane be \(\pi_f=(a,b,c,d)^{\mathsf T}\) and let \(w_f\) be an optional area or curvature weight. The fundamental face quadric is [2]

$$K_f=w_f\pi_f\pi_f^{\mathsf T}.$$

The quadric at a vertex is the sum of incident face quadrics:

$$Q_v=\Sigma_{f\ni v}K_f.$$

Contracting an edge \((u,v)\) to candidate position \(p\) combines the quadrics additively. With homogeneous coordinate \(p_h=(p_x,p_y,p_z,1)^{\mathsf T}\), the geometric cost is

$$E_Q(u,v,p)=p_h^{\mathsf T}(Q_u+Q_v)p_h.$$

The implementation evaluates multiple targets: both endpoints, midpoint, cluster-size-weighted interpolation, the unconstrained three-by-three minimizer when well conditioned, and a stationary point restricted to the segment. Evaluating several targets is inexpensive relative to neighborhood validation and prevents singular or poorly conditioned quadrics from forcing an unsafe position.

The cheapest target is accepted only when every guard passes:

- **Link condition:** an interior manifold edge has exactly the two common neighbors implied by its incident triangles [9].
- **Face validity:** every surviving affected triangle has three distinct indices and area above a scale-aware threshold.
- **Orientation:** a surviving face cannot flip or exceed the configured normal deviation.
- **Combinatorics:** no self-loop, duplicate triangle, disconnected component, or edge incidence other than two is created.
- **Distance:** the propagated cluster radius stays below a conservative tolerance.

### 3.3 Cluster-radius distance certificate

Each live vertex \(v\) represents a cluster \(C(v)\) of original vertices and stores a radius \(r_v\) satisfying

$$\forall x\in C(v):\quad \lVert x-p_v\rVert_2\leq r_v.$$

If \(u\) and \(v\) contract to \(p\), the propagated radius is

$$r_{u\cup v}(p)=\max\!\left(r_u+\lVert p_u-p\rVert_2,\;r_v+\lVert p_v-p\rVert_2\right).$$

The triangle inequality proves the invariant by induction. Rejecting a contraction when this radius exceeds a safe threshold certifies the original-to-simplified direction for the assigned representatives. The simplified-to-original direction is checked against the original vertex set, and the exported mesh is independently validated using a spatial nearest-neighbor structure. A safety margin below the official tolerance absorbs floating-point evaluation, rebasing, and output quantization.

This certificate is intentionally conservative. It can reject feasible contractions because it bounds every represented original vertex by a common ball. That conservatism is valuable during bulk simplification: the exact distance checker is reserved for checkpoints, while no accepted online contraction can silently accumulate unbounded displacement.

### 3.4 Perception-aware collapse objective

Pure QEM preserves supporting planes but does not directly model the flat-normal images. We use the composite proposal cost

$$E=E_Q+\lambda_NE_N+\lambda_CE_C+\lambda_RE_R+\varepsilon_{\mathrm{tie}}.$$

The optional regularizer \(E_R\) represents case-specific curvature or visibility allocation. The infinitesimal deterministic term \(\varepsilon_{\mathrm{tie}}\) stabilizes ordering without materially changing priorities.

#### 3.4.1 Projected-area normal distortion

For every surviving affected face \(f\), let \(n_f\) and \(n'_f(p)\) be its current and candidate normals. Object-space area alone is a poor estimate of raster impact, because a large edge-on face may cover few pixels. We approximate total six-view coverage by a projected area \(A_f^{\mathrm{proj}}\) and define

$$E_N(u,v,p)=\Sigma_{f\in\mathcal A(u,v)}A_f^{\mathrm{proj}}\left[1-\operatorname{clamp}\!\left(n_f\cdot n'_f(p),-1,1\right)\right]^{\gamma}.$$

The six axial cameras make the approximation inexpensive: projected components on the three coordinate planes, with orientation-aware signs, estimate pixel influence. Alternative modes use object-space area or absolute projected components. Projected mode was most reliable on the Bunny-like case, though the exact renderer still decided final operations.

#### 3.4.2 Cluster-normal memory

Repeatedly comparing a candidate with current faces creates reference drift. Every individual step can look acceptable while the aggregate surface slowly rotates away from the original. To preserve original evidence, each live vertex carries an area-vector summary

$$S_v=\Sigma_{f\ni v}2A_fn_f,$$

initialized from original faces and merged exactly as \(S_{u\cup v}=S_u+S_v\). From the merged cluster we derive a normalized target direction \(t_{u\cup v}\). The incremental cluster penalty compares proposed incident normals with this original-support target and subtracts the current error:

$$E_C=\Sigma_f2A_f\left(\left[1-n'_f\cdot t'_{f}\right]^{\beta}-\left[1-n_f\cdot t_f\right]^{\beta}\right).$$

The subtraction is important: an operation that repairs accumulated drift can receive negative credit. The statistic is additive and independent of collapse ordering, so the reference does not degrade with the current mesh. On the 35,292-vertex case, the decisive early setting used \(\lambda_N=0.003\), \(\gamma=0.75\), projected area, \(\lambda_C=0.0001\), and \(\beta=0.5\).

![Figure 2. Cluster-normal memory retains original orientation evidence throughout a collapse tree.](figures/cluster_normal.png)

### 3.5 Exact evaluator and component diagnostics

The local CPU evaluator mirrors the public specification [1]: six positive and negative Cartesian cameras, camera distance 2.5, resolution 1024 by 1024, focal length 800, center at pixel coordinate 512, pixel-center sampling, nearest-triangle z-buffering, flat per-face normal encoding, perspective-correct reciprocal depth, 11-by-11 SSIM windows, foreground masking, and equal normal/depth weight. All release decisions use resolution 1024. Lower resolutions are allowed only as screening because 512-squared evaluations reversed several near-threshold rankings.

Three evaluators support diagnosis. The fast evaluator reports combined VPS. The component evaluator reports normal and depth SSIM for each view, identifying which signal is active. The oracle-normal variant replaces or isolates normal information to distinguish a bad topology from a bad normal surrogate. Together they prevent aggregate score improvements from hiding a catastrophic view or component.

### 3.6 Renderer-aware structural operators

#### 3.6.1 Edge flips

Two adjacent triangles \((a,b,c)\) and \((b,a,d)\) can exchange diagonal \(ab\) for \(cd\). The vertex count and vertex-set Hausdorff distance remain unchanged, but the piecewise-planar normal field can change substantially. We require positive orientation, no duplicate faces, preserved edge incidence, and unchanged closure before measuring the exact render delta. Flips are especially valuable before deletion because they can redistribute normal error and make a one-ring easier to triangulate.

#### 3.6.2 One-ring fan deletion and retriangulation

Deleting a vertex removes its incident fan and leaves a polygonal ring. For low valence, all combinatorially valid triangulations can be enumerated through Catalan recursion. For larger rings, restricted fans and dynamic-programming families reduce the search space. Each proposal is screened for orientation, duplicates, manifoldness, and distance before rendering. A deletion is accepted only when the exact worst relevant view remains above the chosen margin.

#### 3.6.3 Coordinate and tangent fitting

After topology stabilizes, selected vertices are perturbed in coordinate, tangent, or normal directions. The objective is exact combined SSIM subject to triangle-area, orientation, and distance guards. Finite-difference or coordinate-search steps operate only inside bounded trust regions. Foreground cropping reduces cost on the million-vertex branch. Coordinate fitting is a margin-restoration tool, not a substitute for topology search: large gains generally required a flip or deletion first.

### 3.7 Discovering hidden constraints without hidden inputs

Kattis exposed only an overall verdict and score, not per-case render components or hidden meshes. We treated this as a black-box experimental-design problem. No hidden file was downloaded, reconstructed, or placed in the public artifact. The submitted program could, however, compute statistics of its own input at runtime, and the displayed score could reveal a carefully controlled output count.

#### 3.7.1 Isolation and count decoding

A useful probe changes exactly one case \(k\), while the other five outputs are byte-identical to a known Accepted parent. If the displayed score is \(S\) and all other output counts are known, the changed count is recovered as

$$|V'_k|=\operatorname{round}\!\left(|V_k|\left[6\left(1-\frac{S}{100}\right)-\Sigma_{i\neq k}\frac{|V'_i|}{|V_i|}\right]\right).$$

For the 23,201-vertex case, one output vertex changes the score by approximately 0.00071836, so six displayed decimals uniquely determine the integer. This converts an aggregate score into a case-specific measurement while preserving scientific control: only one unknown changes.

#### 3.7.2 Safe diagnostic payloads

Let \(B_k\) be a conservative accepted base count for case \(k\), and let \(h\) be a small nonnegative diagnostic integer. A valid diagnostic branch outputs its normal manifold mesh plus unused duplicate vertices or another structurally harmless carrier, giving

$$|V'_k|=B_k+h.$$

Because faces reference only the certified base geometry, the active surface and its render remain unchanged. The added count reduces the diagnostic score but does not endanger the stored leaderboard best. After Kattis reports the score, the equation above recovers \(h\). Multi-part fingerprints are emitted as base-\(B\) digits,

$$H=d_0+Bd_1+B^2d_2,$$

using separate isolated probes when the safe count slack cannot carry the entire integer.

#### 3.7.3 Rotation-invariant fingerprints

Proxy identity cannot be inferred reliably from vertex count alone. We used invariants computed from normalized geometry: quantized bounding-box ratios, genus and component count, sorted normalized edge-length multisets, edge-graph labels, orientation signatures, and order-sensitive hashes when raw ordering mattered. A representative edge statistic quantizes

$$q_e=\operatorname{round}\!\left(B_q\frac{\lVert p_u-p_v\rVert_2}{D_{\mathrm{AABB}}}\right),$$

then hashes the sorted multiset of \(q_e\). Normalization removes translation and uniform scale; lengths remove rotation; sorting removes edge enumeration order. For the 23,201-vertex branch, three recovered base-15,776 digits were 2,382, 13,367, and 11,736, matching the exact-source proxy and differing from aggregate proxies. Independent pose, edge-graph, orientation, scale, and ordering codes agreed. This supplied evidence for selecting the correct public proxy without exposing the hidden coordinates.

For the 35,292-vertex case, a conservative diagnostic output encoded normalized bounding-box statistics. The recovered values matched the public Bunny proxy, while acceptance behavior still showed a proxy-to-hidden renderer gap. We therefore separated **identity evidence** from **acceptance evidence**: a matching invariant justifies which geometry family to study, but only isolated Kattis probes certify a hidden threshold.

#### 3.7.4 Frontier search and uncertainty control

Once a branch was identified, target counts were reduced by binary or bracketed search. Each probe changed one target or one operator family. A pass moved the certified frontier downward; a failure tightened the lower bound. Near runtime limits, repeated identical submissions distinguished stochastic load from systematic timeout, but repeated failures triggered engineering rather than indefinite resubmission.

Public proxies were rotated through sixteen orientations when practical. A candidate had to pass structural, distance, and exact-render checks across selected rotations, with an explicit margin. Rotation 04 was especially predictive for Lucy and Nefertiti during early tuning, but no single rotation was treated as authoritative. Hidden acceptance boundaries were recorded separately from local proxy thresholds to avoid retroactively labeling proxy scores as official measurements.

### 3.8 Case-specialized solving

The score formula and different perceptual bottlenecks make one universal retained ratio suboptimal. Exact input size and conservative signatures choose a common core with specialized schedules.

| Case | Input vertices | Final vertices | Retained | Dominant bottleneck and final strategy |
|---|---:|---:|---:|---|
| Sphere-like | 4,098 | 25 | 0.6101% | Analytic smooth structure; guarded tiny branch |
| Armadillo | 23,201 | 4,340 | 18.7061% | Flat-normal detail; canonical replay, flips, fitted coordinates |
| Bunny-like | 35,292 | 2,839 | 8.0443% | Normal drift and silhouette; cluster memory, fan deletion, fitting |
| Lucy | 49,987 | 3,030 | 6.0616% | Tight normal/SSIM frontier; repeated flip-delete-fit cycles |
| Slender | 377,084 | 7,400 | 1.9624% | Curvature allocation and runtime; fitted replay, 40 flips, cache |
| Nefertiti | 1,009,118 | 16,500 | 1.6351% | Runtime, memory, facial detail; staged QEM and cropped fitting |

The Sphere-like case admits a tiny guarded structural construction. Armadillo consumes the largest score loss and therefore received the most renderer-aware structural work. Bunny motivated cluster-normal memory. Lucy required careful operation ordering because local Hausdorff headroom did not imply normal-map headroom. Slender was geometrically compressible but exposed a runtime failure until reference caching. Nefertiti required staged passes, bounded heaps, and cropped fitting to stay within compute and memory limits.

### 3.9 Checkpoints, canonical replay, and transaction safety

Offline search can take hours; the submitted program cannot. We compile discoveries into deterministic operation streams. Vertex references are expressed in canonical current-mesh order or through stable local signatures rather than fragile temporary indices. Topology operators, fan choices, and quantized coordinate residuals are bit-packed. Repeated instruction fragments share dictionaries or base-94 payloads.

Replay is transactional. Before a specialized block, the current mesh is either checkpointed or reconstructible from a judge-proven prefix. Every deletion, flip, and coordinate edit checks its expected local topology. A block commits only if its final count and structural invariants match; otherwise it rolls back or retains the accepted checkpoint. This behavior makes the final source a compact constructive certificate rather than an unchecked precomputed mesh dump.

### 3.10 Runtime, memory, and source-budget engineering

Pointer-heavy half-edge structures were avoided on the largest input. Positions, faces, active flags, quadrics, cluster statistics, versions, and initial incidence live in flat arrays. New incidences are appended to compact linked storage. A binary heap stores candidate edges; endpoint versions invalidate stale entries lazily. Periodic pruning bounds heap growth. Input uses a buffered character parser and output uses a one-mebibyte writer.

Rebase checkpoints rebuild exact current adjacency and quadrics, controlling numerical drift and allowing later stages to use smaller identifiers. Foreground crop boxes limit renderer work. Reference caches prevent repeated projection and normalization on Slender and Nefertiti. Static storage reduced allocator overhead and, on Kattis, improved runtime stability compared with otherwise equivalent dynamic variants.

The 131,072-byte source limit became a second optimization problem. Replay payloads were entropy-coded, recurring code fragments were factored, safe whitespace was removed, and unused diagnostic branches were stripped. Every minified candidate was compiled independently and run on the official sample. The final fetched-back source is 130,973 bytes, leaving 99 bytes of headroom.

### 3.11 Fail-closed validation funnel

Candidate validation proceeds from cheapest to most authoritative. First, parse and index checks reject malformed meshes. Structural validation checks positive triangle area, duplicate faces, connectedness, oriented edge multiplicity, and the closed two-manifold condition. Next, an independent spatial search checks both directions of the official vertex-Hausdorff distance. Only then does the exact renderer compute normal, depth, per-view, and combined SSIM. Rotation suites test proxy robustness. Integration compares all untouched outputs byte-for-byte with the accepted parent. Kattis is the final external gate.

![Figure 3. Fail-closed validation funnel used before an informative Kattis probe.](figures/validation_funnel.png)

The funnel is intentionally asymmetric: a candidate can be rejected for uncertainty even if it might pass. Near a hard threshold, false acceptance costs a submission and confounds diagnosis, whereas a conservative local rejection merely postpones a possible gain. This bias was appropriate because many branches had only a few thousandths of SSIM margin.

## 4. Results and Discussions

### 4.1 Official result and exact score reconstruction

Submission 20082703 was reported by Kattis as Accepted, 7/7, score 93.830074. The fetched-back source is 130,973 bytes and has SHA-256 digest `9195d42a73a6b85c8ae30d731f532175bdcd7c2982d421143d631b4c64b1a92c`. The team finished second. These official and reconstructed records are archived in the public artifact [14].

| Case | Original vertices | Final vertices | Retained | Reduction | Score loss |
|---|---:|---:|---:|---:|---:|
| Sphere-like | 4,098 | 25 | 0.610054% | 99.389946% | 0.101676 |
| Armadillo | 23,201 | 4,340 | 18.706090% | 81.293910% | 3.117682 |
| Bunny-like | 35,292 | 2,839 | 8.044316% | 91.955684% | 1.340719 |
| Lucy | 49,987 | 3,030 | 6.061576% | 93.938424% | 1.010263 |
| Slender | 377,084 | 7,400 | 1.962427% | 98.037573% | 0.327071 |
| Nefertiti | 1,009,118 | 16,500 | 1.635091% | 98.364909% | 0.272515 |
| Sum / final | 1,500,780 | 34,134 | \(\Sigma\rho=.3702\) | — | 93.830074 |

Substitution into the official score definition gives

$$100\left[1-\frac{1}{6}\left(\frac{25}{4098}+\frac{4340}{23201}+\frac{2839}{35292}+\frac{3030}{49987}+\frac{7400}{377084}+\frac{16500}{1009118}\right)\right].$$

The exact value is

$$93.83007422510956,$$

which rounds to the displayed score. The global compression ratio, \(1-34134/1500780\), is not the leaderboard formula; the official mean gives every case equal weight [1].

![Figure 4. Final retained ratios, compression, and per-case score contribution.](figures/final_results.png)

### 4.2 Accepted progression

Development contained many small accepted steps and several structural jumps. The selected milestones below are official or archived high-waters; the public artifact contains the data used to generate the progression figure [14].

| Submission or stage | Official score | Main change | Output-count evidence |
|---|---:|---|---|
| Generic high-water | 81.945906 | Guarded generic simplification plateau | Archived lineage |
| Standalone QEM | 86.998654 | Case-routed robust QEM replaces legacy pipeline | Archived lineage |
| 19932621 | 87.913148 | Cluster-normal accumulation; Bunny-like 15.5% | Accepted 7/7 |
| 20020565 / 20021691 lineage | 91.801173 | Renderer-aware cycles and compressed replay | Accepted counts archived |
| 20025898 | 92.932731 | Armadillo, Bunny, Lucy, and Nefertiti checkpoints | [25, 5136, 3400, 3088, 7900, 17660] |
| 20039214 | 93.600980 | Renderer-optimized Armadillo at 4,570 | Accepted 7/7 |
| 20051927 | 93.769981 | Armadillo 4,340 and Bunny 2,848 | Accepted 7/7 |
| 20082128 | 93.812395 | Bunny 2,839, Lucy 3,030, Nefertiti 16,500 | Accepted 7/7 |
| 20082703 | 93.830074 | Slender 7,400 with 40 flips and reference cache | [25, 4340, 2839, 3030, 7400, 16500] |

![Figure 5. Official high-water trajectory from the generic plateau to the final submission.](figures/score_progression.png)

The trajectory supports a central claim: small parameter tuning alone did not explain the final gain. The large jumps came from representation and workflow changes—cluster memory, structural replay, exact proxy identification, renderer-aware deletion, and replay compression. Later progress became incremental because each branch approached a different active constraint.

### 4.3 Hidden frontier reconstruction

The black-box protocol produced useful pass/fail boundaries before the later structural replays. These numbers are not direct hidden-mesh measurements; they are acceptance observations from isolated submissions.

| Branch | Accepted frontier observation | Rejected neighboring observation | Interpretation |
|---|---|---|---|
| Armadillo, early QEM | 34%, 33%, and 32% retained passed | 31% retained failed | Pure QEM normal frontier near 31–32% |
| Bunny-like, early normal QEM | 16% passed | Lower local-only settings failed | Current-face normal cost accumulated drift |
| Bunny-like, cluster memory | 15.5% passed | 14.5% new-cost branch failed hidden | Cluster memory improved hidden frontier but proxy remained optimistic |
| Lucy, early QEM | 9.5% passed | 9.0% local and clustered variants failed | Hidden normal margin tighter than Hausdorff margin |
| Nefertiti, early QEM | 4% passed | 3% combined probe failed | Runtime and perceptual constraints both active |

These boundaries guided allocation rather than defining the final outputs. Offline replay eventually moved far below the early QEM ratios: Armadillo to 18.706%, Bunny-like to 8.044%, Lucy to 6.062%, and Nefertiti to 1.635%. The comparison quantifies the value of structural and coordinate search beyond a better collapse weight.

The proxy audit also changed the research plan. Edge fingerprints showed that the exact 23,201-vertex source, not the agg3 or agg5 alternatives, matched the submitted input statistics. This eliminated a major source of proxy error and justified offline replay on the correct pose. The Bunny identity diagnostic matched the public family, yet its local rotations passed more aggressive counts than hidden Kattis. Identity and acceptance were therefore treated as separate uncertainties.

### 4.4 Case-by-case discussion

#### 4.4.1 Sphere-like branch

The 4,098-vertex smooth input has a simple global shape and weak high-frequency detail. A guarded structural construction reaches 25 vertices while preserving closure and the official render. Its small original size makes each retained vertex relatively valuable, but the branch was solved early and did not dominate later effort.

#### 4.4.2 Armadillo

Armadillo contributes the largest final score loss because it retains 18.706% of its input. Early generic QEM passed around 32% and failed at 31%, demonstrating that flat-normal detail rather than the nominal Hausdorff bound was active. Exact-source fingerprinting removed pose uncertainty. Canonical topology replay, low-valence deletion, edge flips, and tangent/coordinate fitting then reduced accepted counts through 5,136, 4,570, and 4,340. A 4,339 probe also passed, but the final Slender integration was rebased on the byte-identical 4,340 lineage; the one-vertex difference is approximately 0.000718 points and was not worth destabilizing the release.

A notable local result was that the 4,570-vertex candidate slightly exceeded the combined VPS of a 5,136-vertex accepted reference despite using 566 fewer vertices. This is direct evidence that vertex count and perceptual quality are not monotonically coupled when topology and flat normals change together.

#### 4.4.3 Bunny-like branch

Bunny motivated cluster-normal memory. Public rotations could pass lower ratios than hidden Kattis, but the local-current normal term became unreliable over long collapse chains. Carrying original support normals improved the official high-water to 87.913148 and certified 15.5% retained. Later renderer-aware replay reduced accepted counts through 3,900, 3,850, 3,800, 3,775, 3,750, 3,400, 2,856, 2,848, and finally 2,839. The final topology is protected by checksum and fail-closed replay.

The decisive lesson is not that a larger normal weight is always better. Excessive weights froze locally visible faces and wasted vertices. The useful change was preserving the reference distribution across collapse history, allowing a modest normal term to remain meaningful late in simplification.

#### 4.4.4 Lucy

Lucy exposed operation ordering. Direct 9% QEM and early clustered variants failed, whereas 9.5% passed. Subsequent renderer-ranked flips, coordinate fitting, and independent fan deletions advanced through 3,200, 3,150, 3,100, 3,097, 3,088, 3,040, and 3,030 vertices. The final 3,030 proxy measured combined VPS 0.90929170 and Hausdorff-to-tolerance ratio 0.935986. The remaining geometric headroom was larger than the remaining normal-map headroom, confirming that distance alone could not rank the final deletions.

#### 4.4.5 Slender / Yeah Right

The 377,084-vertex model can be reduced aggressively, but generic QEM allocated triangles poorly around curvature. Plane quadrics were weighted by a mild curvature factor,

$$w_f=1+\lambda_{\kappa}\kappa_f^{\eta},\qquad \lambda_{\kappa}=2,\quad \eta=0.035,$$

followed by a fitted collapse replay and 40 topology flips. The 7,800 branch passed; the 7,400 geometry also passed visual checks but initially timed out on test 6 in submission 20082666. A reference cache removed repeated work, and otherwise equivalent submission 20082703 passed 7/7. This is the cleanest systems ablation in the campaign: geometry was already valid, while cache architecture determined acceptance.

#### 4.4.6 Nefertiti

The million-vertex input dominated compute and memory. The accepted schedule begins with a 2.2% QEM stage, rebuilds incidence and quadrics, applies several normal-sensitive pre-passes, and finishes with late passes plus bounded multi-view coordinate fitting. Foreground cropping and heap caps were essential. Accepted counts progressed through 18,669, 17,660, 17,100, 17,000, 16,900, 16,800, 16,700, 16,650, and 16,500. Lower geometrically valid candidates could fail SSIM or runtime. This branch shows why a feasible offline mesh and a deployable contest program are different objects.

### 4.5 Ablation evidence

Only comparisons with an otherwise stable lineage are treated as clean ablations.

- **Cluster-normal memory:** carrying additive original-face normal evidence reduced the hidden Bunny-like retained ratio beyond the preceding current-normal-only branch and raised the accepted high-water to 87.913148. The mechanism isolates reference drift rather than merely a target-count change.
- **Renderer-aware topology:** installing validated Armadillo and Bunny replays produced a discontinuous score jump from the low 92 range to 92.932731 while preserving other accepted outputs. The gain cannot be explained by a QEM weight change alone.
- **Armadillo structural optimization:** reducing 5,136 to 4,570 vertices raised the official score to 93.600980 while maintaining or improving local VPS, demonstrating topology-dependent non-monotonicity.
- **Nefertiti fitting at fixed count:** a broader first normal-fitting pass enabled an accepted 16,700 branch where an unfitted candidate at the same count failed. Count is therefore not a sufficient experimental variable.
- **Slender reference cache:** 7,400-vertex geometry failed runtime without the cache and passed 7/7 with it. This isolates systems architecture as the decisive factor.

### 4.6 Negative results and what they taught us

Negative results shaped the final model and are as informative as accepted milestones.

- Fast generic libraries, MeshLab variants, meshoptimizer, and off-the-shelf simplifiers produced valid meshes but inferior flat-normal SSIM at aggressive ratios. Their objectives did not match the fixed renderer closely enough.
- Isotropic remeshing often preserved smooth depth and silhouettes but produced poor flat-shaded normal maps. Uniform triangle quality is not equivalent to renderer-optimal face orientation.
- Shared-vertex multilayer triangulations and alternative surface decompositions passed structural checks yet changed z-buffer ownership, causing SSIM loss.
- Low-resolution rendering was misleading near the threshold. Resolution 512 sometimes reversed candidate ordering relative to the official 1024, so it was demoted to screening only.
- Coordinate-only optimization produced small improvements but could not cross large plateaus. Topology had to change first; fitting restored margin afterward.
- Proxy success did not guarantee hidden acceptance. Sixteen public rotations reduced orientation overfitting, but isolated Kattis probes remained necessary.
- Heap growth, repeated incidence reconstruction, full-frame fitting, and redundant reference work caused timeouts on the largest cases. Heap pruning, crop bounds, static storage, and caching were necessary even when geometry was unchanged.
- Source minification introduced independent failures, including missing standard headers, altered storage duration, and accidental payload damage. Every release candidate was recompiled, sample-tested, byte-counted, and fetched back after acceptance.

### 4.7 Statistical and experimental interpretation

This was an optimization competition over six fixed instances, not a population study. We therefore do not report confidence intervals over independent meshes. The appropriate evidence is a traceable sequence of deterministic local measurements and official pass/fail observations. Multi-rotation proxy evaluation acts as a stress test, not as sampling from the hidden distribution. Kattis repetitions near runtime boundaries estimate operational stability but do not turn the benchmark into a stochastic trial.

The score is also an incomplete quality measure. It ranks vertex counts only after all hard constraints pass. Two accepted meshes with the same count receive the same score even if one has a larger SSIM margin. Consequently, exact local VPS and distance ratios are used as safety diagnostics, while official counts determine leaderboard value.

### 4.8 Threats to validity and limitations

**Internal validity.** Many late improvements combined code compression, runtime work, and geometry. We label only stable-parent comparisons as clean ablations. Hidden feedback was aggregate, so isolation and byte-identical unchanged outputs were required to identify the active case.

**Construct validity.** The official metric uses six views and flat normal/depth SSIM [1]. It does not measure arbitrary camera paths, materials, texture, or semantic plausibility. A mesh optimized for this evaluator may allocate triangles differently from a general interactive asset.

**External validity.** The final contest program specializes to fixed cases and embeds replay streams. It is not a universal drop-in simplifier. The reusable contributions are the guarded QEM core, cluster-normal memory, exact-render calibration, structural search, fingerprint methodology, and fail-closed deployment workflow.

**Geometric guarantee.** The radius certificate is conservative and task-specific to vertex-set Hausdorff. Applications requiring continuous surface error should use sampling or envelope methods such as [7,8].

**Black-box inference.** Diagnostic fingerprints identify proxy families and parameters computed by the submitted program; they do not reveal the complete hidden input. A hash match can collide, so independent signatures and acceptance behavior are required. The public release excludes hidden inputs and third-party meshes.

### 4.9 Reproducibility and code availability

The complete curated artifact is available in the [public GitHub repository](https://github.com/lequoctran181/IMC_Challenge_2026) [14]. It contains the byte-exact source fetched back from Kattis, machine-readable result metadata, the readable QEM and cluster-normal research core, evaluator sources, score reconstruction, this DOCX/PDF article, figures, build scripts, and a fail-closed release verifier. Official or hidden meshes are not redistributed.

The authoritative release record is:

- Kattis submission: 20082703
- Judgement: Accepted, 7/7
- Official score: 93.830074
- Exact reconstructed score: 93.83007422510956
- Final output counts: [25, 4340, 2839, 3030, 7400, 16500]
- Source size: 130,973 bytes
- Source SHA-256: `9195d42a73a6b85c8ae30d731f532175bdcd7c2982d421143d631b4c64b1a92c`

The release verifier reconstructs the score, checks source and article checksums, confirms the source limit, and compiles the fetched-back C++17 program. Researchers may build the readable simplifier and evaluator on any legally obtained triangular mesh; exact official integrity verification needs no contest mesh.

### 4.10 Model summary and operational workflow

The final system can be summarized as a compiler from a dense mesh to a certified perceptual replay:

- Parse and normalize the input; recognize the case with fail-closed signatures.
- Build quadrics, incidence, original normal summaries, cluster radii, and a versioned edge heap.
- Pop the cheapest candidate, test link, orientation, duplicate, normal, and radius conditions, then commit the contraction.
- Rebase at certified checkpoints to rebuild exact adjacency and control numerical drift.
- Evaluate checkpoints with independent topology, distance, and exact 1024-squared component renderers.
- Search renderer-aware flips, fan deletions, and coordinate residuals in small neighborhoods.
- Canonicalize and compress accepted operations; replay them transactionally from a judge-proven prefix.
- Cache large-case references, compact live vertices and faces, and stream the final OBJ output.
- Fetch every Accepted high-water, record counts and hashes, and use it as the only integration parent.

This workflow is neither purely online nor a hard-coded output. The online component guarantees scalability and a general safe structure. The offline component exploits the fixed evaluator. The replay component makes expensive discovery feasible under contest limits. The evidence protocol makes hidden feedback informative without pretending that proxy measurements are official.

## 5. Conclusion

The challenge was solved not by a single metric but by aligning every layer with the evaluator and deployment constraints. QEM supplied speed and geometric coherence [2]. Link and radius guards protected topology and distance. Projected normal cost and cluster-normal memory addressed the dominant perceptual failure of long collapse sequences. Exact rendering found structural and coordinate edits invisible to local geometry. Controlled black-box probes reduced uncertainty about hidden identities and acceptance frontiers. Compression, caching, and transactional replay converted offline discoveries into a 130,973-byte program that passed within 21 seconds.

The final 93.830074 score, second-place finish, and exact count reconstruction validate the hybrid strategy. The broader lesson is methodological: maintain a conservative invariant during large-scale optimization, measure the real rendered objective at carefully chosen checkpoints, isolate every external experiment, and deploy only deterministic operations backed by traceable evidence. Under a hard perceptual threshold, the strongest solver is not the one with the most elaborate local cost; it is the one that couples a scalable proposal model, an exact measurement model, and a fail-closed evidence and deployment pipeline.

## References

[1] IMC Challenge, “Problem B: Perception-Aware Lossless Simplification of Million-Vertex 3D Meshes for Mobile Platforms,” official Kattis statement, 2026. [Problem statement](https://imc2.kattis.com/contests/imc2-2/problems/simplifygeometry).

[2] M. Garland and P. S. Heckbert, “Surface Simplification Using Quadric Error Metrics,” *Proceedings of SIGGRAPH 1997*, pp. 209–216, 1997. [Author manuscript](https://www.cs.cmu.edu/~garland/Papers/quadrics.pdf), [DOI](https://doi.org/10.1145/258734.258849).

[3] M. Garland and P. S. Heckbert, “Simplifying Surfaces with Color and Texture Using Quadric Error Metrics,” *Proceedings of IEEE Visualization 1998*, pp. 263–269, 1998. [DOI](https://doi.org/10.1109/VISUAL.1998.745312).

[4] H. Hoppe, “Progressive Meshes,” *Proceedings of SIGGRAPH 1996*, pp. 99–108, 1996. [Author project page](https://hhoppe.com/proj/pm/), [DOI](https://doi.org/10.1145/237170.237216).

[5] P. Lindstrom and G. Turk, “Image-Driven Simplification,” *ACM Transactions on Graphics*, vol. 19, no. 3, pp. 204–241, 2000. [Author manuscript](https://faculty.cc.gatech.edu/~turk/my_papers/image_simp_tog2000.pdf), [DOI](https://doi.org/10.1145/353981.353995).

[6] H. Hoppe, “New Quadric Metric for Simplifying Meshes with Appearance Attributes,” *Proceedings of IEEE Visualization 1999*, pp. 59–66, 1999. [Author manuscript](https://www.hhoppe.com/newqem.pdf), [DOI](https://doi.org/10.1109/VISUAL.1999.809869).

[7] J. D. Cohen, A. Varshney, D. Manocha, G. Turk, H. Weber, P. K. Agarwal, F. P. Brooks Jr., and W. V. Wright, “Simplification Envelopes,” *Proceedings of SIGGRAPH 1996*, pp. 119–128, 1996. [DOI](https://doi.org/10.1145/237170.237220).

[8] P. Cignoni, C. Rocchini, and R. Scopigno, “Metro: Measuring Error on Simplified Surfaces,” *Computer Graphics Forum*, vol. 17, no. 2, pp. 167–174, 1998. [DOI](https://doi.org/10.1111/1467-8659.00236).

[9] T. K. Dey, H. Edelsbrunner, S. Guha, and D. V. Nekhayev, “Topology Preserving Edge Contraction,” *Publications de l'Institut Mathématique*, vol. 66, pp. 23–45, 1999. [Open-access record](https://research-explorer.ista.ac.at/record/3582).

[10] Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, “Image Quality Assessment: From Error Visibility to Structural Similarity,” *IEEE Transactions on Image Processing*, vol. 13, no. 4, pp. 600–612, 2004. [Author page](https://ece.uwaterloo.ca/~z70wang/publications/ssim.html), [DOI](https://doi.org/10.1109/TIP.2003.819861).

[11] J. Hasselgren, J. Munkberg, J. Lehtinen, M. Aittala, and S. Laine, “Appearance-Driven Automatic 3D Model Simplification,” *Eurographics Symposium on Rendering*, pp. 85–97, 2021. [Eurographics record](https://diglib.eg.org/items/82ec8df1-d62a-40f0-a007-90374b5140dc), [DOI](https://doi.org/10.2312/sr.20211293).

[12] H.-T. D. Liu, M. Gillespie, B. Chislett, N. Sharp, A. Jacobson, and K. Crane, “Surface Simplification Using Intrinsic Error Metrics,” *ACM Transactions on Graphics*, vol. 42, no. 4, article 118, 2023. [Author project page](https://markjgillespie.com/Research/intrinsic-coarsening/index.html), [DOI](https://doi.org/10.1145/3592403).

[13] H.-T. D. Liu, M. Rahimzadeh, and V. Zordan, “Controlling Quadric Error Simplification with Line Quadrics,” *Computer Graphics Forum*, vol. 44, no. 5, article e70184, 2025. [Eurographics record](https://diglib.eg.org/items/2117a7d7-e66b-417b-8b55-14334e9e237f), [DOI](https://doi.org/10.1111/cgf.70184).

[14] NEU.AddictedTribes, “Certified Perceptual Mesh Simplification: IMC Challenge 2026 Research Artifact,” GitHub, 2026. [Public repository](https://github.com/lequoctran181/IMC_Challenge_2026).

### Appendix A. Core simplification pseudocode

```text
load mesh M
build quadrics Q, original cluster normals S, and radius certificates r
build flat incidence arrays and a versioned candidate heap

while activeVertexCount exceeds checkpointTarget:
    (u, v) <- cheapest current edge candidate
    generate endpoint, midpoint, weighted, segment, and analytic positions
    order positions by QEM, projected-normal, and cluster-normal cost
    for position p in that order:
        if linkCondition(u, v) and validAffectedFaces(p)
           and normalGuard(p) and radiusCertificate(p) <= safeTolerance:
            commit contraction of u into v at p
            merge quadrics, original normal summaries, and radii
            refresh affected heap entries
            break

for each rendererOptimizedTransaction:
    decode flips, fan deletions, and coordinate residuals
    apply only when canonical local preconditions match
    commit only if expected count and structural checks pass
    otherwise retain the judgeProven checkpoint

compact live vertices and faces
stream the OBJ output
```

### Appendix B. Principal parameters in the final lineage

| Branch | Base stage | Normal term | Additional mechanism | Offline tail |
|---|---|---|---|---|
| Sphere-like | Guarded structural reduction | Not required | Exact structural checks | Target 25 |
| Armadillo | Staged QEM checkpoints | \(\lambda_N=.01\), \(\gamma=.75\) | Tangent and coordinate fitting | Canonical replay to 4,340 |
| Bunny-like | Checkpoint near 8.043% | \(\lambda_N=.003\), \(\gamma=.75\) | \(\lambda_C=10^{-4}\), \(\beta=.5\) | Flips, deletions, fit to 2,839 |
| Lucy | 7.448% base then replay | \(\lambda_N=.001\), \(\gamma=1\) | Repeated fit and mesh refresh | Deletion and flip streams to 3,030 |
| Slender | 2.2% curvature-aware QEM | Curvature allocation | \(\lambda_{\kappa}=2\), \(\eta=.035\) | Fitted replay and 40 flips to 7,400 |
| Nefertiti | 2.2% first QEM stage | \(\lambda_N=.0003\), \(\gamma=.75\) | Five pre-passes and two late passes | Cropped fitting to 16,500 |

Parameters were calibrated using exact local rendering and frozen once a branch became judge-proven. They are task-specific settings, not universal constants.

### Appendix C. Verification checklist

- Compile the fetched-back source as C++17 and confirm that its size does not exceed 131,072 bytes.
- Verify the SHA-256 digest before treating the source as authoritative.
- Run the official sample and confirm that the first output line is `8 12`.
- Validate indices, positive triangle areas, duplicate faces, oriented edge multiplicity, connectivity, and closed two-manifold structure.
- Evaluate both directions of the official vertex-Hausdorff distance against five percent of the original bounding-box diagonal.
- Render all six views at 1024 by 1024; do not substitute 512 by 512 for a release decision.
- Record normal, depth, per-view, and combined SSIM rather than relying on one aggregate number.
- Archive and checksum each Accepted source; byte-compare unchanged cases with the integration parent.
