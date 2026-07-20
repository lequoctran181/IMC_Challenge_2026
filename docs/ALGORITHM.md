# Algorithm and mathematical model

This chapter is the concise, web-readable counterpart of the [full Round 2 article](../paper/IMC_Challenge_Round2_NEU_AddictedTribes.pdf). It separates the general method from the case-specific replay data used by submission 20082703.

## 1. Optimization problem

For each input triangular mesh $M=(V,F)$, produce a closed triangular 2-manifold $M'=(V',F')$ that minimizes its vertex count $n(M')$ subject to:

- valid indices, non-degenerate triangles, and exactly two incident faces per edge;
- symmetric vertex-set Hausdorff distance

  $$d_H(M,M') \le 0.05\,\operatorname{diag}(\operatorname{AABB}(M));$$

- perceptual score at least 0.9 under the official six-view normal/depth renderer.

For view $c$, the perceptual score is

$$S_c=\tfrac12\operatorname{SSIM}(I^N_c,I^{N'}_c)+\tfrac12\operatorname{SSIM}(I^D_c,I^{D'}_c),$$

and the test score is the mean of $S_c$ over the six axis-aligned cameras. Among valid outputs, ranking depends only on retained vertex ratios:

Let $N_i$ and $M_i$ denote the input and output vertex counts. Then

$$\operatorname{Score}=100\left(1-\frac16\sum_{i=1}^6\frac{M_i}{N_i}\right).$$

The objective is discontinuous: topology and visibility changes can modify many image pixels at once, while the score rewards each feasible removed vertex according to the case denominator. We therefore use a hierarchy of increasingly expensive models rather than attempting one monolithic optimizer.

## 2. Guarded QEM core

Each active vertex stores a symmetric homogeneous quadric using ten coefficients. For a plane $p=(a,b,c,d)^T$ and weight $w$, the face contribution is $wpp^T$. Contracting edge $(u,v)$ combines quadrics:

$$Q_{uv}=Q_u+Q_v,\qquad E_Q(x)=\tilde{x}^{T}Q_{uv}\tilde{x},\quad \tilde{x}=(x,y,z,1)^T.$$

The implementation tests several candidate positions:

- $x_u$ and $x_v$;
- midpoint $(x_u+x_v)/2$;
- cluster-size-weighted interpolation;
- the unconstrained solution of the $3\times3$ quadratic system when nonsingular;
- a stationary point on the segment $x(t)=x_u+t(x_v-x_u)$ when it lies in $[0,1]$.

The lowest-cost candidate is accepted only if all guards pass:

1. **Link condition:** an interior manifold edge must have exactly two common one-ring neighbors.
2. **Orientation guard:** every affected surviving face keeps a positive area and sufficient agreement with its previous normal.
3. **Combinatorial guard:** no self-loop, duplicate triangle, or edge with incidence other than two may be created.
4. **Distance certificate:** the conservative cluster radius remains within the official tolerance.

For a contraction to $x$, the propagated radius is

$$r_{uv}(x)=\max\{r_u+\|x-x_u\|,\ r_v+\|x-x_v\|\}.$$

By induction, every original vertex represented by the merged cluster lies inside this radius around $x$. This provides a cheap fail-closed certificate; an independent checker still validates exported candidates.

## 3. Projected-area normal penalty

Pure QEM preserves planes but can rotate face normals in visually important regions. For every affected face $f$, let $n_f$ and $n'_f(x)$ be its old and proposed unit normals. We add a nonlinear penalty

$$E_N(x)=\sum_{f\in\mathcal A(u,v)} A_f^{\mathrm{proj}}\bigl(1-\operatorname{clamp}(n_f\cdot n'_f(x),-1,1)\bigr)^p.$$

For the orthographic no-occlusion approximation, the six-camera total is

$$A_f^{\mathrm{proj}}=2A_f\left(|n_x|+|n_y|+|n_z|\right)=2A_f\|n_f\|_1.$$

Opposite camera pairs have equal magnitudes. The implementation uses three face-oriented projected magnitudes with perspective depth and absorbs the constant factor two into $\lambda_N$. This term does not model occlusion, clipping, silhouette ownership, or actual visibility; the local renderer supplies those effects at checkpoints.

The general collapse objective is

$$E(x)=E_Q(x)+\lambda_N E_N(x)+\lambda_C E_C(x)+\lambda_R E_R(x),$$

where $E_C$ is cluster-normal memory and $E_R$ denotes curvature or projected-importance regularizers enabled only in selected regimes.

## 4. Cluster-normal memory

Repeated local normal comparisons gradually forget the original surface. Let $F_0$ be the original faces and $a_f=(p_{f,2}-p_{f,1})\times(p_{f,3}-p_{f,1})=2A_fn_f$. Initially, every vertex receives one $a_f$ for every original vertex-face incidence. If live vertex $v$ represents the original-vertex multiset $C(v)$, then

$$S_v=\sum_{x\in C(v)}\sum_{f\in F_0}\mathbf 1[x\in f]a_f.$$

When $u$ contracts into $v$, the statistic is merged exactly:

$$S_v\leftarrow S_u+S_v.$$

For every support vector above the norm threshold, define $u_v=S_v/\|S_v\|$. For an affected face $g=(v_1,v_2,v_3)$, the implementation forms $T_g=\sum_j u_{v_j}$. In the proposed state, either contracted endpoint is replaced with the normalized merged support vector, producing $T'_g$. With current and proposed face area vectors $b_g,b'_g$, mode 2 uses

$$E_C=\sum_g\|b_g\|\left[\phi(b'_g,T'_g)-\phi(b_g,T_g)\right],\quad
\phi(b,T)=\left(1-\operatorname{clamp}\frac{b\cdot T}{\|b\|\|T\|}\right)^q.$$

The subtraction makes this an incremental loss. A zero or near-zero stored/target vector has no defined direction and causes the candidate cluster term to fail closed. Associativity gives merge-tree invariance only when the final original-vertex cluster multiset is the same; different partitions, topology, surviving identities, or geometry can produce different penalties. This change was decisive on the Bunny-like hidden case and became part of the readable research core in [`src/research/compact_qem_lab.cpp`](../src/research/compact_qem_lab.cpp).

The reported support drift is the active-face, area-weighted mean angle between $b_g$ and $T_g$, sampled initially, every 500 collapses, and at the final target. Undefined targets are excluded and counted. This diagnostic is mechanism-aligned; six-view combined SSIM remains the independent rendered outcome.

![Original face evidence is accumulated through the collapse tree](../paper/figures/cluster_normal.png)

## 5. Specification-matching local renderer as an optimization oracle

The local evaluator implements the public specification:

- cameras at the positive/negative Cartesian axes, distance 2.5;
- 1024×1024 resolution, focal length 800, center $(512,512)$;
- pixel-center rasterization and nearest-triangle depth;
- flat per-face normals encoded as $(n+1)127.5$;
- perspective-correct reciprocal-depth interpolation;
- 11×11 SSIM windows with the official constants and foreground mask.

The renderer is too costly to call for every one of millions of collapses. It is instead used at checkpoints and for small offline neighborhoods. Three structural operators then recover quality that local edge contractions cannot:

### 5.1 Edge flips

For two triangles sharing an edge, the alternate diagonal is evaluated if it preserves orientation, manifoldness, and the distance certificate. The flip is kept only when the local six-view SSIM improves.

### 5.2 Fan deletion and retriangulation

Deleting a low-value vertex creates a polygonal one-ring. For small rings, all valid triangulations are enumerated through Catalan recursion; for larger rings, a restricted fan/dynamic-programming family is used. Each proposal is tested structurally before rendering.

### 5.3 Coordinate fitting

After topology is fixed, small tangent-space and normal-space perturbations optimize pixel-level SSIM. A trust region, triangle-area guard, and cluster-radius bound prevent geometric drift.

These operators are offline. Accepted edits are canonicalized and encoded for deterministic online replay.

## 6. Multi-fidelity solving workflow

```text
parse and normalize
  → build adjacency, quadrics, cluster statistics, visibility weights
  → guarded QEM collapse with versioned priority queue
  → checkpoint at target-count bands
  → fast topology + radius screen
  → symmetric vertex-set Hausdorff and specification-matching 1024² six-view SSIM
  → renderer-aware flips / fan retriangulation / coordinate fitting
  → canonicalize and bit-pack accepted edits
  → replay from checkpoints in the submitted program
  → final independent fail-closed validation
```

![Candidates become progressively more expensive to evaluate](../paper/figures/validation_funnel.png)

## 7. Runtime and source-budget engineering

- **Versioned heap entries** avoid global decrease-key operations; stale entries are discarded on pop.
- **Heap rebuilding** bounds memory after many invalidated candidates.
- **Static CSR adjacency plus append-only deltas** avoids rebuilding million-vertex neighborhoods.
- **Checkpointing** prevents replaying the full contraction history for every local experiment.
- **Canonical ordering** makes topology scripts stable under internal renumbering.
- **Bit-packed replay** stores vertex references, operator types, fan choices, and quantized coordinate deltas near their information-theoretic size.
- **Reference caches** remove redundant rendering work on the largest cases.
- **Fail-closed control flow** retains the last validated or previously Accepted checkpoint whenever time or validation margin becomes uncertain.

The final accepted source is 130,973 bytes, only 99 bytes below the 131,072-byte limit. Its readability is intentionally secondary to conformance; use the research core and this document to understand the method.

## 8. Why the hybrid dominates any one component

QEM provides throughput and geometric coherence. Cluster-normal memory stabilizes appearance over long collapse sequences. The specification-matching local renderer captures occlusion and SSIM effects that no collapse surrogate fully models. Structural edits escape the contraction-only topology family. Scoped validators reject proxy improvements that violate known constraints, and Kattis remains the hidden-test authority. Compression and replay convert expensive offline search into a contest-conforming program within the published limits.

This division of labor—not one parameter setting—is the central model.
