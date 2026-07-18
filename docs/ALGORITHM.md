# Algorithm and mathematical model

This chapter is the concise, web-readable counterpart of the [full Round 2 article](../paper/IMC_Challenge_Round2_NEU_AddictedTribes.pdf). It separates the general method from the case-specific replay data used by the final submission.

## 1. Optimization problem

For each input triangular mesh $M=(V,F)$, produce a closed triangular 2-manifold $M'=(V',F')$ that minimizes $|V'|$ subject to:

- valid indices, non-degenerate triangles, and exactly two incident faces per edge;
- symmetric vertex-set Hausdorff distance

  $$d_H(M,M') \le 0.05\,\operatorname{diag}(\operatorname{AABB}(M));$$

- perceptual score at least 0.9 under the official six-view normal/depth renderer.

For view $c$, the perceptual score is

$$S_c=\tfrac12\operatorname{SSIM}(I^N_c,I^{N'}_c)+\tfrac12\operatorname{SSIM}(I^D_c,I^{D'}_c),$$

and the test score is the mean of $S_c$ over the six axis-aligned cameras. Among valid outputs, ranking depends only on retained vertex ratios:

$$\operatorname{Score}=100\left(1-\frac16\sum_{i=1}^6\frac{|V'_i|}{|V_i|}\right).$$

The objective is discontinuous: topology and visibility changes can modify many image pixels at once, while the final score rewards each removed vertex equally. We therefore use a hierarchy of increasingly expensive models rather than attempting one monolithic optimizer.

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

$A_f^{\mathrm{proj}}$ approximates the total projected area over the six official cameras. Faces that occupy more evaluated pixels therefore receive more protection than equally sized but hidden or edge-on faces.

The general collapse objective is

$$E(x)=E_Q(x)+\lambda_N E_N(x)+\lambda_C E_C(x)+\lambda_R E_R(x),$$

where $E_C$ is cluster-normal memory and $E_R$ denotes curvature/visibility regularizers enabled only in selected regimes.

## 4. Cluster-normal memory

Repeated local normal comparisons gradually forget the original surface. We preserve that history explicitly. Initially, each vertex receives the additive area-weighted normals of its incident original faces:

$$C_v=\sum_{f\ni v} A_f n_f.$$

When $u$ contracts into $v$, the statistic is merged exactly:

$$C_v\leftarrow C_u+C_v.$$

For a proposed target $x$, affected triangles yield a proposed aggregate direction $\widehat C'$. The cluster term measures the loss relative to the stored original aggregate direction $\widehat C=(C_u+C_v)/\|C_u+C_v\|$:

$$E_C(x)=W_{uv}\left[\left(1-\widehat C\cdot\widehat C'\right)^q-\left(1-\widehat C\cdot\widehat C_{\mathrm{current}}\right)^q\right].$$

The subtraction makes this an incremental loss. The key property is associativity: $C_u+C_v$ is independent of contraction order, so the reference evidence does not drift. This change was decisive on the Bunny-like hidden case and became part of the readable research core in [`src/research/compact_qem_lab.cpp`](../src/research/compact_qem_lab.cpp).

![Original face evidence is accumulated through the collapse tree](../paper/figures/cluster_normal.png)

## 5. Exact renderer as an optimization oracle

The local evaluator mirrors the official specification:

- cameras at the positive/negative Cartesian axes, distance 2.5;
- 1024×1024 resolution, focal length 800, center $(512,512)$;
- pixel-center rasterization and nearest-triangle depth;
- flat per-face normals encoded as $(n+1)127.5$;
- perspective-correct reciprocal-depth interpolation;
- 11×11 SSIM windows with the official constants and foreground mask.

The renderer is too costly to call for every one of millions of collapses. It is instead used at checkpoints and for small offline neighborhoods. Three structural operators then recover quality that local edge contractions cannot:

### 5.1 Edge flips

For two triangles sharing an edge, the alternate diagonal is evaluated if it preserves orientation, manifoldness, and the distance certificate. The flip is kept only when exact six-view SSIM improves.

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
  → exact Hausdorff and 1024² six-view SSIM
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
- **Fail-closed control flow** retains the last certified checkpoint whenever time or validation margin becomes uncertain.

The final accepted source is 130,973 bytes, only 99 bytes below the 131,072-byte limit. Its readability is intentionally secondary to conformance; use the research core and this document to understand the method.

## 8. Why the hybrid dominates any one component

QEM provides throughput and geometric coherence. Cluster-normal memory stabilizes appearance over long collapse sequences. The exact renderer captures occlusion and SSIM effects that no local metric fully models. Structural edits escape the contraction-only topology family. Certification prevents improvements on a proxy or one view from becoming invalid submissions. Compression and replay convert the expensive offline search into a legal contest program.

This division of labor—not one parameter setting—is the central model.
