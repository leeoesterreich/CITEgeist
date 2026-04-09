# NB Likelihood Surface Geometry for the CITEgeist Cold-Start Problem

## Executive Summary

This note analyzes the Phase 1 negative-binomial (NB) deconvolution objective implemented in `CITEgeist/model/nb_deconvolution.py`, `CITEgeist/model/nb_likelihood.py`, and `CITEgeist/model/nb_initialization.py`, with emphasis on why QP-initialized proportions outperform cold-start NNLS initialization.

The main conclusion is:

1. The Phase 1 objective is genuinely nonconvex in the joint variables `(p, lambda, r, b, s)`, and different initializations can converge to different stationary points.
2. The largest source of basin sensitivity is not merely NB noise; it is the interaction of:
   - bilinear mean structure `p @ lambda`,
   - sparse / shared marker panel geometry,
   - softmax conditioning in logits,
   - and the moment-matched dispersion `\tilde r_{im}` that depends on both local and global parameters.
3. A good initialization supplies a better effective design matrix for learning `lambda`, so it does provide materially more Fisher information in finite samples.
4. Replacing the NB objective by a simpler Gaussian or Poisson-like approximation would reduce some nonconvexity, but it would not remove the fundamental basin issue because the joint factorization in `(p, lambda)` is still nonconvex.
5. Therefore, the observed gap
   `0.726 - 0.640 = 0.086`
   is best interpreted as a local-attractor problem, not merely a convergence-speed problem.

## 1. Implemented Objective

For spot `i = 1, ..., I`, marker `m = 1, ..., M`, and cell type `t = 1, ..., T`,

```math
S_{im} \sim \operatorname{NB}(\mu_{im}, \tilde r_{im}),
\qquad
\mu_{im} = s_i \left( N_i \sum_{t=1}^T p_{it}\lambda_{tm} + b_m \right).
```

The implementation parameterizes proportions by logits `eta_i \in \mathbb R^T`:

```math
p_{it} = \frac{\exp(\eta_{it}/\tau)}{\sum_u \exp(\eta_{iu}/\tau)}.
```

In the current code path, the approximate spot-level variance is

```math
\operatorname{Var}(S_{im})
= s_i^2 \Bigg[
N_i \sum_t p_{it}\left(\lambda_{tm} + \frac{\lambda_{tm}^2}{r_m}\right)
+
N_i \left(\sum_t p_{it}\lambda_{tm}^2 - \left(\sum_t p_{it}\lambda_{tm}\right)^2\right)
+
\left(b_m + \frac{b_m^2}{r_m}\right)
\Bigg].
```

Then the code moment-matches a spot-marker-specific NB dispersion

```math
\tilde r_{im}
=
\frac{\mu_{im}^2}{\operatorname{Var}(S_{im})-\mu_{im}},
```

with two implementation details that matter geometrically:

1. `Var - mu` is floored below by `0.01 * mu`.
2. `\tilde r_{im}` is capped above at `100`.

Thus the objective used in training is not a standard NB likelihood with a fixed dispersion parameter. It is a profile likelihood through a nonlinear map

```math
(p,\lambda,r,b,s) \mapsto (\mu,\tilde r(\mu,\operatorname{Var})).
```

The optimization routine alternates:

1. Global step: optimize `(lambda, r, b)` with `(p, s)` fixed.
2. Local step: optimize `(eta, s)` with `(lambda, r, b)` fixed.

using finite-step Adam in each block, not exact block maximization.

## 2. Marker-Panel Geometry in the Current Implementation

The current curated panel in `nb_initialization.py` has:

- `T = 7` cell types
- `M = 18` markers
- `27` active type-marker pairs out of `7 x 18 = 126`
- active-mask density `27 / 126 = 0.214`, i.e. `21.4%`

So the implemented panel is even sparser than the rough `~30%` description.

The marker structure is highly uneven:

- 12 markers are exclusive to one type
- 4 markers are shared by exactly 2 types
- 2 markers are shared by 3 or more types

Type-specific exclusive markers are:

- B cells: `CD20`, `CD45RA`, `CD138`
- CD4 T: `CD4`
- CD8 T: `CD8A`
- Macrophages: `CD163`, `CD16`, `CD11c`, `E-Cadherin`
- Endothelial: `CD31`
- Epithelial: `PanCK`
- Fibroblasts: `alphaSMA`

Shared markers create the main confounding directions:

- `CD3E`, `CD45RO`: CD4 vs CD8
- `CD68`, `HLA-DR`, `CD45`: immune cross-talk
- `Vimentin`: macrophage vs endothelial vs fibroblast

This already suggests the expected basin geometry:

- sharp curvature near types with several exclusive markers,
- broad ridges along partially exchangeable type pairs,
- and especially flat valleys when one type has only one exclusive anchor.

## 3. Multimodality in Proportion Space

### 3.1 Fixed globals: what is identifiable?

For fixed `(lambda, r, b, s)`, one spot contributes

```math
\ell_i(p_i) = \sum_{m=1}^M \log f_{\mathrm{NB}}\!\left(S_{im}; \mu_{im}(p_i), \tilde r_{im}(p_i)\right),
\qquad
p_i \in \Delta^{T-1}.
```

The mean map is affine in `p_i`:

```math
\mu_i(p_i) = s_i \left(N_i \Lambda^\top p_i + b\right),
```

where `\Lambda = (\lambda_{tm})_{t,m}`.

Only directions in the simplex tangent space

```math
\mathcal T = \{v \in \mathbb R^T : \mathbf 1^\top v = 0\}
```

are estimable. If two types have nearly identical active-marker profiles, then

```math
\Lambda^\top v \approx 0
```

for some nonzero `v \in \mathcal T`, creating a near-flat ridge in `\ell_i`.

The local rank is bounded by

```math
\operatorname{rank}(D\mu_i|_{\mathcal T})
\le \min(T-1,\operatorname{rank}(\Lambda)).
```

Since `T-1 = 6`, full local identifiability requires six well-separated contrast directions. The panel is likely full rank generically, but not uniformly well conditioned. The weak directions are biologically obvious:

- CD4 vs CD8 when `CD4` and `CD8A` are low or noisy.
- Endothelial vs fibroblast vs macrophage through `Vimentin`.
- Immune mixtures when `CD45` dominates but subtype markers are weak.

### 3.2 Is the likelihood multimodal?

Yes, in the practically relevant sense.

For fixed globals and fixed `\tilde r`, the objective would be much better behaved. For example:

- Poisson log-likelihood in `p` is concave because `k\log(a^\top p + c) - (a^\top p + c)` is concave on the simplex.
- Homoscedastic Gaussian least squares is a convex quadratic in `p`.

But the implemented objective is not in that regime because `\tilde r_{im}` depends nonlinearly on `p` through both `\mu` and `\operatorname{Var}`. That adds nonconcave curvature. Combined with shared-marker ridges, this creates multiple stationary points in proportion space.

For this panel, the most plausible basin picture is:

1. Vertex-adjacent basins for types with multiple exclusive markers.
2. Edge/ridge basins for partially exchangeable pairs such as CD4/CD8.
3. Interior shallow basins where shared markers can be explained by mixed immune or stromal states.

The number of strict local maxima need not be enormous, but the number of practically distinct attraction regions can be. In other words, the surface is likely ridge-dominated rather than combinatorially explosive.

### 3.3 Effect of softmax parameterization

Softmax does not remove nonconvexity. It changes the geometry.

Let

```math
p = \operatorname{softmax}(\eta/\tau), \qquad
J_\eta = \frac{\partial p}{\partial \eta}
= \frac{1}{\tau}\left[\operatorname{diag}(p)-pp^\top\right].
```

Consequences:

1. Interior simplex points correspond one-to-one to equivalence classes of logits `eta + c 1`.
2. Therefore, interior local optima in `p` correspond to interior stationary manifolds in `eta`; softmax does not eliminate them.
3. Boundary optima in constrained `p` become `\|\eta\| \to \infty` limits rather than finite points.
4. Conditioning changes substantially because gradients in logits are

```math
\nabla_\eta \ell = J_\eta^\top \nabla_p \ell.
```

If a true type is initialized with `p_{it} \approx 0`, then its correction signal is multiplied by approximately `p_{it}`. Thus softmax can make basin escape slow once a type is nearly suppressed.

That last point matters for cold start. NNLS initializations that over-collapse a spot toward the wrong face of the simplex can become hard to repair even if the optimum is not far away in `p` distance.

## 4. Fisher Information for Learning `lambda`

The central question is whether better initial proportions provide more gradient signal to learn global parameters. The answer is yes.

### 4.1 Leading-order Fisher information

For an NB observation with fixed dispersion `\tilde r`, the score for the mean is

```math
\frac{\partial \ell}{\partial \mu}
= \frac{\tilde r (S-\mu)}{\mu(\tilde r+\mu)},
```

and the Fisher information for `\mu` is

```math
\mathcal I_\mu
=
\mathbb E\left[\left(\frac{\partial \ell}{\partial \mu}\right)^2\right]
= \frac{\tilde r}{\mu(\tilde r+\mu)}.
```

Ignoring for the moment the derivative of `\tilde r` with respect to parameters, the information for a global emission parameter `\lambda_{tm}` is

```math
\mathcal I(\lambda_{tm})

\approx
\sum_{i=1}^I
\frac{\tilde r_{im}}{\mu_{im}(\tilde r_{im}+\mu_{im})}
\left(\frac{\partial \mu_{im}}{\partial \lambda_{tm}}\right)^2

=
\sum_{i=1}^I
w_{im}\,(s_i N_i p_{it})^2,
```

where

```math
w_{im} := \frac{\tilde r_{im}}{\mu_{im}(\tilde r_{im}+\mu_{im})}.
```

For all types jointly, the information matrix for marker `m` has leading term

```math
\mathcal I_m(\lambda_{\cdot m})
\approx
\sum_{i=1}^I
w_{im}(s_i N_i)^2\, p_i p_i^\top.
```

This is the key result. The information available to learn `\lambda_{\cdot m}` is controlled by the weighted Gram matrix of the current proportions.

### 4.2 Why QP initialization is more informative

If `p_i` are close to ground truth, then:

1. Each type gets exposed on the correct subset of spots.
2. The columns of `P = (p_{it})` are less aliased.
3. The Gram matrix `\sum_i w_{im}(s_i N_i)^2 p_i p_i^\top` has larger informative eigenvalues in subtype-contrast directions.

If NNLS proportions are rough and mix confusable types, then two failures occur:

1. **Attenuation**: true type-specific spots contribute diluted weights to the correct `\lambda_{tm}`.
2. **Collinearity**: confusable types acquire similar columns in `P`, shrinking contrast eigenvalues and making `\lambda` poorly estimable.

Thus QP initialization does not create new information ex nihilo, but in finite sample it does create a materially more informative effective design for the global update. In that sense it provides fundamentally more usable information for learning `\lambda`.

This is especially important in CITEgeist because the global step learns `\lambda` from the current `p`, and the local step cannot recover information that never entered the global block.

### 4.3 What about the `\tilde r` dependence?

The exact Fisher information also contains terms involving

```math
\frac{\partial \tilde r_{im}}{\partial \lambda_{tm}},
```

because `\tilde r_{im}` is itself a function of `\mu` and `\operatorname{Var}`. Those terms increase curvature heterogeneity further, but they do not change the qualitative conclusion above: better proportions improve the conditioning of the global `\lambda` update.

## 5. Can Block Alternation Reach Different Fixed Points?

Yes. There is no guarantee in the current algorithm that all initializations converge to the same fixed point.

### 5.1 Why the standard EM guarantees do not apply directly

Classical EM theory gives monotone likelihood ascent and convergence to a stationary point under exact E/M updates and regularity conditions [1]. ECM extends this to exact conditional-maximization blocks [2].

Those guarantees do not directly cover the present algorithm because:

1. The model is not written as a standard latent-variable complete-data EM with exact conditional expectations.
2. Each block uses a finite number of Adam steps, not an exact maximizer.
3. The global block is trained only on the training subset, while the local block uses all spots.
4. Early stopping is based on held-out likelihood, so the actual stopping point may be pre-asymptotic.

This is closer to inexact nonconvex block coordinate descent than to exact EM.

### 5.2 What can be said theoretically?

For nonconvex block coordinate methods, the generic guarantee is only that limit points are stationary under additional regularity assumptions, often requiring exact or sufficiently accurate block solves and sometimes uniqueness of block minimizers [3,4]. Those assumptions are not strong enough here to imply a unique global fixed point.

Therefore two statements are both true:

1. The algorithm may converge stably.
2. Different initializations may converge to different stationary points.

The second statement is exactly the one relevant to the cold-start gap.

### 5.3 Why different fixed points are plausible here

The joint objective has several coupled symmetries / near-symmetries:

- exchangeability along shared-marker directions,
- scale tradeoffs among `s_i`, `lambda`, and `b`,
- and overdispersion tradeoffs through `r` and `\tilde r`.

Suppose initialization A places spots near the correct immune-subtype faces of the simplex. Then the global step learns subtype-separated `\lambda`. Once that happens, the local step sees stronger subtype gradients, reinforcing the good basin.

Suppose initialization B mixes those same spots across CD4/CD8 or stromal types. Then the global step learns more averaged `\lambda`; afterward the local block sees a flatter surface and remains in the mixed basin.

That is a textbook self-reinforcing local-attractor mechanism.

## 6. Extra Nonconvexity from Moment-Matched Dispersion

### 6.1 Why `\tilde r` coupling matters

In a simpler mean-only model, the marker likelihood depends on `p_i` only through the affine map `\mu_i(p_i)`. Here the likelihood also depends on

```math
\tilde r_{im}(p,\lambda,r,b,s)
=
\frac{\mu_{im}^2}{\operatorname{Var}_{im}(p,\lambda,r,b,s)-\mu_{im}}.
```

The variance contains:

- a within-type term linear in `p`,
- a between-type term
  `\sum_t p_t \lambda_{tm}^2 - (\sum_t p_t\lambda_{tm})^2`,
  which is concave in `p`,
- and inverse-dispersion terms `\lambda_{tm}^2/r_m`, `b_m^2/r_m`.

Hence `\tilde r` contains ratios of quadratic and affine forms in `p`. This is a direct source of additional nonconvex curvature.

### 6.2 Basin consequences

This coupling can create effects absent in Gaussian least squares:

1. A mixed proportion vector can be preferred because it inflates between-type variance and therefore explains noisy counts via a smaller `\tilde r`.
2. The optimizer can trade off fit in the mean against fit in the variance.
3. Wrong early `\lambda` estimates alter the variance surface seen by the local step, deepening the wrong basin.

So yes: the moment-matched dispersion can create additional local optima, or at minimum deepen and widen existing attraction regions.

### 6.3 Would a Gaussian simplification eliminate the gap?

Not in general.

If one replaced the likelihood by a homoscedastic Gaussian loss with fixed variance, then:

- the local `p` subproblem with fixed globals becomes much closer to convex,
- and one important source of nonconvexity disappears.

But the joint problem in `(p,\lambda)` remains bilinear:

```math
\mu_{im} = s_i\left(N_i \sum_t p_{it}\lambda_{tm}+b_m\right),
```

so block alternation can still have multiple stationary points.

Therefore:

1. Gaussianization would likely reduce the size of the gap.
2. It would not guarantee that the gap vanishes.
3. The main unresolved difficulty would still be initialization into the correct `(p,\lambda)` basin.

Strictly speaking, the limit `\tilde r \to \infty` in the NB parameterization is closer to a Poisson-like mean-only regime than to a standard Gaussian regime. That regime is more benign than the current model, but still jointly nonconvex once `\lambda` is unknown.

## 7. What Theory Actually Guarantees

### 7.1 EM / ECM theory

- Wu (1983) shows that EM converges to the set of stationary points under standard regularity assumptions; this does not imply convergence to the same stationary point from every initialization [1].
- Meng and Rubin (1993) show that ECM inherits EM-style monotone ascent when each conditional maximization is carried out exactly [2].

These are local-stationarity results, not basin-uniqueness results.

### 7.2 Local convergence rate depends on overlap

Xu and Jordan (1996) analyze EM for Gaussian mixtures and show that the local behavior depends strongly on overlap between components [5]. Ma, Xu, and Jordan (2000) further show faster local convergence when overlap is small [6].

The closest analog here is marker-profile overlap between cell types:

- well-separated types with several exclusive markers behave like low-overlap components,
- while CD4/CD8 and stromal shared-marker directions behave like high-overlap components.

So the theory predicts exactly what the panel structure suggests: local convergence is fast inside a good basin, but the basin itself can be small when component overlap is high.

### 7.3 Basin-of-attraction analysis

There is no general theorem giving a useful global basin-of-attraction characterization for this NB deconvolution objective. The literature provides:

- stationary-point guarantees under regularity,
- local contraction results under strong separation assumptions,
- and sensitivity-to-initialization observations for mixture models more broadly.

That is enough to reject the hypothesis that all reasonable initializations must converge to the same Phase 1 solution.

## 8. Practical Interpretation of the 0.086 Gap

The empirical facts are:

- QP init: start at `r = 0.639`, finish at `r = 0.726`
- NNLS init: start at `r = 0.318`, finish at `r = 0.640`
- residual gap after Phase 1: `0.086`

The theoretically coherent interpretation is:

1. QP initialization begins inside, or close to, a basin where subtype-separating `\lambda` are learnable.
2. NNLS initialization begins outside that basin for a nontrivial fraction of spots.
3. Once Phase 1 learns slightly averaged globals from the cold start, the local step no longer has enough curvature to migrate to the QP-warm basin.

So the gap is not best viewed as "NB can get there eventually if run longer." Because the algorithm is alternating on a nonconvex surface, more iterations can just stabilize the wrong fixed point.

## 9. What Initialization Class Should Be Sufficient?

The goal is not to match QP numerically. It is to enter the same attraction region.

The theory above suggests that a sufficient class of initializations should satisfy:

1. **Support correctness**: true cell types for a spot should not be initialized at nearly zero probability, because softmax then suppresses the corrective gradient.
2. **Subtype separation on exclusive markers**: spots with clear `CD4`, `CD8A`, `CD31`, `PanCK`, `alphaSMA`, `CD163/CD16/CD11c` signal should be placed near the corresponding simplex face.
3. **Low aliasing across shared-marker pairs**: the initialization should avoid collapsing CD4/CD8 and endothelial/fibro/macrophage into broad mixed states when exclusive anchors are present.
4. **Reasonable global exposure matrix**: across spots, the matrix `P` should have enough column diversity that the Fisher information for `\lambda` is well conditioned.

This points to a concrete initialization class:

- not arbitrary NNLS,
- but constrained or regularized initializations that use the curated active mask and exclusive markers to preserve support and subtype anchors.

Examples include:

- exclusive-marker-first initialization,
- marker-group gating before continuous proportions,
- entropic/simplex smoothing that prevents exact zeros,
- temperature-raised softmax initialization,
- or multi-start schemes centered on several biologically plausible supports.

Theoretically, any initializer that places most spots inside the correct low-overlap neighborhood should be sufficient to reach the QP-warm basin. It does not need QP itself; it needs the same qualitative properties:

- preserve candidate support,
- separate anchor-defined types,
- and avoid early over-collapse.

## 10. Bottom Line

For the current CITEgeist Phase 1 NB solver, the cold-start problem is primarily a basin-of-attraction problem.

The sparse marker panel gives only a few strong subtype-separating directions. QP initialization uses those directions well enough to create a high-information global update for `\lambda`, after which block alternation reinforces the correct basin. NNLS initialization does not. Because the objective is nonconvex and the implemented optimizer is an inexact alternating Adam scheme rather than exact EM/ECM, there is no theoretical reason to expect both initializations to converge to the same fixed point.

Therefore, closing the `0.086` gap without QP is theoretically plausible, but only with an initialization strategy that is basin-aware. A generic unconstrained cold start is not expected to be sufficient.

## References

1. Wu, C. F. J. (1983). On the convergence properties of the EM algorithm. *The Annals of Statistics*, 11(1), 95-103. DOI: `10.1214/aos/1176346060`. URL: http://projecteuclid.org/euclid.aos/1176346060
2. Meng, X.-L., & Rubin, D. B. (1993). Maximum likelihood estimation via the ECM algorithm: A general framework. *Biometrika*, 80(2), 267-278. DOI: `10.1093/biomet/80.2.267`. URL: https://doi.org/10.1093/biomet/80.2.267
3. Grippo, L., & Sciandrone, M. (2000). On the convergence of the block nonlinear Gauss-Seidel method under convex constraints. *Operations Research Letters*, 26(3), 127-136. DOI: `10.1016/S0167-6377(99)00074-7`. URL: https://doi.org/10.1016/S0167-6377(99)00074-7
4. Razaviyayn, M., Hong, M., & Luo, Z.-Q. (2013). A unified convergence analysis of block successive minimization methods for nonsmooth optimization. *SIAM Journal on Optimization*, 23(2), 1126-1153. DOI: `10.1137/120891009`. URL: https://doi.org/10.1137/120891009
5. Xu, L., & Jordan, M. I. (1996). On convergence properties of the EM algorithm for Gaussian mixtures. *Neural Computation*, 8(1), 129-151. DOI: `10.1162/neco.1996.8.1.129`. URL: https://doi.org/10.1162/neco.1996.8.1.129
6. Ma, J., Xu, L., & Jordan, M. I. (2000). Asymptotic convergence rate of the EM algorithm for Gaussian mixtures. *Neural Computation*, 12(12), 2881-2907. DOI: `10.1162/089976600300014764`. URL: https://doi.org/10.1162/089976600300014764
7. McLachlan, G., & Peel, D. (2000). *Finite Mixture Models*. Wiley. URL: https://books.google.com/books/about/Finite_Mixture_Models.html?id=YXqflwEACAAJ
