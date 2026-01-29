# Adaptive Per-Marker Baseline Design

## Problem

VIM (Vimentin) has high signal everywhere in tissue. After CLR + max-norm preprocessing, VIM values are in [0.7, 1.0] range. The current reconstruction model `S[i,m] = beta[m] * Y[i, owner(m)]` has no concept of baseline — it must assign Fibroblast proportion everywhere to reconstruct VIM's ubiquitous high signal. This causes Fibroblast over-prediction (bias = +0.205) and destroys spatial pattern accuracy (Pearson r = 0.156).

## Solution

Add a per-marker intercept (baseline) to the reconstruction model:

```
S[i,m] = alpha[m] + beta[m] * Y[i, owner(m)]
```

- `alpha[m]` captures the "always-on" component of marker m
- `beta[m]` captures the cell-type-dependent variation
- For specific markers (CD68, PanCK): alpha ≈ 0 naturally
- For ubiquitous markers (VIM): alpha learns the floor, beta explains peaks

## EM Integration

### E-step (Gurobi QP, fix alpha/beta, solve for Y)

Loss: `((S[i,m] - alpha[m]) - beta[m] * Y[i,j])^2`

Same QP structure as before, just with baseline-subtracted signal. No solver changes needed.

### M-step (fix Y, update alpha and beta)

For each marker m with combined owner proportions `Z[i] = sum(Y[i, owners])`:

```
beta[m] = Cov(S_m, Z) / Var(Z)
alpha[m] = mean(S_m) - beta[m] * mean(Z)
```

Then clip: `beta[m]` to [beta_min, beta_max], `alpha[m]` to [0, alpha_max].

Add L2 regularization on alpha: `lambda_alpha * alpha[m]^2` pushes alpha toward zero unless data strongly supports a nonzero floor.

With regularization, the M-step for alpha becomes:
```
alpha[m] = (mean(S_m) - beta[m] * mean(Z)) / (1 + lambda_alpha / N)
```

## Pipeline Placement

- **Global EM**: Alpha learned alongside beta in M-step, used in E-step.
- **Local finetuning**: Alpha passed as fixed input (not re-learned per neighborhood). Baseline is a global marker property, not local.

## Parameters

- `alpha_max: float = 0.8` — ceiling for learned baseline
- `lambda_alpha: float = 1.0` — L2 regularization strength toward zero

## Guardrails

- `alpha[m] >= 0` (no negative baseline)
- `alpha[m] <= alpha_max` (don't absorb all signal)
- Beta range [0.1, 2.0] unchanged
- Beta normalization unchanged, applied to residual signal
- Compatible with marker exclusivity weighting
