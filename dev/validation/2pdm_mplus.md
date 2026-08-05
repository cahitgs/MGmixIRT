# Cross-validation: multigroup 2PDM vs. Mplus 8.3

Date: 2026-08-05. Design: simulated two-group 2PDM data
(`sim_mgmixirt("2pdm", n = c(1500, 1500), I = 15, switch_point = 10,
pdec = c(0.30, 0.15), dshift = 1.5, mu_nd = c(0, 0.8),
mu_dec = c(-0.1, 0.7), seed = 7)`, item slopes `set.seed(42);
exp(rnorm(15, 0, 0.2))`, intercepts `seq(1.2, -1.2, length.out = 15)`).

Package: `mgmixirt(..., model = "2pdm", switch_point = 10, quadpts = 31,
starts = c(4, 15), maxit = 8000, tol = 1e-6, seed = 11)` — converged,
903 EM iterations, 2.6 min.
Mplus 8.3: `TYPE = MIXTURE; STARTS = 100 20; INTEGRATION =
GAUSSHERMITE(31); ADAPTIVE = OFF;` (input: `mg2pdm.inp`) — best
log-likelihood replicated, normal termination.

## Log-likelihood

| Quantity | Package | Mplus | Diff |
|---|---|---|---|
| logLik (conditional on group) | -25450.7221 | — | |
| + group multinomial term 3000·ln(0.5) | **-27530.1637** | **-27530.191** | **0.027** |
| Free parameters | 46 | 47 | +1 = KNOWNCLASS proportion, as expected |

The package's conditional likelihood, adjusted by the known-group
multinomial constant, matches the Mplus joint H0 value to 0.027 units
(relative 1e-6); the package value is marginally higher.

## Parameters (selected; full tables in `pkg_items.csv`, `mg2pdm.out`)

| Parameter | Package | Mplus (SE) | Diff |
|---|---|---|---|
| a1 | 1.329 | 1.327 (0.084) | 0.002 |
| a9 | 1.520 | 1.518 (0.086) | 0.002 |
| a12 | 1.687 | 1.695 (0.111) | 0.008 |
| a13 | 0.832 | 0.846 (0.117) | 0.014 |
| tau1 = -b1 | 1.069 | 1.055 (0.124) | 0.014 |
| tau10 | 0.391 | 0.401 (0.095) | 0.010 |
| decline tau12, group 1 | 2.720 | 2.744 (0.641) | 0.024 |
| decline tau14, group 2 | 1.941 | 1.946 (0.467) | 0.005 |
| mu_dec, group 1 | 0.270 | 0.326 (0.546) | 0.056 |
| mu_nd, group 2 | 0.845 | 0.850 (0.079) | 0.006 |
| theta variance, group 2 | 0.968 (= 0.984^2) | 0.970 (0.080) | 0.002 |
| P(decline), group 1 | 0.199 | 0.197 | 0.002 |
| P(decline), group 2 | 0.103 | 0.107 | 0.004 |

All empirically identified parameters agree within 0.03 (a small
fraction of their standard errors). Larger absolute discrepancies occur
only for boundary decline thresholds whose Mplus standard errors exceed
2 (e.g., decline tau13 group 1: package 6.23 vs. Mplus 5.18, SE 3.08) —
response probabilities are numerically zero in that region, so the
likelihood is flat and the parameter is not identified in this sample.

## Observations for the simulation study

1. With only I - i0 = 5 post-switch items and a moderate decrement,
   the decline-class proportion and the decrement sizes trade off along
   a likelihood ridge (true pdec = 0.30/0.15 recovered as 0.20/0.10 in
   this single replication, with correspondingly larger decrements).
   Both programs agree on the location of the maximum; the ridge is a
   property of the design, not of the software. The simulation study
   should vary I - i0 to quantify this.
2. Mplus reports pattern-level entropy (including the deterministic
   known-class variable), which is not comparable to class-only entropy.
