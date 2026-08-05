# Simulation study design (paper Section 4)

Two arms, G = 2 groups throughout, mu_nd = (0, 0.5), sigma = (1, 0.95).

## Arm A — parameter recovery under the true model

Factors (fully crossed, 24 cells, 100 replications each):

| Factor | Levels |
|---|---|
| Generating/fitted model | 2pdm, hybrid, mpdm |
| n per group | 1000, 3000 |
| Test length I | 20, 40 |
| Information regime | high: pi_nd = (.55, .70), omega = (2.5, 4) / low: pi_nd = (.80, .90), omega = (5, 8) |

Fixed structural values: rho = (-0.03, -0.06); HYBRID btilde = (-2.5, -2.5);
MPDM kappa = (0.5, 0.8); 2PDM i0 = 2I/3, dshift = 1.5, decline probability
= 1 - pi_nd of the regime. Item parameters per cell: a ~ exp(N(0, 0.2)),
b equally spaced in [1.4, -1.4] (drawn once per replication seed).

The "low" regime mirrors the empirical estimates of List et al. (2017)
(omega 4.8-10.3, pi_nd 0.62-0.91) and is expected to show the
class-size/magnitude likelihood ridge documented in dev/validation; the
"high" regime quantifies how much information is needed for usable
class-level estimates.

Outcomes: bias/RMSE per parameter (natural scale), cor(a_hat, a_true),
mean |b error|, convergence rate, EM iterations, wall time, fitted-vs-true
log-likelihood difference (positive values confirm the optimizer).

## Arm B — model selection by information criteria

I = 20 cells only (12 cells), 50 replications; each dataset fitted by all
four models (2pl, 2pdm, hybrid, mpdm); 2PDM fitted with i0 fixed at truth
for its own data and at 2I/3 otherwise. Outcomes: AIC/BIC/SABIC selection
rates (analogue of List et al.'s Table 1 with known truth), plus the
2PL-vs-mixture group-mean-difference distortion (their Table 3 effect:
d from the 2PL vs. d from the true model).

## Execution

`run_sim.R` — parallel PSOCK workers, one RDS checkpoint per
(arm, cell, rep, fitted model); safe to interrupt and resume (existing
checkpoints are skipped). Estimation settings: quadpts = 13,
starts = c(3, 10), maxit = 3000, tol = 1e-5 (I = 20) / 1e-4 (I = 40).
Seeds: rep-level seed = 10000 * cell_id + rep.
`collect_sim.R` — aggregates checkpoints into `results/summary_*.rds`.
