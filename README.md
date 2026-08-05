# MGmixIRT

Multigroup mixture item response theory models for performance decline in
low-stakes assessments.

`MGmixIRT` provides the first ready-made implementation of the multigroup
mixture IRT models for performance decline (PD) described by List, Robitzsch,
Lüdtke, Köller, and Nagy (2017):

- **2PDM** — the two-class mixture 2PL model of Bolt, Cohen, and Wollack (2002),
- **HYBRID** — the HYBRID model of Yamamoto (1995),
- **MPDM** — the multiclass performance decline model of Jin and Wang (2014),

each extended to an arbitrary number of manifest groups. Models are estimated
by marginal maximum likelihood (EM algorithm with Gauss-Hermite quadrature).
Class probabilities for the multiclass models follow the parsimonious
shape-parameter formulation of Cao and Stokes (2008).

## Why this package

Until now these models could only be fitted with extensive hand-written Mplus
syntax (the online supplement of List et al. spans about 28 pages of code) or
with licence-restricted software. General-purpose IRT packages either drop the
manifest grouping when a mixture distribution is requested or cannot express
the switching-point structure and the nonlinear constraints on the class
probabilities.

## Installation

```r
# install.packages("remotes")
remotes::install_github("cahitgs/MGmixIRT")
```

## Quick start

```r
library(MGmixIRT)

## simulate two-group data with performance decline (HYBRID mechanism)
sim <- sim_mgmixirt(model = "hybrid", n = c(1500, 1500), I = 20,
                    pi_nd = c(0.65, 0.85), omega = c(5, 8),
                    mu_nd = c(-1, 0), rho = c(-0.01, -0.05), seed = 1)

fit <- mgmixirt(sim$resp, group = sim$group, model = "hybrid")
summary(fit)
```

## References

- Bolt, D. M., Cohen, A. S., & Wollack, J. A. (2002). Item parameter estimation
  under conditions of test speededness. *Journal of Educational Measurement,
  39*, 331–348.
- Cao, J., & Stokes, S. L. (2008). Bayesian IRT guessing models for partial
  guessing behaviors. *Psychometrika, 73*, 209–230.
- Jin, K.-Y., & Wang, W.-C. (2014). Item response theory models for performance
  decline during testing. *Journal of Educational Measurement, 51*, 178–200.
- List, M. K., Robitzsch, A., Lüdtke, O., Köller, O., & Nagy, G. (2017).
  Performance decline in low-stakes educational assessments: Different mixture
  modeling approaches. *Large-scale Assessments in Education, 5*, 15.
- Yamamoto, K. (1995). *Estimating the effects of test length and test time on
  parameter estimation using the HYBRID model* (TOEFL Technical Report TR-10).
  Educational Testing Service.
