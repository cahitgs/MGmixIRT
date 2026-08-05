# MGmixIRT (development version)

* Initial development version: multigroup 2PDM, HYBRID, and MPDM estimation by
  MML-EM with Gauss-Hermite quadrature; data simulation for all three models.
* Vectorised damped-Newton item updates in the M-step (about twice as fast as
  the initial per-item optimiser, identical maxima).
* `mplus_syntax()` generates the full Mplus input for any of the three models
  (List et al., 2017, supplement style), for cross-validation and migration.
* `scores()` returns person-level EAP abilities, posterior class
  probabilities, and modal switching points.
* Cross-validated against Mplus 8.3 on a two-group 2PDM: log-likelihoods
  agree to 0.03 units and all identified parameters to < 0.03
  (see `dev/validation/`).
