## Latent class probability machinery.
##
## For the multiclass models (HYBRID, MPDM) the decline-class probabilities
## follow the shape-parameter formulation of Cao and Stokes (2008), as used
## by Jin and Wang (2014) and List et al. (2017, Eq. 3):
##   pi_i = [i^omega - (i-1)^omega] / (I-1)^omega * (1 - pi_I),  i < I,
## where pi_I is the probability of the no-decline class (delta = I).
## In multinomial-logit form with the no-decline class as reference the
## class logits are tau_i = tau_1 + log(i^omega - (i-1)^omega)
## (List et al. 2017, supplement, Eq. 6), which is the parameterisation
## estimated here: unconstrained (tau_1, log omega) per group.

## Class logits for one group. Classes are ordered as in the model spec:
## for "2pdm" c(decline, no-decline); otherwise delta = 1, ..., I with the
## no-decline class (delta = I) last. The reference (last) logit is 0.
class_logits_g <- function(model, I, tau1 = NULL, logomega = NULL,
                           logit_pdec = NULL) {
  if (model == "2pdm") {
    c(logit_pdec, 0)
  } else {
    i <- seq_len(I - 1L)
    c(tau1 + log_ipow_diff(i, exp(logomega)), 0)
  }
}

probs_from_logits <- function(tau) {
  e <- exp(tau - max(tau))
  e / sum(e)
}

## Solve tau_1 from a target no-decline probability and shape parameter,
## used for start values: pi_1/pi_I = (1 - pi_I)/pi_I / (I-1)^omega.
tau1_from_target <- function(pi_nd, omega, I) {
  log((1 - pi_nd) / pi_nd) - omega * log(I - 1)
}
