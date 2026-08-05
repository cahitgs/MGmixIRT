#' Simulate data from a multigroup mixture performance decline model
#'
#' Generates dichotomous item responses from the multigroup 2PDM, HYBRID,
#' or MPDM (List et al., 2017). Group sizes and group-specific structural
#' parameters are supplied as vectors with one element per group.
#'
#' @param model One of \code{"2pdm"}, \code{"hybrid"}, \code{"mpdm"}, or
#'   \code{"2pl"} (no decline classes).
#' @param n Integer vector of group sample sizes.
#' @param I Number of items.
#' @param a,b Item slopes and intercepts (length \code{I}). Defaults:
#'   slopes drawn from a log-normal distribution, intercepts from N(0, 1).
#' @param switch_point Switching point \code{i0} (2PDM only).
#' @param pdec Decline-class probability per group (2PDM only).
#' @param dshift Intercept decrement per group for items after the switching
#'   point (2PDM only); scalar or vector recycled over groups, or an
#'   \code{(I - i0) x G} matrix.
#' @param pi_nd No-decline class probability per group (HYBRID/MPDM).
#' @param omega Shape parameter per group (HYBRID/MPDM).
#' @param btilde Common aberrant-regime intercept per group (HYBRID).
#' @param kappa Decrement parameter per group (MPDM).
#' @param mu_nd No-decline class ability mean per group.
#' @param mu_dec Decline class ability mean per group (2PDM).
#' @param rho Slope of the class-mean function per group (HYBRID/MPDM).
#' @param sigma Ability standard deviation per group.
#' @param seed Optional RNG seed.
#'
#' @return A list with the response matrix \code{resp}, the group factor
#'   \code{group}, the true switching points \code{delta}, abilities
#'   \code{theta}, and the generating parameters \code{truth}.
#' @export
sim_mgmixirt <- function(model = c("2pdm", "hybrid", "mpdm", "2pl"),
                         n, I, a = NULL, b = NULL,
                         switch_point = NULL, pdec = NULL, dshift = 1,
                         pi_nd = NULL, omega = NULL,
                         btilde = NULL, kappa = NULL,
                         mu_nd = NULL, mu_dec = NULL, rho = NULL,
                         sigma = NULL, seed = NULL) {
  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)
  G <- length(n)
  if (is.null(a)) a <- exp(rnorm(I, 0, 0.25))
  if (is.null(b)) b <- rnorm(I, 0, 1)
  if (is.null(mu_nd)) mu_nd <- rep(0, G)
  if (is.null(sigma)) sigma <- rep(1, G)
  spec <- mgm_spec(model, I, G, i0 = switch_point)

  if (model == "2pdm") {
    if (is.null(pdec)) pdec <- rep(0.2, G)
    if (is.null(mu_dec)) mu_dec <- mu_nd
    if (!is.matrix(dshift))
      dshift <- matrix(dshift, I - spec$i0, G)
  } else if (model != "2pl") {
    if (is.null(pi_nd)) pi_nd <- rep(0.8, G)
    if (is.null(omega)) omega <- rep(5, G)
    if (is.null(rho)) rho <- rep(0, G)
    if (model == "hybrid" && is.null(btilde)) btilde <- rep(-2.2, G)
    if (model == "mpdm" && is.null(kappa)) kappa <- rep(0.5, G)
  }

  resp <- matrix(NA_integer_, sum(n), I)
  group <- factor(rep(seq_len(G), n))
  delta <- integer(sum(n))
  theta <- numeric(sum(n))
  row <- 0L
  for (g in seq_len(G)) {
    if (model == "2pdm") {
      pi_g <- c(pdec[g], 1 - pdec[g])
    } else if (model == "2pl") {
      pi_g <- 1
    } else {
      i <- seq_len(I - 1L)
      pi_dec <- (i^omega[g] - (i - 1)^omega[g]) / (I - 1)^omega[g] *
        (1 - pi_nd[g])
      pi_g <- c(pi_dec, pi_nd[g])
    }
    mu_g <- if (model == "2pdm") c(mu_dec[g], mu_nd[g])
            else if (model == "2pl") mu_nd[g]
            else mu_nd[g] + rho[g] * (I - spec$deltas)
    for (p in seq_len(n[g])) {
      row <- row + 1L
      cls <- sample.int(spec$nclass, 1L, prob = pi_g)
      dlt <- spec$deltas[cls]
      th <- rnorm(1, mu_g[cls], sigma[g])
      eta <- a * th + b
      abr <- seq_len(I) > dlt
      if (any(abr)) {
        if (model == "2pdm") {
          eta[abr] <- eta[abr] - dshift[which(abr) - spec$i0, g]
        } else if (model == "hybrid") {
          eta[abr] <- btilde[g]
        } else {
          eta[abr] <- a[abr] * (th - kappa[g] * (I - dlt)) + b[abr]
        }
      }
      resp[row, ] <- rbinom(I, 1L, plogis(eta))
      delta[row] <- dlt
      theta[row] <- th
    }
  }
  colnames(resp) <- sprintf("i%02d", seq_len(I))
  truth <- list(model = model, a = a, b = b, switch_point = spec$i0,
                pdec = pdec, dshift = if (model == "2pdm") dshift,
                pi_nd = pi_nd, omega = omega, btilde = btilde,
                kappa = kappa, mu_nd = mu_nd, mu_dec = mu_dec,
                rho = rho, sigma = sigma)
  list(resp = resp, group = group, delta = delta, theta = theta,
       truth = truth)
}
