## Model specification and linear-predictor construction.
##
## Notation follows List et al. (2017). A person with switching point
## delta responds to items i <= delta under the regular 2PL regime
##   eta = a_i * theta + b_i          (slope/intercept form)
## and to items i > delta under the model-specific aberrant regime:
##   2pdm:   eta = a_i * theta + b_i - exp(d_ig)   (lower intercepts)
##   hybrid: eta = btilde_g                        (zero slope, common intercept)
##   mpdm:   eta = a_i * (theta - kappa_g * (I - delta)) + b_i
## Class-specific ability means follow mu_{ig} = mu_{Ig} + rho_g * (I - i)
## for the multiclass models (List et al. 2017, Eq. 6); the 2PDM has a
## free decline-class mean per group.

mgm_spec <- function(model, I, G, i0 = NULL) {
  model <- match.arg(model, c("2pdm", "hybrid", "mpdm", "2pl"))
  if (model == "2pdm") {
    if (is.null(i0)) stop("switch_point (i0) is required for model = \"2pdm\"")
    i0 <- as.integer(i0)
    if (i0 < 1L || i0 >= I) stop("switch_point must lie in 1, ..., I-1")
    deltas <- c(i0, I)
  } else if (model == "2pl") {
    deltas <- I                             # everyone in the no-decline class
  } else {
    deltas <- seq_len(I)
  }
  list(model = model, I = I, G = G, i0 = i0,
       deltas = deltas, nclass = length(deltas))
}

## Ability means by class for group g (vector over spec$deltas).
theta_means_g <- function(par, spec, g) {
  if (spec$model == "2pdm") {
    c(par$mu_dec[g], par$mu_nd[g])
  } else if (spec$model == "2pl") {
    par$mu_nd[g]
  } else {
    par$mu_nd[g] + par$rho[g] * (spec$I - spec$deltas)
  }
}

## Cell layout for one group: classes crossed with quadrature nodes,
## class-major ordering (for class c the Q nodes are contiguous).
cells_g <- function(par, spec, g, quad) {
  Q <- length(quad$z)
  delta_cell <- rep(spec$deltas, each = Q)
  z_cell <- rep(quad$z, times = spec$nclass)
  mu <- theta_means_g(par, spec, g)
  sigma <- exp(par$logsigma[g])
  theta_cell <- rep(mu, each = Q) + sigma * z_cell
  list(delta = delta_cell, theta = theta_cell, Q = Q)
}

## Linear predictor matrix (I x C) for group g given cell layout.
eta_g <- function(par, spec, g, cl) {
  I <- spec$I
  eta <- outer(par$a, cl$theta) + par$b
  ab <- outer(seq_len(I), cl$delta, ">")     # aberrant regime indicator
  if (!any(ab)) return(eta)
  if (spec$model == "2pdm") {
    shift <- matrix(0, I, length(cl$delta))
    rows <- (spec$i0 + 1L):I
    dec_cols <- cl$delta < I
    shift[rows, dec_cols] <- exp(par$d[, g])
    eta <- eta - shift * ab
  } else if (spec$model == "hybrid") {
    eta[ab] <- par$btilde[g]
  } else {                                    # mpdm
    kap <- exp(par$logkappa[g])
    shift <- outer(par$a, kap * (spec$I - cl$delta))
    eta[ab] <- eta[ab] - shift[ab]
  }
  eta
}

## Linear predictor for a single item (vector over cells of group g).
eta_item_g <- function(i, a_i, b_i, d_ig = NULL, par, spec, g, cl) {
  eta <- a_i * cl$theta + b_i
  ab <- cl$delta < i                          # aberrant cells for this item
  if (!any(ab)) return(eta)
  if (spec$model == "2pdm") {
    eta[ab] <- eta[ab] - exp(d_ig)
  } else if (spec$model == "hybrid") {
    eta[ab] <- par$btilde[g]
  } else {
    kap <- exp(par$logkappa[g])
    eta[ab] <- eta[ab] - a_i * kap * (spec$I - cl$delta[ab])
  }
  eta
}

## Number of free parameters (conditional on known group membership).
mgm_npar <- function(spec) {
  I <- spec$I; G <- spec$G
  base <- 2L * I + 2L * (G - 1L)              # a, b, mu_nd (non-ref), sigma (non-ref)
  switch(spec$model,
    "2pl"    = base,
    "2pdm"   = base + (I - spec$i0) * G + G + G,   # d, mu_dec, pdec
    "hybrid" = base + G + G + 2L * G,              # btilde, rho, tau1 + omega
    "mpdm"   = base + G + G + 2L * G)              # kappa, rho, tau1 + omega
}
