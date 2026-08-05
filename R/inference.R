## Standard errors and group comparisons.
##
## The observed information matrix is computed numerically on the
## estimation (transformed) scale by differentiating the marginal
## log-likelihood; standard errors for natural-scale group parameters
## and Wald tests of group equality follow by the delta method, as in
## the group comparisons of List et al. (2017, Table 2).

## Pack the free parameters into a named vector (estimation scale).
pack_par <- function(par, spec) {
  I <- spec$I; G <- spec$G
  x <- c(setNames(par$a, paste0("a", seq_len(I))),
         setNames(par$b, paste0("b", seq_len(I))))
  gl <- seq_len(G)
  if (spec$model == "2pl") {
    if (G > 1L)
      x <- c(x, setNames(par$mu_nd[-1L], paste0("mu_nd_g", gl[-1L])),
             setNames(par$logsigma[-1L], paste0("logsigma_g", gl[-1L])))
  } else if (spec$model == "2pdm") {
    dn <- as.vector(outer((spec$i0 + 1L):I, gl,
                          function(i, g) paste0("d", i, "_g", g)))
    x <- c(x, setNames(as.vector(par$d), dn),
           setNames(par$mu_dec, paste0("mu_dec_g", gl)),
           if (G > 1L) setNames(par$mu_nd[-1L], paste0("mu_nd_g", gl[-1L])),
           if (G > 1L) setNames(par$logsigma[-1L],
                                paste0("logsigma_g", gl[-1L])),
           setNames(par$logit_pdec, paste0("logit_pdec_g", gl)))
  } else {
    x <- c(x,
           if (spec$model == "hybrid")
             setNames(par$btilde, paste0("btilde_g", gl))
           else setNames(par$logkappa, paste0("logkappa_g", gl)),
           setNames(par$rho, paste0("rho_g", gl)),
           if (G > 1L) setNames(par$mu_nd[-1L], paste0("mu_nd_g", gl[-1L])),
           if (G > 1L) setNames(par$logsigma[-1L],
                                paste0("logsigma_g", gl[-1L])),
           setNames(par$tau1, paste0("tau1_g", gl)),
           setNames(par$logomega, paste0("logomega_g", gl)))
  }
  x
}

unpack_par <- function(x, spec, template) {
  I <- spec$I; G <- spec$G
  par <- template
  par$a <- unname(x[paste0("a", seq_len(I))])
  par$b <- unname(x[paste0("b", seq_len(I))])
  gl <- seq_len(G)
  if (spec$model == "2pl") {
    ## nothing beyond items and the shared tail handled below
  } else if (spec$model == "2pdm") {
    dn <- as.vector(outer((spec$i0 + 1L):I, gl,
                          function(i, g) paste0("d", i, "_g", g)))
    par$d <- matrix(unname(x[dn]), I - spec$i0, G)
    par$mu_dec <- unname(x[paste0("mu_dec_g", gl)])
    par$logit_pdec <- unname(x[paste0("logit_pdec_g", gl)])
  } else {
    if (spec$model == "hybrid")
      par$btilde <- unname(x[paste0("btilde_g", gl)])
    else par$logkappa <- unname(x[paste0("logkappa_g", gl)])
    par$rho <- unname(x[paste0("rho_g", gl)])
    par$tau1 <- unname(x[paste0("tau1_g", gl)])
    par$logomega <- unname(x[paste0("logomega_g", gl)])
  }
  if (G > 1L) {
    par$mu_nd[-1L] <- unname(x[paste0("mu_nd_g", gl[-1L])])
    par$logsigma[-1L] <- unname(x[paste0("logsigma_g", gl[-1L])])
  }
  par
}

## Natural-scale group parameters as a named vector (for delta method).
natural_par <- function(x, spec, template) {
  par <- unpack_par(x, spec, template)
  G <- spec$G; gl <- seq_len(G)
  out <- c()
  if (spec$model == "2pl") {
    ## only the shared tail below
  } else if (spec$model == "2pdm") {
    out <- c(setNames(plogis(par$logit_pdec), paste0("p_decline_g", gl)),
             setNames(par$mu_dec, paste0("mu_dec_g", gl)))
  } else {
    pi_nd <- vapply(gl, function(g)
      probs_from_logits(class_logits_g(spec$model, spec$I,
                                       tau1 = par$tau1[g],
                                       logomega = par$logomega[g]))[spec$nclass],
      numeric(1))
    out <- c(setNames(pi_nd, paste0("pi_nd_g", gl)),
             setNames(exp(par$logomega), paste0("omega_g", gl)),
             if (spec$model == "hybrid")
               setNames(par$btilde, paste0("btilde_g", gl))
             else setNames(exp(par$logkappa), paste0("kappa_g", gl)),
             setNames(par$rho, paste0("rho_g", gl)))
  }
  c(out, setNames(par$mu_nd, paste0("mu_nd_g", gl)),
    setNames(exp(par$logsigma), paste0("sigma_g", gl)))
}

num_jacobian <- function(f, x, eps = 1e-5) {
  f0 <- f(x)
  J <- matrix(0, length(f0), length(x),
              dimnames = list(names(f0), names(x)))
  for (j in seq_along(x)) {
    h <- eps * max(1, abs(x[j]))
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    J[, j] <- (f(xp) - f(xm)) / (2 * h)
  }
  J
}

#' Observed-information covariance matrix
#'
#' Computes the covariance matrix of the estimated parameters (on the
#' estimation scale) from the numerically differentiated observed
#' information. This requires on the order of \code{npar^2} marginal
#' log-likelihood evaluations and can take a few minutes for long tests.
#'
#' @param object A fitted \code{"mgmixirt"} object.
#' @param ... Unused.
#' @return A covariance matrix with parameter names.
#' @export
vcov.mgmixirt <- function(object, ...) {
  spec <- object$spec
  quad <- gh_quad(object$quadpts)
  Xg <- split.data.frame(object$resp, object$group)
  x0 <- pack_par(object$par, spec)
  fn <- function(x) estep(unpack_par(x, spec, object$par), spec, quad, Xg)$ll
  H <- stats::optimHess(x0, fn)
  info <- -(H + t(H)) / 2
  V <- tryCatch(solve(info), error = function(e) {
    warning("observed information is singular; using a pseudo-inverse. ",
            "Some parameters are empirically unidentified.")
    e <- eigen(info, symmetric = TRUE)
    pos <- e$values > max(e$values) * 1e-10
    e$vectors[, pos, drop = FALSE] %*%
      diag(1 / e$values[pos], sum(pos)) %*%
      t(e$vectors[, pos, drop = FALSE])
  })
  dimnames(V) <- list(names(x0), names(x0))
  V
}

#' Wald tests of group differences
#'
#' Estimates, standard errors, and Wald tests of equality across groups
#' for the natural-scale group parameters (class proportions, shape and
#' magnitude parameters, ability means and standard deviations),
#' following the group comparisons in List et al. (2017, Table 2).
#'
#' @param object A fitted \code{"mgmixirt"} object.
#' @param vcov Optional covariance matrix from \code{\link{vcov.mgmixirt}}
#'   (computed if missing; expensive for long tests).
#' @param ... Unused.
#' @return A data frame with one row per group parameter: estimates and
#'   standard errors per group, and (for G > 1) the Wald chi-square test
#'   of equality across groups with df = G - 1.
#' @export
wald_test <- function(object, ...) UseMethod("wald_test")

#' @rdname wald_test
#' @export
wald_test.mgmixirt <- function(object, vcov = NULL, ...) {
  spec <- object$spec; G <- spec$G
  if (is.null(vcov)) vcov <- vcov.mgmixirt(object)
  x0 <- pack_par(object$par, spec)
  natf <- function(x) natural_par(x, spec, object$par)
  est <- natf(x0)
  J <- num_jacobian(natf, x0)
  Vn <- J %*% vcov %*% t(J)
  se <- sqrt(pmax(diag(Vn), 0))
  fam <- sub("_g[0-9]+$", "", names(est))
  fams <- unique(fam)
  out <- data.frame(parameter = fams)
  for (g in seq_len(G)) {
    idx <- match(paste0(fams, "_g", g), names(est))
    out[[paste0("est_g", g)]] <- est[idx]
    out[[paste0("se_g", g)]] <- se[idx]
  }
  if (G > 1L) {
    Cmat <- cbind(1, diag(-1, G - 1L))
    w <- p <- rep(NA_real_, length(fams))
    for (k in seq_along(fams)) {
      idx <- match(paste0(fams[k], "_g", seq_len(G)), names(est))
      th <- est[idx]
      Vf <- Vn[idx, idx, drop = FALSE]
      CV <- Cmat %*% Vf %*% t(Cmat)
      ct <- as.vector(Cmat %*% th)
      wk <- tryCatch(drop(t(ct) %*% solve(CV, ct)),
                     error = function(e) NA_real_)
      w[k] <- wk
      p[k] <- if (is.na(wk)) NA_real_ else
        stats::pchisq(wk, df = G - 1L, lower.tail = FALSE)
    }
    out$wald_chisq <- w
    out$df <- G - 1L
    out$p_value <- p
  }
  rownames(out) <- NULL
  out
}
