#' @export
logLik.mgmixirt <- function(object, ...) {
  structure(object$logLik, df = object$npar, nobs = object$nobs,
            class = "logLik")
}

#' @export
coef.mgmixirt <- function(object, ...) {
  par <- object$par; spec <- object$spec
  items <- data.frame(a = par$a, b = par$b,
                      row.names = colnames(object$resp))
  if (spec$model == "2pdm") {
    bt <- par$b[(spec$i0 + 1L):spec$I] -
      exp(par$d)                              # class-specific intercepts
    colnames(bt) <- paste0("b_dec.", object$group_levels)
    items <- cbind(items,
                   rbind(matrix(NA_real_, spec$i0, spec$G,
                                dimnames = list(NULL, colnames(bt))), bt))
  }
  gl <- object$group_levels
  groups <- switch(spec$model,
    "2pl" = data.frame(
      row.names = gl,
      mu_nd = par$mu_nd, sigma = exp(par$logsigma)),
    "2pdm" = data.frame(
      row.names = gl,
      p_decline = plogis(par$logit_pdec),
      mu_nd = par$mu_nd, mu_dec = par$mu_dec,
      sigma = exp(par$logsigma)),
    "hybrid" = data.frame(
      row.names = gl,
      pi_nd = vapply(seq_len(spec$G), function(g)
        probs_from_logits(group_class_logits(par, spec, g))[spec$nclass],
        numeric(1)),
      omega = exp(par$logomega), btilde = par$btilde,
      p_aberrant = plogis(par$btilde),
      mu_nd = par$mu_nd, rho = par$rho, sigma = exp(par$logsigma)),
    "mpdm" = data.frame(
      row.names = gl,
      pi_nd = vapply(seq_len(spec$G), function(g)
        probs_from_logits(group_class_logits(par, spec, g))[spec$nclass],
        numeric(1)),
      omega = exp(par$logomega), kappa = exp(par$logkappa),
      mu_nd = par$mu_nd, rho = par$rho, sigma = exp(par$logsigma)))
  list(items = items, groups = groups)
}

#' @export
print.mgmixirt <- function(x, ...) {
  hdr <- c("2pdm" = "two-class mixture PD model (Bolt et al., 2002)",
           "hybrid" = "HYBRID model (Yamamoto, 1995)",
           "mpdm" = "multiclass mixture PD model (Jin & Wang, 2014)",
           "2pl" = "2PL model (no mixture)")
  cat("Multigroup", hdr[[x$model]], "\n")
  if (x$model == "2pdm")
    cat("Switching point i0 =", x$spec$i0, "\n")
  cat(sprintf("%d persons, %d items, %d group(s): %s\n", x$nobs,
              x$spec$I, x$spec$G,
              paste(sprintf("%s (n=%d)", x$group_levels, x$group_sizes),
                    collapse = ", ")))
  ll <- logLik(x)
  cat(sprintf("logLik %.3f on %d parameters | AIC %.2f | BIC %.2f\n",
              as.numeric(ll), attr(ll, "df"), AIC(x), BIC(x)))
  cat(sprintf("Entropy %.3f | EM iterations %d (%s)\n", x$entropy,
              x$iterations,
              if (x$converged) "converged" else "NOT converged"))
  gr <- coef(x)$groups
  cat("\nGroup parameters:\n")
  print(round(gr, 3))
  invisible(x)
}

#' @export
summary.mgmixirt <- function(object, ...) {
  out <- list(fit = object, coefs = coef(object))
  class(out) <- "summary.mgmixirt"
  out
}

#' @export
print.summary.mgmixirt <- function(x, ...) {
  print(x$fit)
  cat("\nItem parameters:\n")
  print(round(x$coefs$items, 3))
  invisible(x)
}
