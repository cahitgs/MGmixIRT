## MML-EM engine.
##
## The marginal likelihood is approximated with Gauss-Hermite quadrature
## tied to the class-specific ability distributions: for group g, class c,
## node q, theta = mu_cg + sigma_g * z_q. The EM algorithm treats the pair
## (class, node) as the latent variable; the M-step maximises the expected
## complete-data log-likelihood by cyclic ascent over parameter blocks
## (items; aberrant-regime parameters; ability distribution; class
## probabilities), each of which increases the observed-data likelihood
## (generalised EM).

group_class_logits <- function(par, spec, g) {
  if (spec$model == "2pl") {
    0
  } else if (spec$model == "2pdm") {
    class_logits_g("2pdm", spec$I, logit_pdec = par$logit_pdec[g])
  } else {
    class_logits_g(spec$model, spec$I, tau1 = par$tau1[g],
                   logomega = par$logomega[g])
  }
}

estep <- function(par, spec, quad, Xg, keep_post = FALSE) {
  G <- spec$G
  stats <- vector("list", G)
  total_ll <- 0
  for (g in seq_len(G)) {
    cl <- cells_g(par, spec, g, quad)
    eta <- eta_g(par, spec, g, cl)
    lp <- log_p_logis(eta)
    lq <- log_q_logis(eta)
    X <- Xg[[g]]
    LL <- X %*% lp + (1 - X) %*% lq
    pi_g <- probs_from_logits(group_class_logits(par, spec, g))
    prior <- rep(log(pi_g), each = cl$Q) + log(rep(quad$w, spec$nclass))
    J <- sweep(LL, 2L, prior, "+")
    llp <- row_logsumexp(J)
    R <- exp(J - llp)
    stats[[g]] <- list(cN = colSums(R), cN1 = crossprod(X, R), cl = cl,
                       post = if (keep_post) R)
    total_ll <- total_ll + sum(llp)
  }
  list(ll = total_ll, stats = stats)
}

## Expected complete-data log-likelihood contribution of the response part
## for one group, given a linear-predictor matrix.
q_resp <- function(eta, cN1, cN) {
  n0 <- matrix(cN, nrow(eta), ncol(eta), byrow = TRUE) - cN1
  sum(cN1 * log_p_logis(eta) + n0 * log_q_logis(eta))
}

## Vectorised damped-Newton update of (a_i, b_i) for all items whose
## linear predictor depends on them alone (all items for HYBRID/MPDM;
## pre-switch items for the 2PDM). Per-item 2x2 Newton systems are solved
## simultaneously; items whose expected complete-data log-likelihood would
## decrease get their step halved individually.
newton_items <- function(par, spec, stats, sweeps = 4L) {
  I <- spec$I; G <- spec$G
  ok <- if (spec$model == "2pdm")
    seq_len(I) <= spec$i0 else rep(TRUE, I)
  covs <- masks <- vector("list", G)
  for (g in seq_len(G)) {
    cl <- stats[[g]]$cl
    C <- length(cl$delta)
    A <- outer(seq_len(I), cl$delta, ">")
    TH <- matrix(cl$theta, I, C, byrow = TRUE)
    if (spec$model == "mpdm") {
      TH <- TH - exp(par$logkappa[g]) * A *
        matrix(spec$I - cl$delta, I, C, byrow = TRUE)
      M <- matrix(1, I, C)
    } else if (spec$model == "hybrid") {
      M <- 1 - A
    } else {
      M <- matrix(1, I, C)
    }
    covs[[g]] <- TH; masks[[g]] <- M
  }
  qitems <- function(a, b) {
    p2 <- par; p2$a <- a; p2$b <- b
    Q <- numeric(I)
    for (g in seq_len(G)) {
      eta <- eta_g(p2, spec, g, stats[[g]]$cl)
      N1 <- stats[[g]]$cN1
      N0 <- matrix(stats[[g]]$cN, I, ncol(eta), byrow = TRUE) - N1
      Q <- Q + rowSums(N1 * log_p_logis(eta) + N0 * log_q_logis(eta))
    }
    Q
  }
  a <- par$a; b <- par$b
  Qcur <- qitems(a, b)
  for (sw in seq_len(sweeps)) {
    ga <- gb <- haa <- hab <- hbb <- numeric(I)
    p2 <- par; p2$a <- a; p2$b <- b
    for (g in seq_len(G)) {
      eta <- eta_g(p2, spec, g, stats[[g]]$cl)
      P <- plogis(eta)
      N1 <- stats[[g]]$cN1
      Ncell <- matrix(stats[[g]]$cN, I, ncol(eta), byrow = TRUE)
      R <- (N1 - Ncell * P) * masks[[g]]
      W <- Ncell * P * (1 - P) * masks[[g]]
      TH <- covs[[g]]
      ga <- ga + rowSums(R * TH);        gb <- gb + rowSums(R)
      haa <- haa + rowSums(W * TH * TH); hab <- hab + rowSums(W * TH)
      hbb <- hbb + rowSums(W)
    }
    det <- haa * hbb - hab^2
    upd <- ok & det > 1e-10
    da <- db <- numeric(I)
    da[upd] <- (hbb[upd] * ga[upd] - hab[upd] * gb[upd]) / det[upd]
    db[upd] <- (haa[upd] * gb[upd] - hab[upd] * ga[upd]) / det[upd]
    fac <- rep(1, I)
    for (h in seq_len(6L)) {
      a2 <- clamp(a + fac * da, 0.01, 15)
      b2 <- clamp(b + fac * db, -30, 30)
      Qnew <- qitems(a2, b2)
      worse <- upd & (Qnew < Qcur - 1e-10) & fac > 0
      if (!any(worse)) break
      fac[worse] <- if (h == 6L) 0 else fac[worse] / 2
    }
    a <- clamp(a + fac * da, 0.01, 15)
    b <- clamp(b + fac * db, -30, 30)
    Qcur <- qitems(a, b)
  }
  par$a <- a; par$b <- b
  par
}

mstep <- function(par, spec, quad, stats) {
  G <- spec$G; I <- spec$I
  model <- spec$model

  ## --- item block ---------------------------------------------------
  par <- newton_items(par, spec, stats)
  if (model == "2pdm") {
    ## post-switch items carry class-specific decrements d_ig: small
    ## nlminb problems per item
    for (i in (spec$i0 + 1L):I) {
      xi0 <- c(par$a[i], par$b[i], par$d[i - spec$i0, ])
      negQ <- function(xi) {
        s <- 0
        for (g in seq_len(G)) {
          eta <- eta_item_g(i, xi[1L], xi[2L], d_ig = xi[2L + g],
                            par = par, spec = spec, g = g,
                            cl = stats[[g]]$cl)
          n1 <- stats[[g]]$cN1[i, ]
          s <- s + sum(n1 * log_p_logis(eta) +
                       (stats[[g]]$cN - n1) * log_q_logis(eta))
        }
        -s
      }
      opt <- nlminb(xi0, negQ,
                    lower = c(0.01, -30, rep(-30, G)),
                    upper = c(15, 30, rep(10, G)))
      par$a[i] <- opt$par[1L]
      par$b[i] <- opt$par[2L]
      par$d[i - spec$i0, ] <- opt$par[-(1:2)]
    }
  }

  ## --- HYBRID common aberrant intercept (closed form) ---------------
  if (model == "hybrid") {
    for (g in seq_len(G)) {
      cl <- stats[[g]]$cl
      ab <- outer(seq_len(I), cl$delta, ">")
      if (any(ab)) {
        S1 <- sum(stats[[g]]$cN1[ab])
        S <- sum(matrix(stats[[g]]$cN, I, length(cl$delta),
                        byrow = TRUE)[ab])
        par$btilde[g] <- qlogis(clamp(S1 / S, 1e-8, 1 - 1e-8))
      }
    }
  }

  ## --- ability distribution block (per group) -----------------------
  for (g in seq_len(G)) {
    fields <- switch(model,
      "2pl"    = if (g > 1L) c("mu_nd", "logsigma"),
      "2pdm"   = c("mu_dec", if (g > 1L) c("mu_nd", "logsigma")),
      "hybrid" = c("rho",    if (g > 1L) c("mu_nd", "logsigma")),
      "mpdm"   = c("rho", "logkappa", if (g > 1L) c("mu_nd", "logsigma")))
    if (length(fields) == 0L) next
    psi0 <- vapply(fields, function(f) par[[f]][g], numeric(1))
    negQ_dist <- function(psi) {
      p2 <- par
      for (k in seq_along(fields)) p2[[fields[k]]][g] <- psi[k]
      cl <- cells_g(p2, spec, g, quad)
      -q_resp(eta_g(p2, spec, g, cl), stats[[g]]$cN1, stats[[g]]$cN)
    }
    lo <- ifelse(fields == "logsigma", log(0.05),
          ifelse(fields == "logkappa", log(1e-4), -10))
    hi <- ifelse(fields == "logsigma", log(10),
          ifelse(fields == "logkappa", log(10), 10))
    opt <- nlminb(psi0, negQ_dist, lower = lo, upper = hi)
    for (k in seq_along(fields)) par[[fields[k]]][g] <- opt$par[k]
  }

  ## --- class probability block --------------------------------------
  if (model == "2pl") return(par)
  for (g in seq_len(G)) {
    Q <- stats[[g]]$cl$Q
    m <- colSums(matrix(stats[[g]]$cN, nrow = Q))
    if (model == "2pdm") {
      p_dec <- clamp(m[1L] / sum(m), 1e-8, 1 - 1e-8)
      par$logit_pdec[g] <- qlogis(p_dec)
    } else {
      negQ_cls <- function(psi) {
        tau <- class_logits_g(model, I, tau1 = psi[1L], logomega = psi[2L])
        -sum(m * log(probs_from_logits(tau)))
      }
      opt <- nlminb(c(par$tau1[g], par$logomega[g]), negQ_cls,
                    lower = c(-30, log(0.05)), upper = c(30, log(50)))
      par$tau1[g] <- opt$par[1L]
      par$logomega[g] <- opt$par[2L]
    }
  }
  par
}

em_run <- function(par, spec, quad, Xg, maxit, tol, verbose = FALSE) {
  ll_old <- -Inf
  trace <- numeric(0)
  converged <- FALSE
  nonmono <- 0L
  for (it in seq_len(maxit)) {
    es <- estep(par, spec, quad, Xg)
    trace[it] <- es$ll
    if (it > 1L) {
      if (es$ll < ll_old - 1e-6) nonmono <- nonmono + 1L
      if (abs(es$ll - ll_old) < tol) {
        converged <- TRUE
        break
      }
    }
    ll_old <- es$ll
    par <- mstep(par, spec, quad, es$stats)
    if (verbose && it %% 10L == 0L)
      message(sprintf("iter %4d  logLik %.4f", it, es$ll))
  }
  list(par = par, ll = trace[length(trace)], trace = trace,
       iterations = length(trace), converged = converged,
       nonmonotone = nonmono)
}

start_par <- function(spec, Xg, quad) {
  I <- spec$I; G <- spec$G
  pbar <- colMeans(do.call(rbind, Xg))
  totals <- lapply(Xg, rowSums)
  m1 <- mean(totals[[1L]]); s1 <- stats::sd(totals[[1L]]) + 1e-9
  mu_nd <- vapply(totals, function(t) (mean(t) - m1) / s1, numeric(1))
  mu_nd[1L] <- 0
  par <- list(a = rep(1, I),
              b = qlogis(clamp(pbar, 0.02, 0.98)),
              mu_nd = mu_nd,
              logsigma = rep(0, G))
  if (spec$model == "2pdm") {
    par$d <- matrix(log(0.7), I - spec$i0, G)
    par$mu_dec <- mu_nd - 0.1
    par$logit_pdec <- rep(qlogis(0.15), G)
  } else if (spec$model != "2pl") {
    par$rho <- rep(-0.02, G)
    par$tau1 <- rep(tau1_from_target(0.8, 5, I), G)
    par$logomega <- rep(log(5), G)
    if (spec$model == "hybrid") par$btilde <- rep(qlogis(0.10), G)
    if (spec$model == "mpdm") par$logkappa <- rep(log(0.5), G)
  }
  par
}

jitter_par <- function(par, spec) {
  G <- spec$G
  par$a <- par$a * exp(rnorm(length(par$a), 0, 0.15))
  par$b <- par$b + rnorm(length(par$b), 0, 0.3)
  if (G > 1L) par$mu_nd[-1L] <- par$mu_nd[-1L] + rnorm(G - 1L, 0, 0.3)
  if (spec$model == "2pdm") {
    par$d <- par$d + rnorm(length(par$d), 0, 0.3)
    par$mu_dec <- par$mu_dec + rnorm(G, 0, 0.3)
    par$logit_pdec <- par$logit_pdec + rnorm(G, 0, 0.5)
  } else if (spec$model != "2pl") {
    par$rho <- par$rho + rnorm(G, 0, 0.02)
    par$tau1 <- par$tau1 + rnorm(G, 0, 0.5)
    par$logomega <- par$logomega + rnorm(G, 0, 0.3)
    if (spec$model == "hybrid") par$btilde <- par$btilde + rnorm(G, 0, 0.4)
    if (spec$model == "mpdm") par$logkappa <- par$logkappa + rnorm(G, 0, 0.3)
  }
  par
}
