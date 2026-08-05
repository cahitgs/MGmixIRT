## Internal numerical utilities.

## Gauss-Hermite quadrature for integration against the standard normal
## density: returns nodes z and weights w with sum(w) = 1, so that
## E[f(Z)] approx sum(w * f(z)) for Z ~ N(0, 1).
## Nodes/weights obtained by the Golub-Welsch algorithm for the
## physicists' Hermite weight exp(-x^2), then rescaled.
gh_quad <- function(Q) {
  Q <- as.integer(Q)
  if (Q < 1L) stop("quadpts must be >= 1")
  if (Q == 1L) return(list(z = 0, w = 1))
  k <- seq_len(Q - 1L)
  J <- matrix(0, Q, Q)
  off <- sqrt(k / 2)
  J[cbind(k, k + 1L)] <- off
  J[cbind(k + 1L, k)] <- off
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  x <- e$values[ord]
  w <- e$vectors[1L, ord]^2
  list(z = sqrt(2) * x, w = w / sum(w))
}

## Row-wise log-sum-exp of a numeric matrix.
row_logsumexp <- function(M) {
  m <- do.call(pmax, c(as.data.frame(M), list(na.rm = TRUE)))
  m + log(rowSums(exp(M - m)))
}

## log(i^omega - (i-1)^omega) for integer i >= 1, numerically stable for
## large omega. Vectorised over i.
log_ipow_diff <- function(i, omega) {
  out <- numeric(length(i))
  pos <- i > 1
  if (any(pos)) {
    ip <- i[pos]
    out[pos] <- omega * log(ip) +
      log1p(-exp(omega * (log(ip - 1) - log(ip))))
  }
  out
}

## Stable log(p) and log(1-p) for a logistic argument.
log_p_logis  <- function(eta) plogis(eta, log.p = TRUE)
log_q_logis  <- function(eta) plogis(-eta, log.p = TRUE)

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
