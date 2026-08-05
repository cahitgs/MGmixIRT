#' Person scores and class assignments
#'
#' Computes expected a posteriori (EAP) ability estimates, their posterior
#' standard deviations, posterior class probabilities, and modal class
#' assignments from a fitted multigroup mixture PD model.
#'
#' Following List et al. (2017), persons can be dichotomised into a no-PD
#' and a PD group by thresholding the posterior probability of the
#' no-decline class at 0.5 (column \code{pd_class}).
#'
#' @param object A fitted \code{"mgmixirt"} object.
#' @param ... Unused.
#'
#' @return A data frame with one row per person (in the original row
#'   order of the response matrix): \code{group}, \code{eap},
#'   \code{eap_sd}, \code{p_nodecline} (posterior probability of the
#'   no-decline class), \code{delta_map} (modal switching point), and
#'   \code{pd_class} (\code{"pd"}/\code{"no-pd"} by the 0.5 rule).
#' @export
scores <- function(object, ...) UseMethod("scores")

#' @export
scores.mgmixirt <- function(object, ...) {
  spec <- object$spec
  quad <- gh_quad(object$quadpts)
  Xg <- split.data.frame(object$resp, object$group)
  es <- estep(object$par, spec, quad, Xg, keep_post = TRUE)
  idx <- split(seq_len(object$nobs), object$group)
  out <- data.frame(group = object$group,
                    eap = NA_real_, eap_sd = NA_real_,
                    p_nodecline = NA_real_, delta_map = NA_integer_)
  for (g in seq_len(spec$G)) {
    R <- es$stats[[g]]$post
    cl <- es$stats[[g]]$cl
    m1 <- as.vector(R %*% cl$theta)
    m2 <- as.vector(R %*% cl$theta^2)
    Qn <- cl$Q
    Pcls <- sapply(seq_len(spec$nclass), function(cc)
      rowSums(R[, ((cc - 1L) * Qn + 1L):(cc * Qn), drop = FALSE]))
    Pcls <- matrix(Pcls, ncol = spec$nclass)
    rows <- idx[[g]]
    out$eap[rows] <- m1
    out$eap_sd[rows] <- sqrt(pmax(m2 - m1^2, 0))
    out$p_nodecline[rows] <- Pcls[, spec$nclass]
    out$delta_map[rows] <- spec$deltas[max.col(Pcls, ties.method = "last")]
  }
  out$pd_class <- ifelse(out$p_nodecline > 0.5, "no-pd", "pd")
  out
}
