#' Plot methods for fitted multigroup mixture PD models
#'
#' Base-graphics displays mirroring the figures of List et al. (2017):
#' \describe{
#'   \item{\code{"classprob"}}{Cumulated decline-class probabilities
#'     across switching points per group (their Fig. 2). For the 2PDM,
#'     a bar chart of the decline-class proportions per group.}
#'   \item{\code{"items"}}{Item slope and intercept estimates by item
#'     position (their Fig. 1, single model).}
#'   \item{\code{"decline"}}{The model-specific representation of the
#'     decline magnitude: class-specific intercepts after the switching
#'     point (2PDM, their Fig. 3), the common aberrant response
#'     probability (HYBRID), or the ability decrement by switching point
#'     (MPDM).}
#' }
#'
#' @param x A fitted \code{"mgmixirt"} object.
#' @param type One of \code{"classprob"}, \code{"items"},
#'   \code{"decline"}.
#' @param ... Passed to the underlying plotting calls.
#' @return Invisibly, \code{x}.
#' @export
plot.mgmixirt <- function(x, type = c("classprob", "items", "decline"),
                          ...) {
  type <- match.arg(type)
  spec <- x$spec; G <- spec$G; I <- spec$I
  gl <- x$group_levels
  co <- coef(x)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  cols <- seq_len(G)

  if (type == "classprob") {
    if (x$model == "2pdm") {
      graphics::barplot(co$groups$p_decline, names.arg = gl,
                        ylab = "P(decline class)", xlab = "Group",
                        ylim = c(0, max(co$groups$p_decline) * 1.3), ...)
    } else {
      s <- seq_len(I - 1)
      cum <- sapply(seq_len(G), function(g) {
        pi_g <- probs_from_logits(group_class_logits(x$par, spec, g))
        cumsum(pi_g[s])
      })
      graphics::matplot(s, cum, type = "s", lty = 1, lwd = 2, col = cols,
                        xlab = expression(paste("Switching point ", delta)),
                        ylab = "Cumulated decline-class probability", ...)
      graphics::legend("topleft", legend = gl, col = cols, lty = 1,
                       lwd = 2, bty = "n")
    }
  } else if (type == "items") {
    graphics::par(mfrow = c(2, 1), mar = c(4, 4, 1.5, 1))
    graphics::plot(seq_len(I), co$items$a, type = "b", pch = 19,
                   xlab = "", ylab = "Slope (a)", ...)
    graphics::plot(seq_len(I), co$items$b, type = "b", pch = 19,
                   xlab = "Item position", ylab = "Intercept (b)", ...)
  } else {                                        # decline
    if (x$model == "2pdm") {
      pos <- seq_len(I)
      graphics::plot(pos, co$items$b, type = "b", pch = 19,
                     xlab = "Item position", ylab = "Intercept",
                     ylim = range(c(co$items$b,
                                    unlist(co$items[, -(1:2)])),
                                  na.rm = TRUE), ...)
      dec <- (spec$i0 + 1L):I
      for (g in seq_len(G))
        graphics::points(dec, co$items[dec, 2L + g], col = g + 1,
                         pch = 17, type = "b", lty = 2)
      graphics::abline(v = spec$i0 + 0.5, lty = 3)
      graphics::legend("bottomleft",
                       legend = c("no-decline", paste("decline,", gl)),
                       col = c(1, seq_len(G) + 1),
                       pch = c(19, rep(17, G)), bty = "n")
    } else if (x$model == "hybrid") {
      graphics::barplot(co$groups$p_aberrant, names.arg = gl,
                        ylab = "P(correct) after switching point",
                        xlab = "Group", ...)
    } else {
      s <- seq_len(I - 1)
      dec <- sapply(seq_len(G), function(g)
        co$groups$kappa[g] * (I - s))
      graphics::matplot(s, dec, type = "l", lty = 1, lwd = 2, col = cols,
                        xlab = expression(paste("Switching point ", delta)),
                        ylab = "Ability decrement", ...)
      graphics::legend("topright", legend = gl, col = cols, lty = 1,
                       lwd = 2, bty = "n")
    }
  }
  invisible(x)
}
