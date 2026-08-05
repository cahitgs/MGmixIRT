#' Fit a multigroup mixture performance decline model
#'
#' Estimates the multigroup 2PDM (Bolt et al., 2002), HYBRID
#' (Yamamoto, 1995), or MPDM (Jin & Wang, 2014) of List et al. (2017) by
#' marginal maximum likelihood (EM algorithm with Gauss-Hermite
#' quadrature). The first level of \code{group} serves as the reference
#' group, whose no-decline ability distribution is standardised to
#' N(0, 1) for identification.
#'
#' @param resp Matrix or data frame of dichotomous (0/1) item responses in
#'   administration order. Missing values are not supported in this
#'   version; recode omitted/not-reached responses before fitting (List et
#'   al., 2017, coded them as incorrect).
#' @param group Factor (or vector coercible to one) of manifest group
#'   membership, one entry per row of \code{resp}. \code{NULL} fits a
#'   single-group model.
#' @param model One of \code{"2pdm"}, \code{"hybrid"}, \code{"mpdm"}.
#' @param switch_point 2PDM only: the pre-specified switching point
#'   \code{i0} (responses up to item \code{i0} are unaffected by decline).
#'   A vector of candidate values triggers a grid search; the model with
#'   the lowest BIC is returned (with the search table attached).
#' @param quadpts Number of Gauss-Hermite quadrature points.
#' @param starts Length-2 integer vector \code{c(n_starts, short_iter)}:
#'   number of random starts and the number of EM iterations for each
#'   short run before the best solution is iterated to convergence.
#' @param maxit,tol EM convergence controls: iteration cap and the
#'   absolute log-likelihood change below which the algorithm stops.
#' @param seed Optional RNG seed for the random starts.
#' @param verbose Print progress messages.
#'
#' @return An object of class \code{"mgmixirt"}.
#' @export
mgmixirt <- function(resp, group = NULL,
                     model = c("2pdm", "hybrid", "mpdm"),
                     switch_point = NULL, quadpts = 15,
                     starts = c(10, 20), maxit = 2000, tol = 1e-5,
                     seed = NULL, verbose = TRUE) {
  cl_call <- match.call()
  model <- match.arg(model)
  X <- as.matrix(resp)
  storage.mode(X) <- "double"
  if (anyNA(X))
    stop("missing responses are not supported in this version; ",
         "recode them before fitting")
  if (!all(X %in% c(0, 1)))
    stop("resp must contain only 0/1 values")
  n <- nrow(X); I <- ncol(X)
  if (is.null(group)) group <- rep(1L, n)
  group <- droplevels(as.factor(group))
  if (length(group) != n) stop("length(group) must equal nrow(resp)")
  G <- nlevels(group)

  if (model == "2pdm" && length(switch_point) > 1L) {
    cand <- as.integer(switch_point)
    fits <- lapply(cand, function(i0)
      mgmixirt(X, group, model = "2pdm", switch_point = i0,
               quadpts = quadpts, starts = starts, maxit = maxit,
               tol = tol, seed = seed, verbose = FALSE))
    bic <- vapply(fits, BIC, numeric(1))
    best <- which.min(bic)
    if (verbose)
      message(sprintf("switch-point search: best i0 = %d (BIC = %.2f)",
                      cand[best], bic[best]))
    out <- fits[[best]]
    out$switch_search <- data.frame(i0 = cand,
                                    logLik = vapply(fits, function(f)
                                      as.numeric(logLik(f)), numeric(1)),
                                    BIC = bic)
    out$call <- cl_call
    return(out)
  }

  spec <- mgm_spec(model, I, G, i0 = switch_point)
  quad <- gh_quad(quadpts)
  Xg <- split.data.frame(X, group)

  if (!is.null(seed)) set.seed(seed)
  base <- start_par(spec, Xg, quad)
  n_starts <- max(1L, as.integer(starts[1L]))
  short <- max(1L, as.integer(starts[2L]))
  cands <- c(list(base),
             replicate(n_starts - 1L, jitter_par(base, spec),
                       simplify = FALSE))
  short_runs <- lapply(cands, function(p)
    em_run(p, spec, quad, Xg, maxit = short, tol = tol))
  best <- which.max(vapply(short_runs, `[[`, numeric(1), "ll"))
  if (verbose && n_starts > 1L)
    message(sprintf("multi-start: best of %d short runs, logLik %.4f",
                    n_starts, short_runs[[best]]$ll))
  run <- em_run(short_runs[[best]]$par, spec, quad, Xg,
                maxit = maxit, tol = tol, verbose = verbose)
  if (!run$converged)
    warning("EM did not converge within maxit iterations")
  if (run$nonmonotone > 0L)
    warning(sprintf("log-likelihood decreased in %d iteration(s); ",
                    run$nonmonotone),
            "consider increasing quadpts")

  fin <- estep(run$par, spec, quad, Xg, keep_post = TRUE)
  post <- vector("list", G)
  for (g in seq_len(G)) {
    R <- fin$stats[[g]]$post
    Qn <- fin$stats[[g]]$cl$Q
    P <- sapply(seq_len(spec$nclass), function(c)
      rowSums(R[, ((c - 1L) * Qn + 1L):(c * Qn), drop = FALSE]))
    post[[g]] <- matrix(P, ncol = spec$nclass)
  }
  Pall <- do.call(rbind, post)
  ent <- 1 - sum(-Pall * log(pmax(Pall, 1e-300))) /
    (nrow(Pall) * log(spec$nclass))

  structure(list(call = cl_call, model = model, spec = spec,
                 par = run$par, logLik = run$ll, npar = mgm_npar(spec),
                 nobs = n, group_levels = levels(group),
                 group_sizes = as.integer(table(group)),
                 posterior = post, entropy = ent,
                 trace = run$trace, iterations = run$iterations,
                 converged = run$converged, quadpts = quadpts,
                 resp = X, group = group),
            class = "mgmixirt")
}
