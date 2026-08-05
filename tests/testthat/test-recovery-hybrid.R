test_that("HYBRID estimator reaches the MLE and recovers stable parameters", {
  skip_on_cran()
  set.seed(42)
  a_true <- exp(rnorm(15, 0, 0.2))
  b_true <- seq(1.2, -1.2, length.out = 15)
  truth <- list(pi_nd = c(0.50, 0.65), omega = c(2.5, 4),
                btilde = c(-2.5, -2.5), mu_nd = c(0, 0.5),
                rho = c(-0.05, -0.08), sigma = c(1, 0.95))
  sim <- sim_mgmixirt("hybrid", n = c(1500, 1500), I = 15, a = a_true,
                      b = b_true, pi_nd = truth$pi_nd,
                      omega = truth$omega, btilde = truth$btilde,
                      mu_nd = truth$mu_nd, rho = truth$rho,
                      sigma = truth$sigma, seed = 501)
  fit <- mgmixirt(sim$resp, sim$group, model = "hybrid", quadpts = 13,
                  starts = c(3, 10), maxit = 2000, tol = 1e-5,
                  seed = 11, verbose = FALSE)
  expect_true(fit$converged)
  expect_true(all(diff(fit$trace) > -1e-4))

  ## the fitted likelihood must dominate the generating parameters
  ## (class-size parameters trade off along a flat ridge in these models,
  ## so point recovery of pi_nd/omega/btilde is a design property studied
  ## in the simulation study, not an estimator property)
  spec <- MGmixIRT:::mgm_spec("hybrid", 15L, 2L)
  quad <- MGmixIRT:::gh_quad(13)
  par_t <- list(a = a_true, b = b_true, mu_nd = truth$mu_nd,
                logsigma = log(truth$sigma), rho = truth$rho,
                btilde = truth$btilde,
                tau1 = mapply(MGmixIRT:::tau1_from_target, truth$pi_nd,
                              truth$omega, MoreArgs = list(I = 15)),
                logomega = log(truth$omega))
  Xg <- split.data.frame(sim$resp, sim$group)
  ll_truth <- MGmixIRT:::estep(par_t, spec, quad, Xg)$ll
  expect_gte(fit$logLik, ll_truth)

  co <- coef(fit)
  expect_gt(cor(co$items$a, a_true), 0.90)
  expect_lt(mean(abs(co$items$b - b_true)), 0.25)
  expect_lt(abs(co$groups$mu_nd[2] - 0.5), 0.30)
  expect_gt(co$groups$sigma[2], 0.80)
  expect_lt(co$groups$sigma[2], 1.10)

  ## posterior classification must separate true decliners
  sc <- scores(fit)
  expect_gt(mean(sc$p_nodecline[sim$delta == 15]) -
            mean(sc$p_nodecline[sim$delta < 15]), 0.10)
})
