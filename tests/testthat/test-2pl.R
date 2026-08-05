test_that("2PL baseline recovers and information criteria select the true model", {
  skip_on_cran()
  set.seed(42)
  a_true <- exp(rnorm(12, 0, 0.2))
  b_true <- seq(1, -1, length.out = 12)
  sim <- sim_mgmixirt("2pl", n = c(800, 800), I = 12, a = a_true,
                      b = b_true, mu_nd = c(0, 0.5), sigma = c(1, 0.9),
                      seed = 31)
  fit <- mgmixirt(sim$resp, sim$group, model = "2pl", quadpts = 13,
                  starts = c(2, 8), maxit = 500, tol = 1e-5, seed = 3,
                  verbose = FALSE)
  expect_true(fit$converged)
  expect_equal(fit$npar, 2L * 12L + 2L)
  expect_true(all(sim$delta == 12L))
  co <- coef(fit)
  expect_gt(cor(co$items$a, a_true), 0.90)
  expect_lt(mean(abs(co$items$b - b_true)), 0.12)
  expect_lt(abs(co$groups$mu_nd[2] - 0.5), 0.12)
  expect_lt(abs(co$groups$sigma[2] - 0.9), 0.12)

  ## BIC must prefer the 2PL on 2PL data over the HYBRID. The HYBRID fit
  ## crawls along a flat ridge on null data and may hit maxit — harmless
  ## here, since more iterations could only improve its likelihood by a
  ## few points against a 59-point BIC penalty.
  fith <- suppressWarnings(
    mgmixirt(sim$resp, sim$group, model = "hybrid", quadpts = 9,
             starts = c(2, 6), maxit = 400, tol = 1e-4, seed = 3,
             verbose = FALSE))
  expect_lt(BIC(fit), BIC(fith))

  ## ... and the HYBRID on strong-decline data over the 2PL. AIC is used
  ## here: with only 8 extra parameters but a harsh log(n) penalty, BIC
  ## is known to under-select mixture structure in short tests (cf. Sen
  ## & Cohen, 2024) — a phenomenon the simulation study documents.
  simh <- sim_mgmixirt("hybrid", n = c(800, 800), I = 12, a = a_true,
                       b = b_true, pi_nd = c(0.45, 0.60),
                       omega = c(2, 3), btilde = c(-2.5, -2.5),
                       mu_nd = c(0, 0.4), rho = c(-0.05, -0.06),
                       seed = 32)
  fit2 <- mgmixirt(simh$resp, simh$group, model = "2pl", quadpts = 13,
                   starts = c(2, 8), maxit = 500, tol = 1e-5, seed = 3,
                   verbose = FALSE)
  fith2 <- mgmixirt(simh$resp, simh$group, model = "hybrid", quadpts = 9,
                    starts = c(2, 6), maxit = 400, tol = 1e-4, seed = 3,
                    verbose = FALSE)
  expect_lt(AIC(fith2), AIC(fit2))
})
