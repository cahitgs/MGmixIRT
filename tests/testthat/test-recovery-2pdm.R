test_that("2PDM recovers generating parameters (two groups)", {
  skip_on_cran()
  set.seed(42)
  a_true <- exp(rnorm(15, 0, 0.2))
  b_true <- seq(1.2, -1.2, length.out = 15)
  ## group 1 is the fitted reference group, so generate with the same
  ## convention (mu_nd[1] = 0, sigma[1] = 1) to keep scales aligned
  sim <- sim_mgmixirt("2pdm", n = c(1500, 1500), I = 15, a = a_true,
                      b = b_true, switch_point = 10,
                      pdec = c(0.30, 0.15), dshift = 1.5,
                      mu_nd = c(0, 0.8), mu_dec = c(-0.1, 0.7),
                      seed = 7)
  fit <- mgmixirt(sim$resp, sim$group, model = "2pdm", switch_point = 10,
                  quadpts = 11, starts = c(4, 15), maxit = 800,
                  tol = 1e-4, seed = 11, verbose = FALSE)
  co <- coef(fit)
  expect_true(fit$converged)
  expect_gt(cor(co$items$a, a_true), 0.85)
  expect_lt(mean(abs(co$items$b - b_true)), 0.25)
  expect_lt(max(abs(co$groups$p_decline - c(0.30, 0.15))), 0.08)
  expect_lt(abs((co$groups$mu_nd[2] - co$groups$mu_nd[1]) - 0.8), 0.2)
  ## log-likelihood must be non-decreasing over EM iterations
  expect_true(all(diff(fit$trace) > -1e-4))
})
