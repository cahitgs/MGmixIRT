test_that("scores() returns coherent person-level output", {
  skip_on_cran()
  sim <- sim_mgmixirt("2pdm", n = c(600, 600), I = 12, switch_point = 8,
                      pdec = c(0.35, 0.2), dshift = 2, seed = 9)
  fit <- mgmixirt(sim$resp, sim$group, model = "2pdm", switch_point = 8,
                  quadpts = 9, starts = c(2, 8), maxit = 300, tol = 1e-4,
                  seed = 3, verbose = FALSE)
  sc <- scores(fit)
  expect_equal(nrow(sc), 1200L)
  expect_true(all(is.finite(sc$eap)))
  expect_true(all(sc$eap_sd > 0))
  expect_true(all(sc$delta_map %in% c(8L, 12L)))
  ## true decliners should receive lower posterior no-decline probability
  expect_gt(mean(sc$p_nodecline[sim$delta == 12]) -
            mean(sc$p_nodecline[sim$delta == 8]), 0.15)
  ## EAP should track the generating abilities
  expect_gt(cor(sc$eap, sim$theta), 0.7)
})
