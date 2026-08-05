test_that("simulated data have the expected structure", {
  sim <- sim_mgmixirt("hybrid", n = c(300, 200), I = 12, seed = 1)
  expect_equal(dim(sim$resp), c(500, 12))
  expect_true(all(sim$resp %in% 0:1))
  expect_equal(levels(sim$group), c("1", "2"))
  expect_true(all(sim$delta %in% 1:12))
})

test_that("decline lowers end-of-test performance", {
  sim <- sim_mgmixirt("hybrid", n = 4000, I = 20,
                      pi_nd = 0.5, omega = 3, btilde = -4,
                      mu_nd = 0, rho = 0, seed = 2)
  dec <- sim$delta < 20
  last5 <- rowMeans(sim$resp[, 16:20])
  expect_gt(mean(last5[!dec]) - mean(last5[dec]), 0.1)
})

test_that("2pdm simulation respects the switching point", {
  sim <- sim_mgmixirt("2pdm", n = c(500, 500), I = 10, switch_point = 6,
                      pdec = c(0.4, 0.2), dshift = 2, seed = 3)
  expect_true(all(sim$delta %in% c(6L, 10L)))
})
