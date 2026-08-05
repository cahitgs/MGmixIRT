test_that("Gauss-Hermite quadrature matches standard normal moments", {
  q <- MGmixIRT:::gh_quad(15)
  expect_equal(sum(q$w), 1, tolerance = 1e-12)
  expect_equal(sum(q$w * q$z), 0, tolerance = 1e-10)
  expect_equal(sum(q$w * q$z^2), 1, tolerance = 1e-8)
  expect_equal(sum(q$w * q$z^4), 3, tolerance = 1e-6)
})

test_that("log_ipow_diff is exact and stable", {
  i <- 1:29
  omega <- 7.33
  expect_equal(MGmixIRT:::log_ipow_diff(i, omega),
               log(i^omega - (i - 1)^omega), tolerance = 1e-10)
  expect_true(all(is.finite(MGmixIRT:::log_ipow_diff(1:30, 40))))
})
