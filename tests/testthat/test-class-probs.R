test_that("omega-parameterised logits reproduce Cao-Stokes probabilities", {
  I <- 10; omega <- 2.5; pi_nd <- 0.8
  tau1 <- MGmixIRT:::tau1_from_target(pi_nd, omega, I)
  tau <- MGmixIRT:::class_logits_g("hybrid", I, tau1 = tau1,
                                   logomega = log(omega))
  pr <- MGmixIRT:::probs_from_logits(tau)
  i <- seq_len(I - 1)
  direct <- (i^omega - (i - 1)^omega) / (I - 1)^omega * (1 - pi_nd)
  expect_equal(sum(pr), 1, tolerance = 1e-12)
  expect_equal(pr[I], pi_nd, tolerance = 1e-10)
  expect_equal(pr[1:(I - 1)], direct, tolerance = 1e-10)
})

test_that("2pdm class logits give the decline probability", {
  tau <- MGmixIRT:::class_logits_g("2pdm", 20, logit_pdec = qlogis(0.3))
  pr <- MGmixIRT:::probs_from_logits(tau)
  expect_equal(pr[1], 0.3, tolerance = 1e-12)
})
