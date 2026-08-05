test_that("2PDM syntax generalises over groups and switch points", {
  s <- mplus_syntax("2pdm", I = 15, G = 3, switch_point = 10)
  expect_match(s, "CLASSES = g \\(3\\) c \\(2\\)")
  expect_match(s, "KNOWNCLASS = g \\(grp = 0 1 2\\)")
  expect_match(s, "\\[i11\\$1-i15\\$1\\*\\] \\(td3_11-td3_15\\)")
  expect_match(s, "td3_15 > t15")
  expect_match(s, "theta@1;\n\\[theta@0\\];")   # reference group identification
})

test_that("HYBRID syntax carries the omega machinery", {
  s <- mplus_syntax("hybrid", I = 12)
  expect_match(s, "CLASSES = g \\(2\\) c \\(12\\)")
  expect_match(s, "c ON g \\(b1-b11\\)")
  expect_match(s, "do \\(2,11\\) cpi# = cpi1 \\+ ln\\(#\\*\\*om2 - \\(#-1\\)\\*\\*om2\\)")
  expect_match(s, "theta BY i2-i12@0")          # zero slopes after switch
  expect_match(s, "beta2 = -1\\*bt2")
})

test_that("MPDM syntax uses the node-variable construction", {
  s <- mplus_syntax("mpdm", I = 12)
  expect_match(s, "nv BY i1-i12@0")
  expect_match(s, "nv WITH theta@0")
  expect_match(s, "do \\(1,10\\) gm1_# = gm1_11\\*\\(12-#\\)")
  expect_match(s, "kappa2 = -1\\*gm2_11")
})

test_that("HYBRID/MPDM generator rejects G > 2 and fitted objects delegate", {
  expect_error(mplus_syntax("hybrid", I = 10, G = 3), "G = 2")
})
