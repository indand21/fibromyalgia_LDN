test_that("placebo exposure is all zeros", {
  bp <- binding_params(0.5, 5, 2, 2, 5)
  pr <- exposure_profile(0, pk_default(), bp)
  expect_true(all(pr$RO == 0 & pr$Cu_NTX == 0))
})

test_that("exposure_profile returns one 24-h interval at steady state", {
  bp <- binding_params(0.5, 5, 2, 2, 5)
  pr <- exposure_profile(4.5, pk_default(), bp)
  expect_equal(range(pr$hour), c(0, 24)); expect_true(max(pr$RO) > 0)
})

test_that("exposure_profile warns when binding is not at steady state", {
  bp <- binding_params(0.5, 5, koff_NTX = 0.001, koff_BN = 0.001, keq = 5)
  expect_warning(exposure_profile(50, pk_default(), bp, n_doses = 3), "steady state")
})

test_that("profile_means averages over the dosing interval only once", {
  pr <- data.frame(hour = c(0, 12, 24), Cu_NTX = c(1, 1, 99), Cu_BN = 0, RO = c(0.2, 0.4, 0.99))
  m <- profile_means(pr, gamma = 1, IC50_TLR4 = 1, Kd_OGFr = 1)
  expect_equal(m$mRO, 0.3); expect_equal(m$free_term, 0.7); expect_equal(m$tlr4, 0.5)
})
