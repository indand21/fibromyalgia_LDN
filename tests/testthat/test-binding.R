const_pk <- function(c_ntx, c_bn, tmax = 200) data.frame(time = seq(0, tmax, 0.5), Cu_NTX = c_ntx, Cu_BN = c_bn)

test_that("single ligand at C = Kd reaches 50% occupancy", {
  bp <- binding_params(Kd_NTX = 0.5, Kd_BN = 5, koff_NTX = 2, koff_BN = 2, keq = 5)
  o <- simulate_occupancy(const_pk(0.5, 0), bp)
  expect_equal(tail(o$RO, 1), 0.5, tolerance = 1e-3)
})

test_that("two competing ligands reach the analytic equilibrium", {
  bp <- binding_params(0.5, 5, 2, 2, 5)
  o <- simulate_occupancy(const_pk(0.3, 4), bp)
  expect_equal(tail(o$RO, 1), equilibrium_occupancy(0.3, 4, 0.5, 5), tolerance = 1e-3)
})

test_that("occupancy is bounded in [0, 1]", {
  bp <- binding_params(0.5, 5, 2, 2, 5)
  o <- simulate_occupancy(const_pk(1e4, 1e4), bp)
  expect_true(all(o$RO >= -1e-9 & o$RO <= 1 + 1e-9))
})
