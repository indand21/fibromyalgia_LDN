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

test_that("exposure_profile extends the dosing run (up to max_doses) until occupancy is at steady state", {
  bp <- binding_params(0.5, 5, koff_NTX = 0.005, koff_BN = 0.005, keq = 5)
  p <- modifyList(pk_default(), list(Q = 1e-6, CL = 60, V1 = 400, CLm = 60, Vm = 800))
  expect_warning(exposure_profile(1, p, bp, n_doses = 4), "steady state")
  expect_silent(pr <- exposure_profile(1, p, bp, n_doses = 4, max_doses = 64))
  ref <- exposure_profile(1, p, bp, n_doses = 64)
  expect_equal(mean(pr$RO), mean(ref$RO), tolerance = 0.01)
  expect_equal(attr(pr, "n_doses") >= 8, TRUE)
})

# ---- R1 fix round 4 item 1: the stopping rule must certify the EXTRAPOLATED remaining gap, not one interval step ----
test_that("exposure_ss_gap extrapolates the geometric tail, and guards degenerate cases", {
  # means 0.4, 0.45, 0.475: steps 0.05 then 0.025, so r = 0.5 and the remaining gap is 0.025 * 0.5 / 0.5 = 0.025
  expect_equal(exposure_ss_gap(c(0.4, 0.45, 0.475)), 0.025)
  expect_equal(exposure_ss_gap(c(0.4, 0.49, 0.499)), 0.09 * 0.1 * 0.1 / 0.9, tolerance = 1e-9)
  # at the realistic stage-J ratio (about 0.85) the remaining gap is nearly 6 times the last step, which is the whole point
  m <- c(0.70, 0.70 + 0.01, 0.70 + 0.01 + 0.0085)
  expect_equal(exposure_ss_gap(m), 0.0085 * 0.85 / 0.15, tolerance = 1e-9)
  expect_gt(exposure_ss_gap(m), 5 * (m[3] - m[2]))
  expect_identical(exposure_ss_gap(c(0.4, 0.5, 0.65)), Inf)            # r >= 1: cannot certify
  expect_equal(exposure_ss_gap(c(0.5, 0.5, 0.5)), 0)                   # already flat
  expect_equal(exposure_ss_gap(c(0.40, 0.50, 0.4999)), 1e-4)           # non-monotone: fall back to the last step
})

test_that("a slow-koff (stage-J-like) profile matches a 240-dose reference within 0.2 percentage points", {
  bp <- binding_params(Kd_NTX = 0.83, Kd_BN = 3.72, koff_NTX = 0.044, koff_BN = 0.0042, keq = 2)
  p <- modifyList(pk_default(), list(Q = 1e-6, CL = 380, V1 = 1400, CLm = 68, Vm = 1100, fp = 0.95, fu_NTX = 0.79, fu_BN = 0.79))
  ref <- suppressWarnings(exposure_profile(4.5, p, bp, n_doses = 240))
  got <- exposure_profile(4.5, p, bp, max_doses = EXPOSURE_MAX_DOSES)
  expect_lt(abs(mean(got$RO) - mean(ref$RO)), 0.002)
  expect_lt(abs(max(got$RO) - max(ref$RO)), 0.002)
  expect_lt(abs(min(got$RO) - min(ref$RO)), 0.002)
  expect_gt(attr(got, "n_doses"), 14)                                   # the old 14-dose default was not at steady state
  expect_lte(attr(got, "ss_gap"), EXPOSURE_SS_TOL)
})

test_that("a fast-koff profile is still certified at the 14-dose default (no extra cost)", {
  bp <- binding_params(0.5, 5, 2, 2, 5)
  pr <- exposure_profile(4.5, pk_default(), bp, max_doses = EXPOSURE_MAX_DOSES)
  expect_equal(attr(pr, "n_doses"), 14)
  expect_lte(attr(pr, "ss_gap"), EXPOSURE_SS_TOL)
})

test_that("the n_doses attribute is an integer on every path, so callers can vapply over it", {
  fast <- binding_params(0.5, 5, 2, 2, 5)
  slow <- binding_params(Kd_NTX = 0.83, Kd_BN = 3.72, koff_NTX = 0.044, koff_BN = 0.0042, keq = 2)   # stage-J-like: extends
  p <- modifyList(pk_default(), list(Q = 1e-6, CL = 380, V1 = 1400, CLm = 68, Vm = 1100, fp = 0.95))
  prs <- list(zero = exposure_profile(0, p, fast), no_ext = exposure_profile(4.5, p, fast),
              capped = exposure_profile(4.5, p, fast, max_doses = EXPOSURE_MAX_DOSES),
              extended = exposure_profile(4.5, p, slow, max_doses = EXPOSURE_MAX_DOSES))
  for (nm in names(prs)) { expect_type(attr(prs[[nm]], "n_doses"), "integer"); expect_type(attr(prs[[nm]], "ss_gap"), "double") }
  expect_gt(attr(prs$extended, "n_doses"), 14)                      # the extension path really was taken
  expect_type(EXPOSURE_MAX_DOSES, "integer")
  expect_equal(vapply(prs, function(z) attr(z, "n_doses"), 0L)[c("zero", "no_ext")], c(zero = 0L, no_ext = 14L))
})
