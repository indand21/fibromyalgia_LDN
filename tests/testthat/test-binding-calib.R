# R1 (review C1, C4, M8): joint binding calibration (stage J), precision weights, keq choice, Laplace draws.
fast_pk <- function() modifyList(pk_default(), list(Q = 1e-6, CL = 60, V1 = 400, CLm = 60, Vm = 800, fu_NTX = 0.8, fu_BN = 0.8))
bp_base <- function() binding_params(Kd_NTX = 0.5, Kd_BN = 2, koff_NTX = 20, koff_BN = 2, keq = 2)
pet_rows <- function() data.frame(
  study = c(rep("Rab", 4), rep("Lee", 3)), dose_mg = c(2, 5, 15, 50, 50, 50, 50), regimen = "single",
  time_h = c(4, 4, 4, 4, 24, 48, 72), reference = c(rep("absolute", 4), rep(PET_NORM_REF, 3)),
  occupancy_pct = NA_real_, sd = NA_real_, stringsAsFactors = FALSE)

test_that("precision weights use the reported sd when finite and 10 percentage points otherwise", {
  expect_equal(pet_precision_sd(c(4.9, NA, Inf, 12)), c(4.9, 10, 10, 12))
  expect_equal(PET_DEFAULT_SD, 10)
  rows <- data.frame(occupancy_pct = c(50, 50), sd = c(5, NA))
  expect_equal(pet_wsse_from_pred(c(60, 60), rows), (10 / 5)^2 + (10 / 10)^2)
})

test_that("normalised-reference rows are scored relative to the model's own 1-hour prediction", {
  bp <- scale_binding(bp_base(), kd_scale = 5, koff_scale = 0.05)
  r <- data.frame(dose_mg = 5, regimen = "single", time_h = c(24, 24, 1), reference = c(PET_NORM_REF, "absolute", "absolute"),
                  stringsAsFactors = FALSE)
  p <- predict_pet(bp, r, fast_pk())
  expect_equal(p[1], 100 * p[2] / p[3], tolerance = 1e-6)
  expect_gt(abs(p[1] - p[2]), 20)                    # the two estimands genuinely differ in this scenario
  rows <- r[1, ]; rows$occupancy_pct <- p[1]; rows$sd <- NA
  expect_equal(pet_wsse(bp, rows, fast_pk()), 0, tolerance = 1e-10)
  rows$reference <- "absolute"                       # the same number scored as absolute occupancy is a poor fit
  expect_equal(pet_wsse(bp, rows, fast_pk()), ((p[2] - p[1]) / PET_DEFAULT_SD)^2, tolerance = 1e-6)
  expect_gt(pet_wsse(bp, rows, fast_pk()), 4)
})

test_that("predict_pet gives the same answer whether rows are predicted together or one at a time", {
  rows <- pet_rows(); bp <- scale_binding(bp_base(), 3, 0.1)
  together <- predict_pet(bp, rows, fast_pk())
  alone <- vapply(seq_len(nrow(rows)), function(i) predict_pet(bp, rows[i, ], fast_pk()), 0)
  expect_equal(together, alone, tolerance = 1e-4)
})

test_that("scale_binding scales both Kd and both koff and keeps kon = koff / Kd", {
  b <- scale_binding(bp_base(), kd_scale = 10, koff_scale = 0.1, keq = 8)
  expect_equal(c(b$Kd_NTX, b$Kd_BN, b$koff_NTX, b$koff_BN, b$keq), c(5, 20, 2, 0.2, 8))
  expect_equal(b$kon_NTX, b$koff_NTX / b$Kd_NTX); expect_equal(b$kon_BN, b$koff_BN / b$Kd_BN)
  expect_equal(c(b$kd_scale, b$koff_scale), c(10, 0.1))
})

test_that("the joint fit recovers known Kd and koff scales from synthetic PET rows within 10%", {
  rows <- pet_rows(); pk <- fast_pk()
  rows$occupancy_pct <- predict_pet(scale_binding(bp_base(), 10, 0.05), rows, pk)
  fit <- fit_binding_joint(bp_base(), rows, pk, keq = 2, start = c(log(3), log(0.2)))
  expect_equal(fit$kd_scale, 10, tolerance = 0.10)
  expect_equal(fit$koff_scale, 0.05, tolerance = 0.10)
  expect_lt(fit$wsse, 1e-3)
  expect_named(fit$se_log, c("log_kd_scale", "log_koff_scale"))
  expect_equal(fit$bp$Kd_NTX, 0.5 * fit$kd_scale)
  expect_equal(fit$keq, 2)
})

test_that("Laplace draws have the requested moments and their quantiles bracket the point estimate", {
  mode <- c(a = log(10), b = log(0.05)); cov <- matrix(c(0.04, 0.01, 0.01, 0.09), 2)
  d <- laplace_draws(mode, cov, n = 4000, seed = 42)
  expect_equal(dim(d), c(4000L, 2L)); expect_equal(colnames(d), c("a", "b"))
  expect_equal(unname(colMeans(d)), unname(mode), tolerance = 0.02)
  expect_equal(unname(cov(d)), cov, tolerance = 0.1)
  f <- function(x) exp(x[["a"]]) * exp(x[["b"]])
  q <- summarise_draws(apply(d, 1, f))
  # R1 fix round 4: `n` is the number of NON-missing draws the quantiles are based on, not the number requested
  expect_named(q, c("median", "lo95", "hi95", "n"))
  expect_equal(q[["n"]], 4000)
  expect_equal(summarise_draws(c(1, 2, 3, NA, NaN))[["n"]], 3)
  expect_equal(summarise_draws(c(1, 2, 3, NA))[["median"]], 2)
  expect_lt(q[["lo95"]], f(mode)); expect_gt(q[["hi95"]], f(mode))
  expect_identical(laplace_draws(mode, cov, 10, seed = 1), laplace_draws(mode, cov, 10, seed = 1))
  expect_error(laplace_draws(mode, diag(c(1, -1)), 10), "positive definite")
})

test_that("occupancy_ed50 returns the single dose giving 50% occupancy at the scan time", {
  bp <- scale_binding(bp_base(), 10, 1); pk <- fast_pk()
  ed <- occupancy_ed50(bp, pk, time_h = 4)
  r <- data.frame(dose_mg = ed, regimen = "single", time_h = 4, reference = "absolute")
  expect_equal(predict_pet(bp, r, pk), 50, tolerance = 1e-3)
})

