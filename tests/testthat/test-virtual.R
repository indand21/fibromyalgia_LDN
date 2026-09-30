cache <- list("0" = data.frame(hour = 0:24, Cu_NTX = 0, Cu_BN = 0, RO = 0),
              "4.5" = data.frame(hour = 0:24, Cu_NTX = 1, Cu_BN = 0, RO = 0.4))
th <- list(kout = 0.1, kpbo = 0.1, Emax = 0.5, kA = 0.2, gamma = 1)

test_that("simulated NRS stays in [0, 10] and baselines respect enrolment bounds", {
  x <- simulate_individuals("H1", th, parse_schedule("0:4.5"), 84, cache, 0.2, 5000, 0.8, 0.8, 7, 3, 1.5, seed = 1)
  expect_true(all(x$pain_end >= 0 & x$pain_end <= 10)); expect_true(all(x$pain0 >= 4 & x$pain0 <= 10))
})

test_that("responder rates use each subject's own baseline", {
  r <- responder_rates(c(10, 5, 8), c(7, 3.4, 8))
  expect_equal(unname(r), c(2/3, 2/3, 0))
})

test_that("zero variability gives identical subjects", {
  x <- simulate_individuals("H1", th, parse_schedule("0:4.5"), 84, cache, 0.2, 10, 0, 0, 6, 0, 0, seed = 1)
  expect_equal(length(unique(round(x$pain_end, 10))), 1)
})

test_that("arm_shift moves pain_end by pain0 times the shift", {
  a <- simulate_individuals("H0", th, parse_schedule("0:4.5"), 84, cache, 0.2, 10, 0, 0, 6, 0, 0, seed = 1)
  b <- simulate_individuals("H0", th, parse_schedule("0:4.5"), 84, cache, 0.2, 10, 0, 0, 6, 0, 0, seed = 1,
                            arm_shift = 0.1)
  expect_equal(b$pain_end - a$pain_end, a$pain0 * 0.1)
})

test_that("a vector pmax_trial gives each subject its own placebo ceiling", {
  x <- simulate_individuals("H0", th, parse_schedule("0:0"), 84, cache, c(rep(0.1, 5), rep(0.4, 5)), 10,
                            0, 0, 6, 0, 0, seed = 1)
  expect_equal(length(unique(round(x$pain_end[1:5], 10))), 1)
  expect_equal(length(unique(round(x$pain_end[6:10], 10))), 1)
  expect_lt(x$pain_end[6], x$pain_end[1])
})

test_that("virtual_trials returns one row per draw and replicate with the expected columns", {
  fit <- list(hyp = "H0", fixed = list(),
              draws = cbind(mu_p = c(-1.5, -1.4), log_tau = c(-1.5, -1.6), log_kpbo = log(0.1),
                            log_kout = log(0.1), log_sigma_rel = log(0.05)))
  design <- list(trial = "T", schedule_ldn = "0:4.5", day = 14, n_ldn = 20, n_pbo = 20, base_mean = 6, base_sd = 1)
  om <- list(omega_p = 0.5, omega_e = 0, sd_res = 0.8)
  v <- virtual_trials(fit, design, cache, om, n_rep = 5, n_draws = 2)
  expect_equal(nrow(v), 10)
  expect_true(all(c("draw", "rep", "diff", "sd_ldn", "sd_pbo", "resp15_ldn", "resp15_pbo", "resp30_ldn",
                    "resp30_pbo", "resp50_ldn", "resp50_pbo") %in% names(v)))
  expect_true(all(is.finite(v$diff)))
})

test_that("calibrate_omegas stops clearly when sd_pct is missing for the needed arms", {
  fit <- list(hyp = "H1", fixed = list(), draws = cbind(mu_p = 0))
  arms <- data.frame(trial = "Younger 2013", arm = "placebo", sd_pct = NA_real_)
  expect_error(calibrate_omegas(fit, arms, cache), "sd_pct is missing")
  expect_error(calibrate_omegas(fit, arms[, c("trial", "arm")], cache), "sd_pct column")
})

test_that("calibrate_omegas warns and flags at_bound when the target SD is unattainable", {
  dr <- cbind(mu_p = -1.5, log_kpbo = log(0.1), log_kout = log(0.1)); dr <- cbind(dr, `eta_Younger 2013` = 0)
  fit <- list(hyp = "H0", fixed = list(), draws = dr)
  arms <- data.frame(trial = "Younger 2013", arm = "placebo", schedule = "0:0", day = 14, baseline_nrs = 6,
                     sd_pct = 0.1)
  expect_warning(om <- calibrate_omegas(fit, arms, cache, sd_res = 2, n = 500), "did not match")
  expect_true(om$at_bound)
  expect_equal(om$target_sd_pbo, 0.1); expect_gt(om$achieved_sd_pbo, 1.1)
  expect_equal(om$base_sd_source, "pbo:assumed 1.5")
})

test_that("calibrate_omegas warns when several placebo rows qualify", {
  dr <- cbind(mu_p = -1.5, log_kpbo = log(0.1), log_kout = log(0.1)); dr <- cbind(dr, `eta_Younger 2013` = 0)
  fit <- list(hyp = "H0", fixed = list(), draws = dr)
  arms <- data.frame(trial = "Younger 2013", arm = "placebo", schedule = "0:0", day = 14, baseline_nrs = 6,
                     baseline_sd = 1.2, sd_pct = c(30, 30))
  w <- testthat::capture_warnings(om <- calibrate_omegas(fit, arms, cache, sd_res = 0.8, n = 500))
  expect_true(any(grepl("more than one Younger 2013 placebo", w)))
  expect_equal(om$base_sd_source, "pbo:data")
})
