cache <- list("0" = data.frame(hour = 0:24, Cu_NTX = 0, Cu_BN = 0, RO = 0),
              "4.5" = data.frame(hour = 0:24, Cu_NTX = 1, Cu_BN = 0, RO = 0.4))
th <- list(kout = 0.1, kpbo = 0.1, Emax = 0.5, kA = 0.2, gamma = 1)

test_that("simulated NRS stays in [0, 10] and baselines respect the trial's enrolment bounds", {
  # R2 fix round 1: the bounds are now per trial (default [0, 10]); this trial screens at NRS >= 4
  x <- simulate_individuals("H1", th, parse_schedule("0:4.5"), 84, cache, 0.2, 5000, 0.8, 0.8, 7, 3, 1.5, seed = 1,
                            base_bounds = c(4, 10))
  expect_true(all(x$pain_end >= 0 & x$pain_end <= 10)); expect_true(all(x$pain0 >= 4 & x$pain0 <= 10))
  y <- simulate_individuals("H1", th, parse_schedule("0:4.5"), 84, cache, 0.2, 5000, 0.8, 0.8, 7, 3, 1.5, seed = 1)
  expect_true(all(y$pain0 >= 0 & y$pain0 <= 10)); expect_true(any(y$pain0 < 4))   # physical bounds by default
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

# ---- R2 item 3 (review I7): the model-structure residual is shared by the two arms of a replicate trial ----
vt_fit <- function(log_sigma_rel) list(hyp = "H0", fixed = list(),
  draws = cbind(mu_p = c(-1.5, -1.4), log_tau = c(-1.5, -1.6), log_kpbo = log(0.1),
                log_kout = log(0.1), log_sigma_rel = log_sigma_rel))
vt_design <- list(trial = "T", schedule_ldn = "0:4.5", day = 14, n_ldn = 60, n_pbo = 60,
                  base_mean = 6, base_sd = 1.2)
vt_om <- list(omega_p = 0.3, omega_e = 0, sd_res = 0.5)

test_that("virtual_trials: the shared residual all but cancels in the replicate arm difference", {
  a <- virtual_trials(vt_fit(log(1e-8)), vt_design, cache, vt_om, n_rep = 300, n_draws = 2)
  b <- virtual_trials(vt_fit(log(0.25)), vt_design, cache, vt_om, n_rep = 300, n_draws = 2)
  # the shared shift multiplies each subject's own baseline, so it leaves the difference only through the
  # (small) baseline imbalance between the two arms -- not as the root-2 inflation an independent draw gave.
  indep <- sqrt(2) * 0.25 * vt_design$base_mean          # what two independent arm draws would have added
  expect_lt(abs(stats::sd(b$diff) - stats::sd(a$diff)), 0.05 * indep)
  expect_equal(stats::sd(b$diff), stats::sd(a$diff), tolerance = 0.05)
  # the absolute arm level does move with it
  expect_false(isTRUE(all.equal(b$sd_ldn, a$sd_ldn, tolerance = 1e-8)))
})

test_that("virtual_trials: each replicate trial gets its own shared residual", {
  set.seed(1)
  v <- virtual_trials(vt_fit(log(0.25)), vt_design, cache, vt_om, n_rep = 40, n_draws = 1)
  expect_equal(nrow(v), 40)
  expect_gt(stats::sd(v$diff), 0)
})

# ---- R2 item 4 (review I3), corrected in fix round 1 (IMPORTANT 1) ----
# Baselines are drawn from a truncated normal on the trial's bounds, but the REPORTED mean and SD describe the
# enrolled (already screened) population. The parent is therefore moment-matched so that the TRUNCATED
# distribution reproduces the reported moments; truncating N(reported) would apply the screen a second time.
test_that("simulated baselines reproduce the REPORTED moments within the trial's bounds, with no point mass", {
  for (tg in list(c(6.25, 1.45, 4), c(7.565, 1.605, 4), c(5.08, 1.5, 0))) {
    x <- simulate_individuals("H0", th, parse_schedule("0:0"), 84, cache, 0.1115, 200000,
                              2.0, 0, tg[1], tg[2], 0.48, seed = 1, base_bounds = c(tg[3], 10))
    expect_true(all(x$pain0 > tg[3] & x$pain0 < 10))
    expect_equal(mean(x$pain0 == tg[3]), 0); expect_equal(mean(x$pain0 == 10), 0)
    expect_equal(mean(x$pain0), tg[1], tolerance = 0.005)
    expect_equal(stats::sd(x$pain0), tg[2], tolerance = 0.01)
    bp <- attr(x, "baseline_parent")
    expect_true(bp$matched); expect_equal(c(bp$lower, bp$upper), c(tg[3], 10))
  }
})

test_that("an infeasible target falls back to the reported moments on the bounds, and says so", {
  x <- simulate_individuals("H0", th, parse_schedule("0:0"), 84, cache, 0.1115, 20000,
                            2.0, 0, 5.08, 1.5, 0.48, seed = 3, base_bounds = c(4, 10))
  bp <- attr(x, "baseline_parent")
  expect_false(bp$matched); expect_equal(c(bp$mu, bp$sigma), c(5.08, 1.5))
  expect_true(all(x$pain0 > 4 & x$pain0 < 10))
})

test_that("per-trial baseline bounds are read from the curated file and looked up by trial", {
  bt <- read_baseline_bounds(file.path("..", "..", "data", "params", "baseline_bounds.csv"))
  expect_equal(baseline_bounds_for(bt, "younger2013"), c(0.5, 10))     # floor from the fix-round-2 stability rule
  expect_equal(baseline_bounds_for(bt, "FINAL"), c(4, 10))                 # case-insensitive
  expect_equal(baseline_bounds_for(bt, "innova"), c(4, 10))
  expect_error(baseline_bounds_for(bt, "nosuchtrial"), "no baseline bounds")
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(trial = c("a", "b"), lower_nrs = c(3, 0), upper_nrs = c(10, 10),
                              source = "s", location = "l", verified = "yes"), tmp, row.names = FALSE)
  t2 <- read_baseline_bounds(tmp)
  expect_equal(baseline_bounds_for(t2, "a"), c(3, 10)); expect_equal(baseline_bounds_for(t2, "B"), c(0, 10))
  utils::write.csv(data.frame(trial = "a", lower_nrs = 11, upper_nrs = 10, source = "s", location = "l",
                              verified = "yes"), tmp, row.names = FALSE)
  expect_error(read_baseline_bounds(tmp), "lower_nrs")
})

test_that("simulate_individuals reports how often the physical pain_end clip binds", {
  x <- simulate_individuals("H0", th, parse_schedule("0:0"), 84, cache, 0.1115, 20000,
                            2.0, 0, 5.08, 1.5, 0.48, seed = 2)
  cl <- attr(x, "pain_end_clip")
  expect_true(all(c("low", "high") %in% names(cl)))
  expect_true(all(cl >= 0 & cl <= 1))
  expect_true(all(x$pain_end >= 0 & x$pain_end <= 10))
})

test_that("calibrate_omegas reads the trial's bounds from the curated table and reports the parent", {
  dr <- cbind(mu_p = -1.5, log_kpbo = log(0.1), log_kout = log(0.1)); dr <- cbind(dr, `eta_younger2013` = 0)
  fit <- list(hyp = "H0", fixed = list(), draws = dr)
  arms <- data.frame(trial = "younger2013", arm = "placebo", schedule = "0:0", day = 14, baseline_nrs = 5.08,
                     sd_pct = 27.9)
  bt <- read_baseline_bounds(file.path("..", "..", "data", "params", "baseline_bounds.csv"))
  om <- suppressWarnings(calibrate_omegas(fit, arms, cache, sd_res = 0.48, n = 500, bounds = bt))
  expect_true(startsWith(om$baseline_bounds, "[0.5, 10]"))
  expect_true(om$baseline_parent_matched)
  expect_equal(om$baseline_parent_mu, 5.08, tolerance = 0.01)
  bad <- bt; bad$trial[bad$trial == "younger2013"] <- "other"
  expect_error(suppressWarnings(calibrate_omegas(fit, arms, cache, sd_res = 0.48, n = 500, bounds = bad)), "no baseline bounds")
})

# ---- R2 fix round 1: virtual_trials uses the design's bounds and reports the pain_end clip (Minor 5) ----
test_that("virtual_trials applies design$base_bounds and reports the pain_end clip per draw and arm", {
  d4 <- modifyList(vt_design, list(base_bounds = c(4, 10)))
  v <- virtual_trials(vt_fit(log(0.05)), d4, cache, vt_om, n_rep = 20, n_draws = 2)
  cl <- attr(v, "pain_end_clip")
  expect_s3_class(cl, "data.frame")
  expect_equal(nrow(cl), 2 * 2)                               # 2 draws x 2 arms
  expect_true(all(c("draw", "arm", "low", "high") %in% names(cl)))
  expect_setequal(unique(cl$arm), c("ldn", "placebo"))
  expect_true(all(cl$low >= 0 & cl$low <= 1 & cl$high >= 0 & cl$high <= 1))
  bp <- attr(v, "baseline_parent")
  expect_equal(c(bp$lower, bp$upper), c(4, 10))
  v0 <- virtual_trials(vt_fit(log(0.05)), vt_design, cache, vt_om, n_rep = 5, n_draws = 1)
  expect_equal(c(attr(v0, "baseline_parent")$lower, attr(v0, "baseline_parent")$upper), c(0, 10))
})

# ---- R2 fix round 2, ruling 1: the Younger floor is chosen so that the calibration does not depend on the seed ----
test_that("at the curated Younger floor, calibrate_omegas' omega_p is stable across seeds", {
  dr <- cbind(mu_p = logit(0.1114), log_kpbo = log(0.106), log_kout = log(0.1)); dr <- cbind(dr, `eta_younger2013` = 0)
  fit <- list(hyp = "H0", fixed = list(), draws = dr)
  arms <- data.frame(trial = "younger2013", arm = "placebo", schedule = "0:0", day = 28, baseline_nrs = 50.8,
                     scale = "0-100", sd_pct = 27.9)
  bt <- read_baseline_bounds(file.path("..", "..", "data", "params", "baseline_bounds.csv"))
  expect_gt(baseline_bounds_for(bt, "younger2013")[1], 0)          # a floor above zero (% change undefined at 0)
  op <- vapply(c(11, 12, 13, 14, 15, 16), function(s)
    suppressWarnings(calibrate_omegas(fit, arms, cache, sd_res = 0.48, n = 4000, bounds = bt, seed = s))$omega_p, 0)
  expect_true(all(abs(op / stats::median(op) - 1) <= 0.10))
})
