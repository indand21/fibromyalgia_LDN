priors <- data.frame(parameter = c("pbo_rel_mean", "pbo_logit_sd", "kpbo_per_day"),
                     value = c(0.12, 0.35, 0.1), sd = c(NA, NA, 0.5))
flat <- list("0" = data.frame(hour = 0:24, Cu_NTX = 0, Cu_BN = 0, RO = 0),
             "4.5" = data.frame(hour = 0:24, Cu_NTX = 1, Cu_BN = 0, RO = 0.4))
arms <- data.frame(trial = "T", arm = c("pbo", "ldn"), schedule = c("0:0", "0:4.5"), day = 84,
                   outcome = "pct_reduction", baseline_nrs = 6, mean = c(12, 25), se = c(2, 2))

test_that("parameter vectors map to valid theta for every hypothesis", {
  for (h in c("H0", "H1", "H2", "H3")) {
    p <- init_par(h, "T"); th <- theta_from(p, h, list(IC50_TLR4 = 1000, Kd_OGFr = NA))
    expect_true(th$kout > 0 && th$kpbo > 0)
    if (h != "H0") expect_true(th$Emax > 0 && th$Emax < 1)
  }
})

test_that("log posterior is finite at the initial values", {
  for (h in c("H0", "H1", "H2", "H3")) {
    p <- init_par(h, "T")
    expect_true(is.finite(log_prior(p, h, priors, "T")))
    expect_true(is.finite(log_lik_arms(p, h, arms, flat, list(IC50_TLR4 = 1000, Kd_OGFr = NA))))
  }
})

test_that("pointwise log-likelihood has one entry per arm", {
  p <- init_par("H1", "T")
  expect_length(log_lik_arms(p, "H1", arms, flat, list(IC50_TLR4 = 1000, Kd_OGFr = NA), pointwise = TRUE), 2)
})

test_that("H3 with a fixed finite Kd uses it and needs no log10_Kd_OGFr", {
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = 50)
  p <- init_par("H3", "T"); p <- p[setdiff(names(p), "log10_Kd_OGFr")]
  expect_equal(theta_from(p, "H3", fx)$Kd_OGFr, 50)
  expect_true(is.finite(log_prior(p, "H3", priors, "T")))
  expect_true(is.finite(log_lik_arms(p, "H3", arms, flat, fx)))
  p2 <- init_par("H3", "T"); p2[["log10_Kd_OGFr"]] <- 2
  expect_equal(theta_from(p2, "H3", list(IC50_TLR4 = 1000, Kd_OGFr = NA))$Kd_OGFr, 100)
})

test_that("sigma_rel widens the likelihood: lowers a perfect fit, raises a bad fit", {
  p <- init_par("H0", "T")
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  pred <- arm_outcome(simulate_arm("H0", theta_from(p, "H0", fx), parse_schedule("0:0"), 84, flat, pmax_trial(p, "T")),
                      84, "pct_reduction", 6)
  perfect <- data.frame(trial = "T", arm = "pbo", schedule = "0:0", day = 84, outcome = "pct_reduction",
                        baseline_nrs = 6, mean = pred, se = 2)
  bad <- transform(perfect, mean = pred + 30)
  lo <- p; lo[["log_sigma_rel"]] <- log(0.001); hi <- p; hi[["log_sigma_rel"]] <- log(0.3)
  expect_lt(log_lik_arms(hi, "H0", perfect, flat, fx), log_lik_arms(lo, "H0", perfect, flat, fx))
  expect_gt(log_lik_arms(hi, "H0", bad, flat, fx), log_lik_arms(lo, "H0", bad, flat, fx))
})

test_that("outcome_scale returns 100 for pct_reduction and baseline for change_nrs", {
  expect_equal(outcome_scale("pct_reduction", 6), 100)
  expect_equal(outcome_scale("change_nrs", 6.5), 6.5)
  expect_equal(outcome_scale(c("pct_reduction", "change_nrs"), 7), c(100, 7))
})

test_that("log_prior errors clearly on missing priors", {
  expect_error(log_prior(init_par("H0", "T"), "H0", priors[-1, ], "T"), "pbo_rel_mean")
})

test_that("log_lik_arms errors on a non-finite se and is finite on the finite-se subset", {
  p <- init_par("H1", "T"); fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  bad <- arms; bad$se[2] <- NA
  expect_error(log_lik_arms(p, "H1", bad, flat, fx), "non-finite se")
  expect_true(is.finite(log_lik_arms(p, "H1", bad[is.finite(bad$se), ], flat, fx)))
})
