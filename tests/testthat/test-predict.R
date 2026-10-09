test_that("lpd_normal equals log of mean density", {
  expect_equal(lpd_normal(0, c(0, 0), 1), stats::dnorm(0, 0, 1, log = TRUE))
  expect_equal(lpd_normal(1, c(0, 2), 1), log(mean(stats::dnorm(1, c(0, 2), 1))))
})

pred_fixture <- function(log_sigma_rel = NULL) {
  d <- cbind(mu_p = logit(0.1), log_tau = log(0.3), log_kpbo = log(0.1), log_kout = log(0.1))
  if (!is.null(log_sigma_rel)) d <- cbind(d, log_sigma_rel = log_sigma_rel)
  list(hyp = "H0", fixed = list(), draws = d)
}
pred_cache <- list("0" = data.frame(hour = 0:24, Cu_NTX = 0, Cu_BN = 0, RO = 0))
pred_rows <- data.frame(trial = "X", arm = "pbo", schedule = "0:0", day = 84, endpoint = "change_nrs", baseline_nrs = 6)

test_that("predict_trial draws a new-trial placebo effect and returns one prediction per draw", {
  p <- predict_trial(pred_fixture(), pred_rows, pred_cache, n_draws = 50, seed = 1)
  expect_equal(nrow(p), 50); expect_true(sd(p$pred) > 0)
})

test_that("predict_trial adds the model-structure residual when log_sigma_rel is present", {
  p0 <- predict_trial(pred_fixture(-Inf), pred_rows, pred_cache, n_draws = 200, seed = 1)
  p1 <- predict_trial(pred_fixture(log(0.1)), pred_rows, pred_cache, n_draws = 200, seed = 1)
  expect_gt(sd(p1$pred), sd(p0$pred))
  # residual SD is on the change_nrs scale: sigma_rel * baseline_nrs
  expect_equal(sd(p1$pred - p0$pred), 0.1 * 6, tolerance = 0.2)
})

# ---- R1 fix round 4 item 4: predict_pet must honour the curated dose count ----
test_that("predict_pet uses the n_doses column when present and the documented default otherwise", {
  bp <- binding_params(0.8, 3.7, 0.044, 0.0042, 2)
  p <- modifyList(pk_default(), list(Q = 1e-6, CL = 380, V1 = 1400, CLm = 68, Vm = 1100, fp = 0.95))
  r <- data.frame(dose_mg = 50, regimen = "daily", time_h = 24, reference = "absolute", n_doses = c(4, 7))
  pr <- predict_pet(bp, r, p)
  expect_lt(pr[1], pr[2])                                          # fewer doses, less accumulation
  expect_equal(pr[2], predict_pet(bp, r[2, setdiff(names(r), "n_doses")], p))   # default for a daily row is 7
  s <- data.frame(dose_mg = 50, regimen = "single", time_h = 4, reference = "absolute")
  expect_equal(predict_pet(bp, s, p), predict_pet(bp, cbind(s, n_doses = 1), p))
  expect_equal(unname(PET_DEFAULT_N_DOSES[c("single", "daily")]), c(1, 7))
})

# ---- R2 item 3 (review I7): one shared model-structure residual per trial per draw ----
two_arm_rows <- data.frame(trial = "X", arm = c("ldn", "pbo"), schedule = c("0:0", "0:0"), day = 84,
                           endpoint = "change_nrs", baseline_nrs = c(6, 6))

test_that("the shared residual cancels in the difference between two arms of the same trial", {
  lo <- predict_trial(pred_fixture(log(1e-6)), two_arm_rows, pred_cache, n_draws = 400, seed = 3)
  hi <- predict_trial(pred_fixture(log(0.3)), two_arm_rows, pred_cache, n_draws = 400, seed = 3)
  dlo <- lo$pred[lo$row == 1] - lo$pred[lo$row == 2]
  dhi <- hi$pred[hi$row == 1] - hi$pred[hi$row == 2]
  expect_equal(dhi, dlo)                                   # identical: the shared term is subtracted away
  # while each arm's own level still carries it
  expect_gt(stats::sd(hi$pred[hi$row == 1]), stats::sd(lo$pred[lo$row == 1]))
})

test_that("a single arm's predictive variance still carries the full residual", {
  p0 <- predict_trial(pred_fixture(-Inf), two_arm_rows, pred_cache, n_draws = 3000, seed = 5)
  p1 <- predict_trial(pred_fixture(log(0.1)), two_arm_rows, pred_cache, n_draws = 3000, seed = 5)
  for (rw in 1:2)
    expect_equal(stats::sd(p1$pred[p1$row == rw] - p0$pred[p0$row == rw]), 0.1 * 6, tolerance = 0.1)
})

test_that("two different trials get independent shared residuals", {
  a <- predict_trial(pred_fixture(log(0.2)), two_arm_rows[1, ], pred_cache, n_draws = 500, seed = 11)
  b <- predict_trial(pred_fixture(log(0.2)), two_arm_rows[1, ], pred_cache, n_draws = 500, seed = 12)
  expect_false(isTRUE(all.equal(a$pred, b$pred)))
  # and within one call the per-draw shared term varies from draw to draw
  expect_gt(stats::sd(a$pred), 0)
})

# ---- R2 item 6 (review I8): Monte Carlo standard error of a held-out elpd ----
test_that("elpd_mcse reproduces the single-observation delta-method formula and the elpd itself", {
  set.seed(8); draws <- stats::rnorm(4000, -0.4, 0.6)
  r <- elpd_mcse(obs = -0.34, pred = cbind(draws), se = 0.21)
  w <- stats::dnorm(-0.34, draws, 0.21)
  expect_equal(r$elpd, lpd_normal(-0.34, draws, 0.21))
  expect_equal(r$mcse, stats::sd(w) / (mean(w) * sqrt(length(draws))))
  expect_equal(r$n_draws, 4000); expect_equal(r$n_obs, 1)
})

test_that("elpd_mcse matches a bootstrap over draws and falls as 1/sqrt(S)", {
  set.seed(9)
  d1 <- stats::rnorm(2000, -0.4, 0.6); d2 <- stats::rnorm(2000, 0.2, 0.7)
  pm <- cbind(d1, d2); obs <- c(-0.34, 0.31); se <- c(0.21, 0.26)
  r <- elpd_mcse(obs, pm, se)
  bs <- replicate(400, { i <- sample.int(2000, replace = TRUE); elpd_mcse(obs, pm[i, ], se)$elpd })
  expect_equal(r$mcse, stats::sd(bs), tolerance = 0.25)
  r4 <- elpd_mcse(obs, pm[rep(seq_len(2000), 4), ], se)
  expect_equal(r4$mcse, r$mcse / 2, tolerance = 0.05)
})

test_that("elpd_mcse accounts for draws shared between observations", {
  set.seed(10); d <- stats::rnorm(3000, -0.4, 0.6)
  shared <- elpd_mcse(c(-0.34, -0.34), cbind(d, d), c(0.21, 0.21))
  one <- elpd_mcse(-0.34, cbind(d), 0.21)
  expect_equal(shared$mcse, 2 * one$mcse)        # perfectly correlated terms add linearly
  expect_equal(shared$elpd, 2 * one$elpd)
})

test_that("elpd_mcse validates its inputs", {
  expect_error(elpd_mcse(c(1, 2), cbind(stats::rnorm(10)), 1), "obs")
  expect_error(elpd_mcse(1, cbind(stats::rnorm(1)), 1), "draw")
})

test_that("elpd_mcse accepts a plain draw vector, a matrix and a list of draw vectors", {
  set.seed(11); d <- stats::rnorm(500, -0.4, 0.6)
  a <- elpd_mcse(-0.34, d, 0.21)
  b <- elpd_mcse(-0.34, cbind(d), 0.21)
  cc <- elpd_mcse(-0.34, list(d), 0.21)
  expect_equal(a$mcse, b$mcse); expect_equal(a$mcse, cc$mcse); expect_equal(a$n_obs, 1)
})

# ---- R2 fix round 1, Minor 3: contrast MCSE from the per-draw PAIRED difference ----
test_that("elpd_mcse_diff is zero for identical predictions and matches a paired bootstrap", {
  set.seed(12); S <- 3000
  common <- stats::rnorm(S)                                   # shared predictive seed -> correlated draws
  pa <- cbind(-0.4 + 0.6 * common, 0.1 + 0.7 * stats::rnorm(S))
  pb <- cbind(-0.3 + 0.6 * common + 0.1 * stats::rnorm(S), 0.2 + 0.7 * stats::rnorm(S))
  obs <- c(-0.34, 0.31); se <- c(0.31, 0.42)
  expect_equal(elpd_mcse_diff(obs, pa, pa, se)$mcse, 0)
  r <- elpd_mcse_diff(obs, pa, pb, se)
  expect_equal(r$diff, elpd_mcse(obs, pa, se)$elpd - elpd_mcse(obs, pb, se)$elpd)
  bs <- replicate(400, { i <- sample.int(S, replace = TRUE)
    elpd_mcse(obs, pa[i, ], se)$elpd - elpd_mcse(obs, pb[i, ], se)$elpd })
  expect_equal(r$mcse, stats::sd(bs), tolerance = 0.25)
  # positively correlated draws: the paired MCSE is below the independence formula
  indep <- sqrt(elpd_mcse(obs, pa, se)$mcse^2 + elpd_mcse(obs, pb, se)$mcse^2)
  expect_lt(r$mcse, indep)
  expect_error(elpd_mcse_diff(obs, pa, pb[-1, ], se), "same number of draws")
})

# ---- R2 fix round 1, Minor 6: unequal baselines (FINAL 6.3 vs 6.2, INNOVA 7.66 vs 7.47) ----
test_that("with unequal arm baselines the shared residual scales with each arm and leaves only z*sigma*(b1 - b2)", {
  rows <- data.frame(trial = "X", arm = c("ldn", "pbo"), schedule = "0:0", day = 84, endpoint = "change_nrs",
                     baseline_nrs = c(7.66, 7.47))
  lo <- predict_trial(pred_fixture(-Inf), rows, pred_cache, n_draws = 300, seed = 4)
  hi <- predict_trial(pred_fixture(log(0.3)), rows, pred_cache, n_draws = 300, seed = 4)
  s1 <- hi$pred[hi$row == 1] - lo$pred[lo$row == 1]; s2 <- hi$pred[hi$row == 2] - lo$pred[lo$row == 2]
  expect_equal(s1 / 7.66, s2 / 7.47)                        # ONE z per draw, scaled by each arm's own baseline
  ddiff <- (hi$pred[hi$row == 1] - hi$pred[hi$row == 2]) - (lo$pred[lo$row == 1] - lo$pred[lo$row == 2])
  expect_equal(ddiff, s1 / 7.66 * (7.66 - 7.47))           # what survives in the difference: z*sigma*(b1 - b2)
  expect_lt(stats::sd(ddiff), 0.3 * 0.25)                    # a small fraction of sigma*baseline
})

# ---- R3a item 1 (controller ruling): held-out intervals for the OBSERVED value include the observation error ----
test_that("heldout_intervals: the predictive interval is wider than the true-value interval by the obs_se amount", {
  set.seed(21); draws <- stats::rnorm(40000, -0.2, 0.3)
  iv <- heldout_intervals(obs = 0.1, draws = draws, obs_se = 0.4, seed = 5)
  expect_equal(c(iv$true_lo90, iv$true_hi90), unname(stats::quantile(draws, c(0.05, 0.95))))
  w_true <- iv$true_hi90 - iv$true_lo90; w_pred <- iv$pred_hi90 - iv$pred_lo90
  expect_gt(w_pred, w_true)
  expect_equal(w_true, 2 * stats::qnorm(0.95) * 0.3, tolerance = 0.02)
  expect_equal(w_pred, 2 * stats::qnorm(0.95) * sqrt(0.3^2 + 0.4^2), tolerance = 0.02)   # sd adds in quadrature
  expect_true(nzchar(iv$interval_definition))
  # reproducible from the seed, and it leaves the caller's RNG stream alone
  set.seed(1); a <- stats::runif(1); set.seed(1)
  iv2 <- heldout_intervals(obs = 0.1, draws = draws, obs_se = 0.4, seed = 5)
  expect_identical(stats::runif(1), a)
  expect_identical(iv, iv2)
})

test_that("heldout_intervals: inside90 is judged against the PREDICTIVE interval", {
  set.seed(22); draws <- stats::rnorm(20000, 0, 0.02)          # a near point mass for the true difference
  iv <- heldout_intervals(obs = -0.34, draws = draws, obs_se = 0.31, seed = 1)
  expect_true(-0.34 < iv$true_lo90)                              # outside the true-value interval ...
  expect_true(iv$inside90)                                       # ... but inside the predictive one
  far <- heldout_intervals(obs = 3, draws = draws, obs_se = 0.31, seed = 1)
  expect_false(far$inside90)
})

test_that("heldout_intervals: without an observation SE only the true-value interval is formed, and it says so", {
  draws <- stats::rnorm(1000)
  for (se in c(NA_real_, 0)) {
    iv <- heldout_intervals(obs = 0.2, draws = draws, obs_se = se, seed = 1)
    expect_true(is.na(iv$pred_lo90) && is.na(iv$pred_hi90) && is.na(iv$inside90))
    expect_true(is.finite(iv$true_lo90) && is.finite(iv$true_hi90))
    expect_match(iv$interval_definition, "no observation SE")
  }
})
