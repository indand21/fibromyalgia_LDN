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
