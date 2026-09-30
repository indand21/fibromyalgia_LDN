test_that("stage selection prefers the stage with the better held-out RMSE, not the better overall or calibration fit", {
  obs_ho <- c(90, 80, 46, 30, 95); obs_cal <- c(27, 45, 58, 61, 91, 97)
  # A2 fits the calibration rows perfectly but is worse on held-out rows; A is better held-out
  A_ho <- obs_ho + c(-5, 10, -8, 5, 3); A2_ho <- obs_ho + c(-50, -60, -40, -25, -20)
  A_cal <- obs_cal + 60; A2_cal <- obs_cal
  d <- select_binding_stage(obs_ho, A_ho, A2_ho, c(obs_ho, obs_cal), c(A_ho, A_cal), c(A2_ho, A2_cal))
  expect_equal(d$stage_adopted, "A")
  expect_lt(d$stage_A_rmse_heldout, d$stage_A2_rmse_heldout)
  expect_lt(d$stage_A2_rmse_all, d$stage_A_rmse_all)   # A2 "wins" overall yet is NOT adopted
  expect_equal(c(d$n_heldout, d$n_all), c(5, 11))
  expect_match(d$reason, "stage A kept")
})

test_that("stage selection adopts A2 only when it improves held-out RMSE", {
  obs <- c(10, 20, 30)
  d <- select_binding_stage(obs, obs + 20, obs + 2, obs, obs + 20, obs + 2)
  expect_equal(d$stage_adopted, "A2")
  d0 <- select_binding_stage(obs, obs + 5, obs + 5, obs, obs + 5, obs + 5)   # tie: keep A
  expect_equal(d0$stage_adopted, "A")
})
