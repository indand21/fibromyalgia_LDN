# R1 (review C1, M8): the binary A-vs-A2 held-out RMSE rule is replaced by a documented adoption of the joint
# calibration (stage J) with a keq sensitivity. The two tests of the old select_binding_stage() rule were removed because
# that rule is what the review found wrong (unweighted RMSE on the Lee duration rows selected for kinetics, not potency).

test_that("keq choice keeps the default unless another value lowers the weighted SSE by more than 3.84", {
  tab <- data.frame(keq = c(0.5, 2, 8), wsse = c(10, 12, 11))
  expect_equal(select_keq(tab)$keq, 2)                       # 12 - 10 = 2 < 3.84
  tab$wsse <- c(5, 12, 11)
  s <- select_keq(tab); expect_equal(s$keq, 0.5); expect_match(s$reason, "3.84")
  expect_equal(select_keq(data.frame(keq = c(0.5, 2, 8), wsse = c(9, 12, 8.1)))$keq, 8)
})

test_that("binding_stage_decision adopts J and records A held-out, A2, J and keq results with a reason", {
  keq_tab <- data.frame(keq = c(0.5, 2, 8), wsse = c(30, 28, 29), kd_scale = 10, koff_scale = 0.05)
  d <- binding_stage_decision(
    A = list(rmse = 31.4, wsse = 400, n = 11, ed50 = 0.37),
    A2 = list(kd_scale = 14.9, rmse_calib = 7.2, rmse_all = 31.6, wsse = 300, ed50 = 5.4),
    J = list(kd_scale = 12, kd_scale_se_log = 0.1, koff_scale = 0.05, koff_scale_se_log = 0.3, wsse = 28, rmse = 8, ed50 = 5.5,
             hessian_pd = TRUE),
    keq_tab = keq_tab, keq_choice = select_keq(keq_tab), ed50_obs = c(value = 5.6, lo95 = 3.65, hi95 = 7.54))
  expect_equal(nrow(d), 1); expect_equal(d$stage_adopted, "J")
  for (col in c("stage_A_rmse_heldout", "stage_A_wsse_heldout", "stage_A2_kd_scale", "stage_J_kd_scale", "stage_J_koff_scale",
                "keq_adopted", "keq_wsse_0.5", "keq_wsse_2", "keq_wsse_8", "reason"))
    expect_true(col %in% names(d), info = col)
  expect_equal(d$keq_adopted, 2)
  expect_match(d$reason, "post-hoc"); expect_match(d$reason, "held-out")
})
