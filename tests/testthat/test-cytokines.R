cp <- data.frame(cytokine = c("TNF", "IL6"), quantity = "thalf", value = c(1, 2), unit = "h")
prof <- data.frame(hour = 0:24, Cu_NTX = 1, Cu_BN = 0, RO = 0.4)

test_that("H0 and H1 predict no cytokine change", {
  for (h in c("H0", "H1")) expect_equal(cytokine_prediction(h, list(Emax = 0.5), prof, cp, 56)$pct_change_pred, c(0, 0))
})

test_that("H2 steady-state change equals -100 x mean TLR4 occupancy (upper bound, full TRIF dependence)", {
  r <- cytokine_prediction("H2", list(Emax = 0.5, IC50_TLR4 = 1), prof, cp, 56)
  expect_equal(r$pct_change_pred, c(-50, -50), tolerance = 1e-6)
})

test_that("H3 signal is Emax x mean OGFr occupancy", {
  r <- cytokine_prediction("H3", list(Emax = 0.4, Kd_OGFr = 1), prof, cp, 56)
  expect_equal(r$pct_change_pred, c(-20, -20), tolerance = 1e-6)
})

test_that("half-lives in minutes and hours give the same kinetics; early days approach steady state gradually", {
  cpm <- data.frame(cytokine = "TNF", quantity = "thalf", value = 60, unit = "min")
  cph <- data.frame(cytokine = "TNF", quantity = "thalf", value = 1, unit = "h")
  th <- list(Emax = 0.5, IC50_TLR4 = 1)
  expect_equal(cytokine_prediction("H2", th, prof, cpm, 0.05)$pct_change_pred,
               cytokine_prediction("H2", th, prof, cph, 0.05)$pct_change_pred)
  expect_equal(cytokine_prediction("H2", th, prof, cph, 0.05)$pct_change_pred, -50 * (1 - exp(-log(2) * 1.2)), tolerance = 1e-8)
})

test_that("an unknown half-life unit is an error", {
  bad <- data.frame(cytokine = "TNF", quantity = "thalf", value = 1, unit = "day")
  expect_error(cytokine_prediction("H2", list(Emax = 0.5, IC50_TLR4 = 1), prof, bad, 56), "unit")
})

test_that("real curated cytokine parameters yield finite predictions, with IL-10 from the mean of bounds", {
  real <- read_dataset(file.path("..", "..", "data", "params", "cytokine_params.csv"))
  r <- cytokine_prediction("H2", list(Emax = 0.5, IC50_TLR4 = 1), prof, real, 56)
  expect_gte(nrow(r), 3)
  expect_true(all(is.finite(r$pct_change_pred)))
  expect_true(all(c("TNF-alpha", "IL-6", "IL-10") %in% r$cytokine))
  expect_equal(r$halflife_basis[r$cytokine == "IL-10"], "mean of bounds")
  expect_equal(r$halflife_basis[r$cytokine == "IL-6"], "thalf")
})

test_that("bounds average and missing half-lives error", {
  b <- data.frame(cytokine = "X", quantity = c("thalf_lower", "thalf_upper"), value = c(1, 3), unit = "h")
  th <- list(Emax = 0.5, IC50_TLR4 = 1)
  expect_equal(cytokine_prediction("H2", th, prof, b, 0.05)$pct_change_pred, -50 * (1 - exp(-log(2) / 2 * 1.2)), tolerance = 1e-8)
  expect_error(cytokine_prediction("H2", th, prof, b[1, ], 56), "no usable half-life")
})
