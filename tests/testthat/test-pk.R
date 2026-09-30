test_that("one-compartment parent without metabolite matches the Bateman function", {
  p <- modifyList(pk_default(), list(ka = 1, fp = 0, fm = 0, CL = 50, V1 = 500, Q = 0))
  s <- simulate_pk(10, seq(0, 24, 0.5), p)
  k <- 50 / 500; dose <- mg_to_nmol(10, "NTX")
  expected <- dose / 500 * 1 / (1 - k) * (exp(-k * s$time) - exp(-1 * s$time))
  expect_equal(s$Cp_NTX, expected, tolerance = 1e-5)
})

test_that("PK is dose-linear", {
  a <- simulate_pk(5, seq(0, 48, 0.25)); b <- simulate_pk(10, seq(0, 48, 0.25))
  expect_equal(max(b$Cp_NTX) / max(a$Cp_NTX), 2, tolerance = 1e-6)
  expect_equal(max(b$Cp_BN) / max(a$Cp_BN), 2, tolerance = 1e-6)
})

test_that("unbound concentrations apply fu", {
  s <- simulate_pk(4.5, seq(0, 24, 1)); p <- pk_default()
  expect_equal(s$Cu_NTX, p$fu_NTX * s$Cp_NTX)
})

test_that("nca_summary recovers a known half-life and AUC", {
  t <- seq(0, 48, 0.25); c <- 100 * exp(-log(2) / 4 * t)
  r <- nca_summary(t, c)
  expect_equal(unname(r["thalf"]), 4, tolerance = 1e-6)
  expect_equal(unname(r["AUCinf"]), 100 * 4 / log(2), tolerance = 1e-3)
})

test_that("fit_pk recovers CL and V1 from synthetic statistics", {
  truth <- modifyList(pk_default(), list(CL = 150, V1 = 300))
  obs <- expand.grid(dose_mg = c(50, 100), analyte = "NTX", stat = c("Cmax", "AUCinf"), stringsAsFactors = FALSE)
  obs$mean <- pk_predict_stats(truth, obs); obs$sd <- 0.2 * obs$mean
  fit <- fit_pk(obs, start = pk_default(), pars = c("CL", "V1"))
  expect_equal(fit$par$CL, 150, tolerance = 0.02); expect_equal(fit$par$V1, 300, tolerance = 0.02)
})

test_that("pk_predict_stats rejects unknown stat or analyte", {
  obs <- data.frame(dose_mg = 50, analyte = "NTX", stat = "AUClast")
  expect_error(pk_predict_stats(pk_default(), obs), "AUClast")
  obs <- data.frame(dose_mg = 50, analyte = "XYZ", stat = "Cmax")
  expect_error(pk_predict_stats(pk_default(), obs), "XYZ")
})
