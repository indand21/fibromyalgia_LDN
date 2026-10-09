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

# ---- R1 (review C3): Tmax is in the PK objective ----
test_that("fit_pk includes Tmax in the objective (with a CV floor of 0.25) and uses it to identify ka", {
  truth <- modifyList(pk_default(), list(ka = 1.5, CL = 30, V1 = 400, Q = 1e-6))   # k = 0.075/h << ka: thalf carries no ka information
  obs <- data.frame(dose_mg = 100, analyte = c("NTX", "NTX"), stat = c("Tmax", "thalf"), stringsAsFactors = FALSE)
  obs$mean <- pk_predict_stats(truth, obs); obs$sd <- c(NA, 0.1 * obs$mean[2])
  fit <- fit_pk(obs, start = modifyList(truth, list(ka = 0.5)), pars = "ka")
  expect_true("Tmax" %in% fit$obs$stat)
  expect_equal(fit$cv[fit$obs$stat == "Tmax"], 0.25)
  expect_equal(fit$par$ka, 1.5, tolerance = 0.05)
})

test_that("predicted Tmax is continuous in ka (not stuck on the simulation grid)", {
  obs <- data.frame(dose_mg = 100, analyte = "NTX", stat = "Tmax", stringsAsFactors = FALSE)
  t1 <- pk_predict_stats(modifyList(pk_default(), list(ka = 1.00)), obs)
  t2 <- pk_predict_stats(modifyList(pk_default(), list(ka = 1.02)), obs)
  expect_lt(t2, t1)
})

test_that("pk_hessian_summary reports eigenvalues, positive definiteness and SEs", {
  h <- diag(c(4, 2)); dimnames(h) <- list(c("a", "b"), c("a", "b"))
  s <- pk_hessian_summary(list(hessian = h, pars = c("a", "b")))
  expect_true(s$pd); expect_equal(sort(s$eigen), c(1, 2)); expect_equal(unname(s$se), sqrt(c(0.5, 1)))
  s2 <- pk_hessian_summary(list(hessian = diag(c(4, -2)), pars = c("a", "b")))
  expect_false(s2$pd); expect_true(all(is.na(s2$se)))
})

# ---- R1 fix round 1: pk_predict_stats honours the dosing regimen (Verebey 1976 chronic rows are once-daily dosing) ----
test_that("daily-dosing Cmax exceeds single-dose Cmax for an accumulating analyte", {
  p <- modifyList(pk_default(), list(Q = 1e-6, CL = 30, V1 = 400, CLm = 20, Vm = 600))   # t1/2 ~9 h (NTX), ~21 h (BN)
  obs <- data.frame(dose_mg = 100, analyte = "BN", stat = "Cmax", regimen = c("single", "daily"), stringsAsFactors = FALSE)
  pr <- pk_predict_stats(p, obs)
  expect_gt(pr[2], pr[1] * 1.2)
  # independent check by linear superposition of the single-dose curve: C_n(t) = sum_j C_1(t + 24 j) over the n doses used
  n <- pk_n_daily_doses(p)
  s <- simulate_pk(100, seq(0, 24 * (n + 4), 0.02), p); c1 <- nM_to_ngml(s$Cp_BN, "BN")
  tt <- seq(0, 24, 0.02); css <- vapply(tt, function(t) sum(stats::approx(s$time, c1, t + 24 * (0:(n - 1)))$y), 0)
  expect_equal(pr[2], unname(refine_peak(tt, css)["Cmax"]), tolerance = 1e-5)
})

# ---- R1 fix round 4 item 3b: the superposition path must reproduce the event-driven multi-dose solve ----
test_that("superposition and the event-driven multi-dose solve agree to 1e-5 relative on every daily statistic", {
  for (p in list(modifyList(pk_default(), list(Q = 1e-6, CL = 30, V1 = 400, CLm = 20, Vm = 600)),
                 modifyList(pk_default(), list(fp = 0.95, CL = 19.2, V1 = 71.3, Q = 15.5, V2 = 104.6,
                                               CLm = 68.2, Vm = 1100, ka = 1.86)))) {
    obs <- expand.grid(dose_mg = 100, analyte = c("NTX", "BN"), stat = c("Cmax", "Tmax", "thalf", "C24", "C72"),
                       regimen = c("single", "daily"), stringsAsFactors = FALSE)
    a <- pk_predict_stats(p, obs, superposition = TRUE)
    b <- pk_predict_stats(p, obs, superposition = FALSE)
    expect_equal(a, b, tolerance = 1e-5)
  }
})

test_that("the superposition path costs less than the event-driven path on the real training rows", {
  p <- modifyList(pk_default(), list(fp = 0.95, CL = 19.2, V1 = 71.3, Q = 15.5, V2 = 104.6, CLm = 68.2, Vm = 1100, ka = 1.86))
  obs <- data.frame(dose_mg = 100, analyte = c("NTX", "BN", "NTX", "BN"), stat = c("Cmax", "Cmax", "Cmax", "Cmax"),
                    regimen = c("single", "single", "daily", "daily"), stringsAsFactors = FALSE)
  el <- function(sp) { t0 <- Sys.time(); for (i in 1:3) pk_predict_stats(p, obs, superposition = sp)
                       as.numeric(Sys.time() - t0, units = "secs") }
  expect_lt(el(TRUE), el(FALSE))
})

test_that("daily and single rows at the same dose now give different predictions, but the same half-life", {
  obs <- data.frame(dose_mg = 100, analyte = c("NTX", "NTX", "NTX", "NTX"), stat = c("Cmax", "Cmax", "thalf", "thalf"),
                    regimen = c("single", "daily", "single", "daily"), stringsAsFactors = FALSE)
  p <- modifyList(pk_default(), list(Q = 1e-6, CL = 30, V1 = 400))
  pr <- pk_predict_stats(p, obs)
  expect_false(isTRUE(all.equal(pr[1], pr[2])))
  expect_equal(pr[3], pr[4], tolerance = 1e-3)             # terminal decay after the last dose: same linear kinetics
})

test_that("an unknown regimen errors; a missing regimen column means single dose", {
  obs <- data.frame(dose_mg = 100, analyte = "NTX", stat = "Cmax", regimen = "weekly", stringsAsFactors = FALSE)
  expect_error(pk_predict_stats(pk_default(), obs), "weekly")
  o1 <- data.frame(dose_mg = 100, analyte = "NTX", stat = "Cmax", stringsAsFactors = FALSE)
  expect_equal(pk_predict_stats(pk_default(), o1), pk_predict_stats(pk_default(), transform(o1, regimen = "single")))
})

# ---- R1 fix round 2: C24 (concentration 24 h after the last dose) ----
test_that("C24 equals the simulated concentration at 24 h (analytic one-compartment case, single and daily)", {
  p <- modifyList(pk_default(), list(ka = 1, fp = 0, fm = 0, CL = 50, V1 = 500, Q = 1e-9))
  k <- 50 / 500; dose <- mg_to_nmol(100, "NTX")
  bateman <- function(t) dose / 500 * 1 / (1 - k) * (exp(-k * t) - exp(-1 * t))
  obs <- data.frame(dose_mg = 100, analyte = "NTX", stat = "C24", regimen = c("single", "daily"), stringsAsFactors = FALSE)
  pr <- pk_predict_stats(p, obs)
  expect_equal(pr[1], nM_to_ngml(bateman(24), "NTX"), tolerance = 1e-4)
  css24 <- sum(bateman(24 * (1:60)))                         # steady state, 24 h after a dose = sum over all past doses
  expect_equal(pr[2], nM_to_ngml(css24, "NTX"), tolerance = 1e-3)
})

test_that("C24 is in the objective with the concentration CV floor", {
  obs <- data.frame(dose_mg = 100, analyte = c("NTX", "NTX"), stat = c("C24", "Tmax"), regimen = "single", mean = c(2, 1),
                    sd = c(NA, NA), stringsAsFactors = FALSE)
  fit <- fit_pk(obs, start = modifyList(pk_default(), list(Q = 1e-6)), pars = "ka")
  expect_true("C24" %in% fit$objective_stats)
  expect_equal(fit$cv, c(PK_CV_FLOOR, PK_TMAX_CV_FLOOR))
})

# ---- R1 fix round 3: C72 and parent half-lives ----
test_that("C72 equals the simulated concentration at 72 h after the last dose (analytic one-compartment case)", {
  p <- modifyList(pk_default(), list(ka = 1, fp = 0, fm = 0, CL = 50, V1 = 500, Q = 1e-9))
  k <- 50 / 500; dose <- mg_to_nmol(100, "NTX")
  bateman <- function(t) dose / 500 * 1 / (1 - k) * (exp(-k * t) - exp(-1 * t))
  obs <- data.frame(dose_mg = 100, analyte = "NTX", stat = "C72", regimen = c("single", "daily"), stringsAsFactors = FALSE)
  pr <- pk_predict_stats(p, obs)
  expect_equal(pr[1], nM_to_ngml(bateman(72), "NTX"), tolerance = 1e-4)
  expect_equal(pr[2], nM_to_ngml(sum(bateman(72 + 24 * (0:60))), "NTX"), tolerance = 1e-3)
})

test_that("pk_parent_half_lives gives the one- and two-compartment disposition half-lives", {
  p1 <- modifyList(pk_default(), list(CL = 50, V1 = 500, Q = 1e-9))
  expect_equal(pk_parent_half_lives(p1)[["terminal"]], log(2) / 0.1, tolerance = 1e-5)
  p2 <- modifyList(pk_default(), list(CL = 50, V1 = 500, Q = 100, V2 = 1000))
  K <- matrix(c(-(0.1 + 0.2), 0.2, 0.1, -0.1), 2)   # d(A1, A2)/dt
  ev <- sort(-eigen(K)$values)
  hl <- pk_parent_half_lives(p2)
  expect_equal(hl[["terminal"]], log(2) / ev[1], tolerance = 1e-6)
  expect_equal(hl[["distribution"]], log(2) / ev[2], tolerance = 1e-6)
})

test_that("a failed or non-finite simulation gives NA statistics rather than an error (so the optimiser can reject it)", {
  expect_true(all(is.na(refine_peak(c(0, 1, 2, 3), c(0, NA, NaN, 1)))))
  expect_true(all(is.na(nca_summary(c(0, 1, 2, 3), c(0, NA, 1, 0.5)))))
})


# ---- R4 / R4 fix round 1: dose-dependent first pass ----
# R4 fix round 1 rewrote this block: the dose dependence is now ANCHORED at the reference dose (100 mg, Verebey), so
# fp(100 mg) = fp exactly for every D50 (review I1), the low-dose parameter is fp0 (fp as dose -> 0; review I2: the old
# name is gone), a third structure (parent-only presystemic loss, floss0) was added (review I6), and the exposure-cache test
# was rewritten so that it fails when the dose dependence is switched off (review I10).
pk_dd <- function(fp = 0.80, fp0 = 0.98, D50 = 30) modifyList(pk_default(), list(fp = fp, fp0 = fp0, fp_D50 = D50))
pk_po <- function(fp = 0.80, floss0 = 0.8, D50 = 30) modifyList(pk_default(), list(fp = fp, floss0 = floss0, fp_D50 = D50))

test_that("fp_at_dose reduces exactly to the linear model when the dose dependence is switched off", {
  for (p in list(pk_dd(D50 = Inf), pk_dd(fp0 = 0.80), modifyList(pk_dd(), list(fp0 = NULL)), pk_default()))
    for (d in c(0.01, 4.5, 50, 100, 300)) expect_identical(fp_at_dose(d, p), p$fp)
  for (p in list(pk_po(D50 = Inf), modifyList(pk_po(), list(floss0 = NULL)), pk_default()))
    for (d in c(0.01, 4.5, 100)) { expect_identical(ploss_at_dose(d, p), 0); expect_identical(fp_at_dose(d, p), p$fp) }
})

test_that("the dose dependence is anchored at 100 mg: fp(100 mg) = fp and no parent is lost at 100 mg, for every D50", {
  for (D50 in c(10, 30, 100)) for (fp0 in c(0.5, 0.95, 0.99, 0.9999)) {
    p <- pk_dd(fp = 0.95, fp0 = fp0, D50 = D50)
    expect_equal(1 - fp_at_dose(100, p), 0.05, tolerance = 1e-12)
    expect_equal(pk_bioavailability(100, p), 0.05, tolerance = 1e-12)
  }
  for (D50 in c(10, 30, 100)) for (fl in c(0.1, 0.9)) {
    p <- pk_po(fp = 0.95, floss0 = fl, D50 = D50)
    expect_equal(ploss_at_dose(100, p), 0, tolerance = 1e-12)
    expect_equal(pk_bioavailability(100, p), 0.05, tolerance = 1e-12)
  }
  expect_equal(FP_DREF_MG, 100)
})

test_that("fp0 is the conversion fraction as dose -> 0, fp(D) is monotone between the anchor and fp0, and stays in (0, 1)", {
  p <- pk_dd(fp = 0.80, fp0 = 0.98)
  expect_equal(fp_at_dose(1e-9, p), 0.98, tolerance = 1e-6)
  v <- vapply(c(0.1, 1, 4.5, 10, 30, 100), fp_at_dose, 0, p = p)
  expect_true(all(diff(v) < 0))                                  # more conversion (less parent) at low dose
  q <- pk_dd(fp = 0.95, fp0 = 0.5, D50 = 100)                    # opposite direction: escape falls with dose above 100 mg
  w <- vapply(c(1, 100, 300, 1000, 1e5), fp_at_dose, 0, p = q)
  expect_true(all(w > 0 & w < 1))
})

test_that("parent-only loss: floss0 at dose -> 0, zero at and above the 100 mg anchor, and fp itself unchanged", {
  p <- pk_po(fp = 0.95, floss0 = 0.7, D50 = 30)
  expect_equal(ploss_at_dose(1e-9, p), 0.7, tolerance = 1e-6)
  expect_true(all(diff(vapply(c(0.1, 1, 4.5, 30, 100), ploss_at_dose, 0, p = p)) < 0))
  expect_equal(vapply(c(100, 200, 1000), ploss_at_dose, 0, p = p), c(0, 0, 0))
  expect_identical(fp_at_dose(4.5, p), 0.95)
  expect_equal(pk_bioavailability(4.5, p), 0.05 * (1 - ploss_at_dose(4.5, p)))
})

test_that("simulate_pk: fp(D) gives less parent AND more 6-beta-naltrexol at low dose; parent-only loss gives less parent and no more 6-beta-naltrexol; both equal the linear model at 100 mg", {
  lin <- modifyList(pk_default(), list(fp = 0.80)); t <- seq(0, 6, 0.25)
  dd <- pk_dd(fp = 0.80, fp0 = 0.98); po <- pk_po(fp = 0.80, floss0 = 0.8)
  s_lin <- simulate_pk(4.5, t, lin); s_dd <- simulate_pk(4.5, t, dd); s_po <- simulate_pk(4.5, t, po)
  expect_lt(max(s_dd$Cp_NTX), 0.5 * max(s_lin$Cp_NTX)); expect_gt(max(s_dd$Cp_BN), max(s_lin$Cp_BN))
  expect_lt(max(s_po$Cp_NTX), 0.5 * max(s_lin$Cp_NTX)); expect_lte(max(s_po$Cp_BN), max(s_lin$Cp_BN))
  for (q in list(dd, po)) {
    a <- simulate_pk(100, t, q); b <- simulate_pk(100, t, lin)
    expect_equal(a$Cp_NTX, b$Cp_NTX, tolerance = 1e-10); expect_equal(a$Cp_BN, b$Cp_BN, tolerance = 1e-10)
  }
  s0 <- simulate_pk(0, t, dd); expect_true(all(s0$Cp_NTX == 0) && all(s0$Cp_BN == 0))
})

test_that("superposition remains valid at a fixed dose under both dose-dependent structures", {
  for (p in list(pk_dd(), pk_po())) {
    p$Q <- 1e-6
    obs <- data.frame(dose_mg = 4.5, analyte = c("NTX", "BN"), stat = "Cmax", regimen = "daily", stringsAsFactors = FALSE)
    expect_equal(pk_predict_stats(p, obs, superposition = TRUE), pk_predict_stats(p, obs, superposition = FALSE), tolerance = 1e-4)
  }
})

test_that("the C0.25 / C0.5 / C0.75 / C1 stats equal the simulated concentration at 15, 30, 45 and 60 min (analytic one-compartment)", {
  p <- modifyList(pk_default(), list(ka = 1, fp = 0, fm = 0, CL = 50, V1 = 500, Q = 1e-9))
  k <- 50 / 500; dose <- mg_to_nmol(100, "NTX")
  bateman <- function(t) dose / 500 * 1 / (1 - k) * (exp(-k * t) - exp(-1 * t))
  obs <- data.frame(dose_mg = 100, analyte = "NTX", stat = c("C0.25", "C0.5", "C0.75", "C1"), regimen = "single", stringsAsFactors = FALSE)
  expect_equal(pk_predict_stats(p, obs), nM_to_ngml(bateman(c(0.25, 0.5, 0.75, 1)), "NTX"), tolerance = 1e-4)
  expect_true(all(c("C0.25", "C0.5", "C0.75", "C1", "Css13h") %in% PK_STATS))
})

test_that("Css13h stat equals the mean of the simulated curve over 1-3 h (analytic one-compartment, single dose)", {
  p <- modifyList(pk_default(), list(ka = 1, fp = 0, fm = 0, CL = 50, V1 = 500, Q = 1e-9))
  k <- 50 / 500; dose <- mg_to_nmol(100, "NTX")
  ant_b <- function(t) dose / 500 * 1 / (1 - k) * (-exp(-k * t) / k + exp(-1 * t) / 1)
  obs <- data.frame(dose_mg = 100, analyte = "NTX", stat = "Css13h", regimen = "single", stringsAsFactors = FALSE)
  expect_equal(pk_predict_stats(p, obs), nM_to_ngml((ant_b(3) - ant_b(1)) / 2, "NTX"), tolerance = 1e-3)
})

test_that("exposure_cache applies fp(D) per dose: the low-dose parent is far below the linear cache, the 100 mg profile is unchanged", {
  bc_fake <- structure(list(Kd_NTX = 0.47, Kd_BN = 2.12, koff_NTX = 21, koff_BN = 2, kon_NTX = 21 / 0.47, kon_BN = 2 / 2.12, keq = 2),
                       class = "binding_params")
  lin <- modifyList(pk_default(), list(fp = 0.80)); dd <- pk_dd(fp = 0.80, fp0 = 0.98)
  c_dd <- exposure_cache(c(1.5, 100), dd, bc_fake, n_doses = 14, max_doses = 14)
  c_lin <- exposure_cache(c(1.5, 100), lin, bc_fake, n_doses = 14, max_doses = 14)
  # fails if the cache ignores the dose dependence (then the two caches are identical)
  expect_lt(max(c_dd[["1.5"]]$Cu_NTX), 0.5 * max(c_lin[["1.5"]]$Cu_NTX))
  expect_equal(c_dd[["100"]]$Cu_NTX, c_lin[["100"]]$Cu_NTX, tolerance = 1e-8)
})

test_that("predict_pet uses fp(D) for each row's own dose", {
  bp <- list(Kd_NTX = 0.47, Kd_BN = 2.12, koff_NTX = 21, koff_BN = 2, kon_NTX = 21 / 0.47, kon_BN = 2 / 2.12, keq = 2)
  lin <- modifyList(pk_default(), list(fp = 0.80)); dd <- pk_dd(fp = 0.80, fp0 = 0.98)
  rows <- data.frame(dose_mg = c(4.5, 100), regimen = "single", time_h = 4, reference = "absolute", n_doses = 1, stringsAsFactors = FALSE)
  o_lin <- predict_pet(bp, rows, lin); o_dd <- predict_pet(bp, rows, dd)
  expect_lt(o_dd[1], o_lin[1]); expect_equal(o_dd[2], o_lin[2], tolerance = 1e-8)
})

test_that("pk_to_theta / pk_from_theta round-trip fp0 and floss0 on the logit scale (review M12)", {
  p <- modifyList(pk_default(), list(fp = 0.95, fp0 = 0.993, floss0 = 0.4))
  pars <- c("ka", "fp", "fp0", "floss0", "CL")
  th <- pk_to_theta(p, pars)
  expect_equal(th[["fp0"]], logit(0.993)); expect_equal(th[["floss0"]], logit(0.4)); expect_equal(th[["ka"]], log(p$ka))
  back <- pk_from_theta(th, pk_default(), pars)
  for (nm in pars) expect_equal(back[[nm]], p[[nm]], tolerance = 1e-12)
  expect_true(all(c("fp", "fp0", "floss0") %in% PK_LOGIT_PARS)); expect_equal(length(PK_LOGIT_PARS), 3L)
  expect_false(any(grepl("_hi$", names(pk_default()))))
})

# ---- R4 fix round 1 (review C2): weights from the reported dispersion, an assumed CV where none is reported, row weights ----
test_that("fit_pk weights each row by obj_weight / CV^2, with CV = sd / mean, else cv_assumed, else the floor", {
  obs <- data.frame(dose_mg = 100, analyte = "NTX", stat = c("Cmax", "C24", "Tmax", "C0.5"), regimen = "single",
                    mean = c(40, 2, 1, 10), sd = c(20, NA, NA, 8), cv_assumed = c(NA, 0.8, NA, NA), obj_weight = c(1, 1, 1, 0.25),
                    stringsAsFactors = FALSE)
  expect_equal(pk_obs_cv(obs), c(0.5, 0.8, PK_TMAX_CV_FLOOR, 0.8))
  expect_equal(pk_obs_weight(obs), c(1, 1, 1, 0.25))
  fit <- fit_pk(obs, start = modifyList(pk_default(), list(Q = 1e-6)), pars = "ka")
  pred <- pk_predict_stats(fit$par, obs)
  expect_equal(fit$value, sum(c(1, 1, 1, 0.25) * ((log(pred) - log(obs$mean)) / c(0.5, 0.8, 0.25, 0.8))^2), tolerance = 1e-8)
  expect_equal(sum(fit$row_wss), fit$value, tolerance = 1e-8)
  o2 <- obs[, c("dose_mg", "analyte", "stat", "regimen", "mean", "sd")]            # legacy files: no new columns
  expect_equal(pk_obs_cv(o2), c(0.5, PK_CV_FLOOR, PK_TMAX_CV_FLOOR, 0.8)); expect_equal(pk_obs_weight(o2), rep(1, 4))
})

# ---- R4 fix round 1 (review I7): the adoption rule, applied mechanically ----
test_that("pk_adoption_rule adopts a one-extra-parameter structure only on dWSS > 3.84, a PD Hessian and no Verebey row worse by > 1.5-fold", {
  cand <- function(dw, pd, worse) data.frame(dose_model = c("linear", "fp_dose", "parent_loss"), delta_wss = c(0, dw),
                                             hessian_pd = c(TRUE, pd), max_verebey_worsening = c(1, worse))
  expect_equal(pk_adoption_rule(cand(c(10, 5), c(TRUE, TRUE), c(1.1, 1.2)))$adopted, "fp_dose")       # both qualify: larger dWSS
  expect_equal(pk_adoption_rule(cand(c(5, 10), c(TRUE, TRUE), c(1.1, 1.2)))$adopted, "parent_loss")
  expect_equal(pk_adoption_rule(cand(c(3.8, 3.0), c(TRUE, TRUE), c(1, 1)))$adopted, "linear")          # below 3.84
  expect_equal(pk_adoption_rule(cand(c(50, 2), c(FALSE, TRUE), c(1, 1)))$adopted, "linear")            # Hessian not PD
  expect_equal(pk_adoption_rule(cand(c(50, 40), c(TRUE, TRUE), c(1.6, 1.2)))$adopted, "parent_loss")   # Verebey row worse 1.6x
  r <- pk_adoption_rule(cand(c(50, 40), c(TRUE, TRUE), c(1.6, 1.51)))
  expect_equal(r$adopted, "linear"); expect_equal(r$qualifies, c(FALSE, FALSE, FALSE)); expect_match(r$reason, "linear")
  expect_equal(PK_ADOPT_DWSS, stats::qchisq(0.95, 1)); expect_equal(PK_ADOPT_VEREBEY_FOLD, 1.5)
})

test_that("pk_verebey_worsening is the largest ratio of fold errors (alternative over linear) over the reference rows", {
  obs <- c(10, 10, 10); lin <- c(10, 5, 20); alt <- c(15, 5, 10)
  expect_equal(pk_verebey_worsening(obs, alt, lin), 1.5)
  expect_equal(pk_verebey_worsening(obs, lin, lin), 1)
})

test_that("PK_DOSE_MODELS names exactly one extra parameter per alternative structure and none for linear", {
  expect_equal(names(PK_DOSE_MODELS), c("linear", "fp_dose", "parent_loss"))
  expect_null(PK_DOSE_MODELS$linear$par); expect_equal(PK_DOSE_MODELS$fp_dose$par, "fp0"); expect_equal(PK_DOSE_MODELS$parent_loss$par, "floss0")
})

# ---- R4 fix round 1 (review C3): every PK output comes from the adopted fit ----
test_that("pk_write_outputs writes every PK output from one fit and pk_output_check catches a stale output", {
  obs <- data.frame(study = "S", dose_mg = c(4.5, 100), regimen = "single", analyte = "NTX", stat = "Cmax", unit = "ng/mL",
                    mean = c(1, 40), sd = c(0.5, 20), stringsAsFactors = FALSE)
  st <- modifyList(pk_default(), list(fp = 0.95, fp0 = 0.97, fp_D50 = 30, Q = 1e-6))
  fit <- fit_pk(obs, st, pars = c("CL", "fp0"))
  fit$structure <- "one-compartment"; fit$dose_model <- "fp_dose"; fit$fixed <- c("fp", "fp_D50")
  fit$hessian_summary <- pk_hessian_summary(fit)
  d <- tempfile("pkout"); dir.create(d)
  pk_write_outputs(fit, obs, d, notes = "test note")
  write.csv(data.frame(structure = "one-compartment", dose_model = "fp_dose", selected = TRUE), file.path(d, "pk_structure_summary.csv"), row.names = FALSE)
  write.csv(data.frame(study = "S", pred_selected = pk_predict_stats(fit$par, obs)), file.path(d, "pk_structure_comparison.csv"), row.names = FALSE)
  expect_equal(pk_output_check(d, fit, obs), character(0))
  tab <- read.csv(file.path(d, "pk_fit_table.csv")); expect_true(all(c("CL", "fp0", "fp", "fp_D50") %in% tab$parameter))
  expect_match(paste(readLines(file.path(d, "pk_fit_notes.txt")), collapse = " "), "fp0")
  # a stale table (the linear fit's, without fp0) is caught
  write.csv(tab[tab$parameter != "fp0", ], file.path(d, "pk_fit_table.csv"), row.names = FALSE)
  expect_true(any(grepl("pk_fit_table.csv", pk_output_check(d, fit, obs))))
  pk_write_outputs(fit, obs, d, notes = "test note")
  v <- read.csv(file.path(d, "pk_fit_vs_obs.csv")); v$pred <- v$pred * 1.1; write.csv(v, file.path(d, "pk_fit_vs_obs.csv"), row.names = FALSE)
  expect_true(any(grepl("pk_fit_vs_obs.csv", pk_output_check(d, fit, obs))))
  pk_write_outputs(fit, obs, d, notes = "test note")
  write.csv(data.frame(structure = "one-compartment", dose_model = "linear", selected = TRUE), file.path(d, "pk_structure_summary.csv"), row.names = FALSE)
  expect_true(any(grepl("pk_structure_summary.csv", pk_output_check(d, fit, obs))))
  unlink(d, recursive = TRUE)
})

# ---- R4 fix round 1 (review C1): sensitivity rows are selected by a stable key, never by label text ----
test_that("scripts 10 and 11 select PK sensitivity rows by the key column that script 03 writes, never by label text", {
  s01 <- paste(readLines(file.path("..", "..", "scripts", "01_fit_pk.R")), collapse = "\n")
  s03 <- paste(readLines(file.path("..", "..", "scripts", "03_test_pet.R")), collapse = "\n")
  s10 <- paste(readLines(file.path("..", "..", "scripts", "10_figures.R")), collapse = "\n")
  s11 <- paste(readLines(file.path("..", "..", "scripts", "11_tables.R")), collapse = "\n")
  for (s in list(s03, s10, s11)) {
    expect_false(grepl("grepl\\([^)]*\\$analysis", s))         # no label-text selection of occupancy-sensitivity rows
    expect_false(grepl("\\$analysis ==", s))
    expect_false(grepl("\\$analysis !=", s))
  }
  for (s in list(s10, s11)) expect_match(s, "osn$key", fixed = TRUE)
  expect_match(s03, "key = ", fixed = TRUE)
  expect_match(s01, "PK_SENS_KEYS", fixed = TRUE)
  # the keys are shared constants, so a label can be reworded without breaking a consumer
  expect_equal(unname(PK_SENS_KEYS[c("fp090", "fp099", "terminal")]), c("fp090", "fp099", "terminal"))
  expect_equal(unname(PK_SENS_KEYS[c("D50_10", "D50_100")]), c("D50_10", "D50_100"))
  expect_equal(OCC_SENS_KEYS, c(primary = "primary", koff_eq = "koff_eq"))
})
