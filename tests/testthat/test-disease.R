flat_cache <- function(ro = 0.4, cu = 1) list(
  "0" = data.frame(hour = seq(0, 24, 1), Cu_NTX = 0, Cu_BN = 0, RO = 0),
  "4.5" = data.frame(hour = seq(0, 24, 1), Cu_NTX = cu, Cu_BN = 0, RO = ro))
th <- list(kout = 0.1, kpbo = 0.1, Emax = 0.5, kA = 0.2, gamma = 1, IC50_TLR4 = 1, Kd_OGFr = 1)

test_that("H0 and placebo arms change only through the placebo term", {
  s <- simulate_arm("H0", th, parse_schedule("0:4.5"), 84, flat_cache(), Pmax = 0.12)
  expect_equal(s$S, rep(1, 85))
  expect_equal(s$rel[85], 1 - 0.12 * (1 - exp(-0.1 * 84)))
  p <- simulate_arm("H1", th, parse_schedule("0:0"), 84, flat_cache(), Pmax = 0.12)
  expect_equal(p$S, rep(1, 85))
})

test_that("H2 steady state equals 1 - Emax * mean occupancy of TLR4", {
  s <- simulate_arm("H2", th, parse_schedule("0:4.5"), 400, flat_cache(cu = 1), Pmax = 0)
  expect_equal(tail(s$S, 1), 1 - 0.5 * 0.5, tolerance = 1e-6)
})

test_that("H1 steady state equals 1 - Emax * mRO * mean((1-RO)^gamma)", {
  s <- simulate_arm("H1", th, parse_schedule("0:4.5"), 600, flat_cache(ro = 0.4), Pmax = 0)
  expect_equal(tail(s$S, 1), 1 - 0.5 * 0.4 * 0.6, tolerance = 1e-5)
})

test_that("dose lookup respects titration boundaries", {
  cache <- c(flat_cache(), list("1.5" = data.frame(hour = seq(0, 24, 1), Cu_NTX = 0, Cu_BN = 0, RO = 0)))
  s <- simulate_arm("H1", th, parse_schedule("0:1.5;7:4.5"), 14, cache, Pmax = 0)
  expect_equal(s$S[1:8], rep(1, 8))   # days 0-6 on the RO = 0 profile
  expect_lt(s$S[15], 1)               # effect after day 7
})

test_that("daily-step H1 matches a fine-grid integration within 1%", {
  hours <- seq(0, 24, 0.25); ro <- 0.5 + 0.4 * sin(2 * pi * hours / 24)
  cache <- list("4.5" = data.frame(hour = hours, Cu_NTX = 0, Cu_BN = 0, RO = ro))
  s <- simulate_arm("H1", th, parse_schedule("0:4.5"), 84, cache, Pmax = 0)
  dt <- 0.25 / 24; A <- 0; S <- 1; f <- stats::approxfun(hours, ro)
  for (i in seq_len(84 / dt)) {
    r <- f(((i - 1) * dt * 24) %% 24)
    A <- A + dt * th$kA * (r - A)
    S <- S + dt * th$kout * ((1 - th$Emax * A * (1 - r)^th$gamma) - S)
  }
  message("fine-grid S(84) = ", format(S, digits = 10), "; daily-step S(84) = ",
          format(tail(s$S, 1), digits = 10), "; rel diff = ", format(abs(tail(s$S, 1) - S) / S, digits = 4))
  expect_equal(tail(s$S, 1), S, tolerance = 0.01)
})

test_that("arm_outcome maps rel to the reported scale", {
  sim <- data.frame(day = 0:2, S = 1, Pbo = 0, rel = c(1, 0.9, 0.8))
  expect_equal(arm_outcome(sim, 2, "pct_reduction", 6), 20)
  expect_equal(arm_outcome(sim, 2, "change_nrs", 6), -1.2)
})

test_that("vectorised arm equals scalar arm for each subject", {
  v <- simulate_arm_vec("H1", th, parse_schedule("0:4.5"), 84, flat_cache(), Pmax_vec = c(0.1, 0.2), Emax_vec = c(0.3, 0.6))
  a <- simulate_arm("H1", modifyList(th, list(Emax = 0.6)), parse_schedule("0:4.5"), 84, flat_cache(), Pmax = 0.2)
  expect_equal(v[, 2], a$rel, tolerance = 1e-12)
})

test_that("NA in unused theta components does not warn or error", {
  th2 <- modifyList(th, list(Kd_OGFr = NA_real_, IC50_TLR4 = NA_real_))
  for (h in c("H0", "H1")) expect_no_warning(simulate_arm(h, th2, parse_schedule("0:4.5"), 10, flat_cache(), Pmax = 0.1))
  th3 <- modifyList(th, list(Kd_OGFr = NA_real_))
  expect_no_warning(simulate_arm("H2", th3, parse_schedule("0:4.5"), 10, flat_cache(), Pmax = 0.1))
})

test_that("H3 steady state equals 1 - Emax * mean(Cu/(Kd+Cu))", {
  s <- simulate_arm("H3", th, parse_schedule("0:4.5"), 400, flat_cache(cu = 3), Pmax = 0)
  expect_equal(tail(s$S, 1), 1 - 0.5 * 3 / (1 + 3), tolerance = 1e-6)
})

test_that("required parameters are enforced per hypothesis", {
  sch <- parse_schedule("0:4.5")
  expect_error(simulate_arm("H3", modifyList(th, list(Kd_OGFr = NA_real_)), sch, 10, flat_cache(), 0), "H3 requires finite theta")
  th_no <- th; th_no$IC50_TLR4 <- NULL
  expect_error(simulate_arm("H2", th_no, sch, 10, flat_cache(), 0), "H2 requires finite theta")
  th_ne <- th; th_ne$Emax <- NULL
  expect_error(simulate_arm("H1", th_ne, sch, 10, flat_cache(), 0), "Emax")
  expect_error(simulate_arm("H1", modifyList(th, list(Emax = 1.2)), sch, 10, flat_cache(), 0), "Emax")
  expect_error(simulate_arm("H1", modifyList(th, list(kout = NA_real_)), sch, 10, flat_cache(), 0), "kout")
  expect_error(simulate_arm("H9", th, sch, 10, flat_cache(), 0), "unknown hypothesis")
  expect_no_error(simulate_arm("H0", th_ne, sch, 10, flat_cache(), 0))
})
