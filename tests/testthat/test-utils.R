test_that("unit conversions use correct molecular weights", {
  expect_equal(ngml_to_nM(1, "NTX"), 1000 / 341.4, tolerance = 1e-10)
  expect_equal(round(ngml_to_nM(1, "NTX"), 3), 2.929)
  expect_equal(nM_to_ngml(ngml_to_nM(7.3, "BN"), "BN"), 7.3)
  expect_equal(mg_to_nmol(4.5, "NTX"), 4.5e6 / 341.4)
})

test_that("logit and inv_logit are inverses", {
  expect_equal(inv_logit(logit(c(0.01, 0.5, 0.97))), c(0.01, 0.5, 0.97))
})

# ---- R2 item 4 (review I3): truncated, not censored, normal sampling ----
tm_ref <- function(mean, sd, lower, upper) {   # independent reference for the analytic moments
  a <- (lower - mean) / sd; b <- (upper - mean) / sd
  Z <- stats::pnorm(b) - stats::pnorm(a); da <- stats::dnorm(a); db <- stats::dnorm(b)
  m <- mean + sd * (da - db) / Z
  v <- sd^2 * (1 + (a * da - b * db) / Z - ((da - db) / Z)^2)
  c(mean = m, sd = sqrt(v))
}

test_that("rtnorm_trunc samples a properly normalised truncated normal", {
  set.seed(1); x <- rtnorm_trunc(200000, 5.08, 1.5, 4, 10)
  expect_true(all(x > 4 & x < 10))
  expect_equal(mean(x == 4), 0)                   # no point mass at the bound (the censoring defect)
  tm <- tm_ref(5.08, 1.5, 4, 10)
  expect_equal(mean(x), tm[["mean"]], tolerance = 0.01)
  expect_equal(stats::sd(x), tm[["sd"]], tolerance = 0.01)
})

test_that("rtnorm_trunc handles an unbounded side, a zero SD and a reproducible seed", {
  set.seed(2); y <- rtnorm_trunc(5000, 0, 1, -Inf, Inf)
  expect_equal(mean(y), 0, tolerance = 0.05); expect_equal(stats::sd(y), 1, tolerance = 0.05)
  expect_equal(rtnorm_trunc(4, 6, 0, 4, 10), rep(6, 4))
  expect_equal(rtnorm_trunc(3, 20, 0, 4, 10), rep(10, 3))        # a degenerate mean outside the interval
  set.seed(3); a <- rtnorm_trunc(10, 5, 1.5, 4, 10)
  set.seed(3); b <- rtnorm_trunc(10, 5, 1.5, 4, 10)
  expect_identical(a, b)
  expect_error(rtnorm_trunc(2, 5, 1, 10, 4), "lower")
})

# ---- R2 fix round 1, IMPORTANT 1: the REPORTED baseline moments describe the ENROLLED population ----
# so the parent normal must be chosen such that its truncation to the trial's bounds HAS those moments.
test_that("trunc_moments agrees with an independent reference and handles infinite bounds", {
  expect_equal(unname(trunc_moments(5.08, 1.5, 4, 10)), unname(tm_ref(5.08, 1.5, 4, 10)))
  expect_equal(unname(trunc_moments(7.4, 1.9, 4, 10)), unname(tm_ref(7.4, 1.9, 4, 10)))
  expect_equal(unname(trunc_moments(0, 1, -Inf, Inf)), c(0, 1))
  expect_equal(trunc_moments(0, 1, 0, Inf)[["mean"]], sqrt(2 / pi))          # half-normal mean
})

test_that("match_trunc_parent returns a parent whose truncation has the REPORTED moments", {
  for (tg in list(c(6.25, 1.45, 4, 10), c(7.565, 1.605, 4, 10), c(5.08, 1.5, 0, 10))) {
    p <- match_trunc_parent(tg[1], tg[2], tg[3], tg[4])
    expect_true(p$matched)
    tmm <- trunc_moments(p$mu, p$sigma, tg[3], tg[4])
    expect_equal(tmm[["mean"]], tg[1], tolerance = 1e-6)
    expect_equal(tmm[["sd"]], tg[2], tolerance = 1e-6)
    expect_gte(p$sigma, tg[2])                     # a truncated normal is narrower than its parent
  }
  y <- match_trunc_parent(5.08, 1.5, 0, 10)       # [0, 10] is ~3.4 SD away: almost nothing truncated
  expect_equal(y$mu, 5.08, tolerance = 0.01); expect_equal(y$sigma, 1.5, tolerance = 0.01)
})

test_that("match_trunc_parent falls back to the reported moments, flagged, when no parent exists", {
  # a mean only 1.08 above the floor cannot carry an SD of 1.5 under any truncated normal on [4, 10]
  p <- match_trunc_parent(5.08, 1.5, 4, 10)
  expect_false(p$matched)
  expect_equal(c(p$mu, p$sigma), c(5.08, 1.5))
  expect_match(p$note, "could not be matched")
  q <- match_trunc_parent(6, 0, 4, 10)            # degenerate SD: no matching needed
  expect_true(q$matched); expect_equal(c(q$mu, q$sigma), c(6, 0))
  expect_error(match_trunc_parent(11, 1, 4, 10), "outside")
})
