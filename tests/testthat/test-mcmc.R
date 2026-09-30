test_that("adaptive Metropolis recovers a correlated Gaussian", {
  S <- matrix(c(1, 0.8, 0.8, 1), 2); Si <- solve(S)
  lp <- function(x) -0.5 * drop(t(x - c(1, -2)) %*% Si %*% (x - c(1, -2)))
  r <- am_sample(lp, c(a = 0, b = 0), n_iter = 30000, n_burn = 5000, seed = 1)
  expect_equal(unname(colMeans(r$draws)), c(1, -2), tolerance = 0.1)
  expect_equal(unname(apply(r$draws, 2, sd)), c(1, 1), tolerance = 0.1)
  expect_gt(r$accept, 0.1); expect_lt(r$accept, 0.6)
})

test_that("run_chains reports rhat and ess", {
  lp <- function(x) -0.5 * sum(x^2)
  r <- run_chains(lp, c(a = 0), n_chains = 4, n_iter = 6000, n_burn = 2000)
  expect_true(all(c("rhat", "ess_bulk") %in% names(r$summary)))
  expect_lt(max(r$summary$rhat), 1.05)
})

test_that("run_chains returns a chain id per draw row", {
  r <- run_chains(function(x) -0.5 * sum(x^2), c(a = 0), n_chains = 3, n_iter = 500, n_burn = 100)
  expect_equal(length(r$chain), nrow(r$draws)); expect_equal(sort(unique(r$chain)), 1:3)
})
