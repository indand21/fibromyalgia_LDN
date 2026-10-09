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

# ---- R2 item 5: chains are seeded PER CHAIN, so a parallel fit reproduces the sequential one exactly ----
test_that("run_chain depends only on the seed and the chain index, not on the ambient RNG state or kind", {
  lp <- function(x) -0.5 * sum(x^2)
  a <- run_chain(lp, c(a = 0), n_iter = 400, n_burn = 100, seed = 2026, k = 2)
  set.seed(999); invisible(stats::runif(37))
  b <- run_chain(lp, c(a = 0), n_iter = 400, n_burn = 100, seed = 2026, k = 2)
  set.seed(5, kind = "L'Ecuyer-CMRG")
  c3 <- run_chain(lp, c(a = 0), n_iter = 400, n_burn = 100, seed = 2026, k = 2)
  set.seed(1, kind = "Mersenne-Twister")
  expect_identical(a$draws, b$draws); expect_identical(a$draws, c3$draws)
  d <- run_chain(lp, c(a = 0), n_iter = 400, n_burn = 100, seed = 2026, k = 3)
  expect_false(identical(a$draws, d$draws))
})

test_that("run_chains with a parallel mapper is identical to the sequential run", {
  lp <- function(x) -0.5 * sum(x^2)
  seq_fit <- run_chains(lp, c(a = 0, b = 0), n_chains = 4, n_iter = 600, n_burn = 200)
  root <- normalizePath(file.path("..", ".."), winslash = "/")
  par_fit <- with_cluster(2, function(cl) run_chains(lp, c(a = 0, b = 0), n_chains = 4, n_iter = 600, n_burn = 200,
                            map_chains = function(X, FUN) parallel::parLapply(cl, X, FUN)), root = root)
  expect_identical(seq_fit$draws, par_fit$draws)
  expect_identical(seq_fit$accept, par_fit$accept)
  expect_identical(seq_fit$chain, par_fit$chain)
})

test_that("combine_chains assembles a chain id per draw and a posterior summary", {
  lp <- function(x) -0.5 * sum(x^2)
  ch <- lapply(1:3, function(k) run_chain(lp, c(a = 0), n_iter = 400, n_burn = 100, seed = 7, k = k))
  f <- combine_chains(ch)
  expect_equal(nrow(f$draws), 3 * 300)
  expect_equal(sort(unique(f$chain)), 1:3)
  expect_true(all(c("rhat", "ess_bulk") %in% names(f$summary)))
})

# ---- R2 fix round 1, Minor 1: arguments given as VARIABLES must reach the workers ----
# A promise for an argument whose expression is a global variable points at the caller's global environment, which
# a PSOCK worker does not have; run_chains() must force its arguments before handing the closure to map_chains.
test_that("run_chains with a parallel mapper works when its arguments are global variables", {
  root <- normalizePath(file.path("..", ".."), winslash = "/")
  assign("r2fx_lp", function(x) -0.5 * sum(x^2), envir = globalenv())
  assign("r2fx_init", c(a = 0, b = 0), envir = globalenv())
  assign("r2fx_n", 300L, envir = globalenv())
  on.exit(rm(list = c("r2fx_lp", "r2fx_init", "r2fx_n"), envir = globalenv()), add = TRUE)
  seq_fit <- eval(quote(run_chains(r2fx_lp, r2fx_init, n_chains = 2, n_iter = r2fx_n, n_burn = 100)), globalenv())
  par_fit <- with_cluster(2, function(cl) {
    mapper <- function(X, FUN) parallel::parLapply(cl, X, FUN)
    eval(substitute(run_chains(r2fx_lp, r2fx_init, n_chains = 2, n_iter = r2fx_n, n_burn = 100, map_chains = M),
                    list(M = mapper)), globalenv())
  }, root = root)
  expect_identical(seq_fit$draws, par_fit$draws)
})

# ---- R2 fix round 1, Minor 2: the seed pins the normal and sample kinds too ----
test_that("run_chain is identical whatever normal.kind and sample.kind are in force", {
  lp <- function(x) -0.5 * sum(x^2)
  a <- run_chain(lp, c(a = 0), n_iter = 300, n_burn = 100, seed = 2026, k = 1)
  old <- RNGkind()
  RNGkind(normal.kind = "Box-Muller"); suppressWarnings(RNGkind(sample.kind = "Rounding"))
  b <- suppressWarnings(run_chain(lp, c(a = 0), n_iter = 300, n_burn = 100, seed = 2026, k = 1))  # restoring "Rounding" warns
  expect_identical(RNGkind()[2:3], c("Box-Muller", "Rounding"))     # the caller's kinds are restored
  suppressWarnings(do.call(RNGkind, as.list(old))); set.seed(1)
  expect_identical(a$draws, b$draws)
})
