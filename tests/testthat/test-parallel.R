root <- normalizePath(file.path("..", ".."), winslash = "/")

test_that("with_cluster runs a trivial job on 2 workers with the R/ sources loaded", {
  r <- with_cluster(2, function(cl) parallel::parLapply(cl, 1:4, function(i) c(i * 2, exists("simulate_pk", envir = globalenv(), inherits = FALSE))),
                    seed = 1, root = root)
  expect_equal(vapply(r, `[`, 0, 1), c(2, 4, 6, 8))
  expect_true(all(vapply(r, `[`, 0, 2) == 1))
})

test_that("with_cluster is reproducible under a fixed seed and stops its workers", {
  job <- function(cl) unlist(parallel::parLapply(cl, 1:6, function(i) stats::runif(1)))
  a <- with_cluster(2, job, seed = 2026, root = root)
  b <- with_cluster(2, job, seed = 2026, root = root)
  c <- with_cluster(2, job, seed = 7, root = root)
  expect_identical(a, b); expect_false(identical(a, c))
  cl_seen <- NULL
  with_cluster(2, function(cl) { cl_seen <<- cl; 1 }, seed = 1, root = root)
  expect_error(parallel::clusterEvalQ(cl_seen, 1))          # cluster was stopped on exit
  expect_error(with_cluster(2, function(cl) stop("boom"), seed = 1, root = root), "boom")
})

test_that("default_n_workers is between 1 and 10", {
  n <- default_n_workers(); expect_gte(n, 1); expect_lte(n, 10)
})

# ---- R1 fix round 4 item 2: the RNG kind must not leak out of with_cluster (R2 parallelises MCMC) ----
# clusterSetRNGStream() sets RNGkind("L'Ecuyer-CMRG") and restores only .Random.seed. When .Random.seed does not yet exist
# (a session that has not drawn a random number -- exactly script 03 before its first draw) it is REMOVED instead, so the
# L'Ecuyer kind is left in force and every later set.seed() draws from the wrong generator. Hence the fresh state here.
fresh_rng <- function() if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())

test_that("a seeded draw is identical before and after a with_cluster call, from a fresh RNG state", {
  mode <- c(a = 0, b = 1); cov <- diag(c(0.04, 0.09))
  fresh_rng(); before <- laplace_draws(mode, cov, 5, seed = 7)
  fresh_rng(); invisible(with_cluster(2, function(cl) 1L, seed = 123, root = root))
  expect_identical(RNGkind()[1], "Mersenne-Twister")
  expect_identical(laplace_draws(mode, cov, 5, seed = 7), before)
})

test_that("with_cluster restores the RNG kind and seed when both were already set", {
  set.seed(11, kind = "L'Ecuyer-CMRG"); before <- .Random.seed
  invisible(with_cluster(2, function(cl) 1L, seed = 5, root = root))
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
  expect_identical(.Random.seed, before)
  set.seed(1, kind = "Mersenne-Twister")
})

test_that("laplace_draws is reproducible whatever RNG kind is in force", {
  mode <- c(a = 0, b = 1); cov <- diag(c(0.04, 0.09))
  set.seed(3, kind = "Mersenne-Twister"); a <- laplace_draws(mode, cov, 4, seed = 7)
  set.seed(3, kind = "L'Ecuyer-CMRG"); b <- laplace_draws(mode, cov, 4, seed = 7)
  expect_identical(a, b)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")        # the caller's kind is left alone
  set.seed(1, kind = "Mersenne-Twister")
})

test_that("task_seeds gives per-task streams, so a parallel job is independent of the worker count", {
  s <- task_seeds(2026, 6)
  expect_length(s, 6)
  expect_true(all(vapply(s, function(x) x[1] %% 100 == 7, TRUE)))     # L'Ecuyer-CMRG streams
  expect_false(identical(s[[1]], s[[2]]))
  job <- function(cl) parallel::parLapply(cl, task_seeds(2026, 6), function(st) with_task_seed(st, stats::runif(2)))
  a <- with_cluster(2, job, root = root)
  b <- with_cluster(4, job, root = root)
  expect_identical(a, b)
  expect_false(identical(a[[1]], a[[2]]))
})

test_that("pk_select_structure gives the same answer with a sequential and a cluster-backed mapper", {
  obs <- data.frame(dose_mg = 100, analyte = "NTX", stat = c("Cmax", "Tmax", "thalf"), regimen = "single",
                    mean = c(40, 1, 10), sd = c(10, NA, 2), stringsAsFactors = FALSE)
  st <- modifyList(pk_default(), list(fp = 0.95))
  a <- pk_select_structure(obs, st, free = c("ka", "CL", "V1"))
  b <- with_cluster(2, function(cl) pk_select_structure(obs, st, free = c("ka", "CL", "V1"),
                                                       map_fits = function(X, FUN) parallel::parLapply(cl, X, FUN)),
                    root = root)
  expect_equal(a$fit$value, b$fit$value)
  expect_equal(unlist(a$fit$par[a$fit$pars]), unlist(b$fit$par[b$fit$pars]))
  expect_equal(a$one$value, b$one$value); expect_equal(a$two$value, b$two$value)
})
