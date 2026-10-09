# scripts/12_derived_quantities.R and its helpers (R/derived.R).

test_that("script 12 runs in the primary pass right after the tables, and only there", {
  expect_equal(which(PRIMARY_STEPS == "12_derived_quantities.R"), which(PRIMARY_STEPS == "11_tables.R") + 1L)
  expect_false("12_derived_quantities.R" %in% SENSITIVITY_STEPS)
  expect_true(file.exists(file.path("..", "..", "scripts", "12_derived_quantities.R")))
})

test_that("script 12 never references test data or the freeze manifest, and writes only derived_quantities.csv", {
  l <- readLines(file.path("..", "..", "scripts", "12_derived_quantities.R"))
  code <- sub("#.*$", "", l)
  expect_false(any(grepl("data/test", l, fixed = TRUE) | grepl("data\\test", l, fixed = TRUE) | grepl("read_test_data", l, fixed = TRUE)))
  expect_false(any(grepl("freeze_manifest|freeze_stage|manifest_path", code)))
  w <- code[grepl("write[.]csv|saveRDS|writeLines|cat\\(.*file *=|file[.]copy|sink\\(", code)]
  expect_length(w, 1)
  expect_match(w, "derived_quantities.csv", fixed = TRUE)
})

test_that("dq_row has the derived_quantities.csv columns", {
  r <- dq_row("x", 1.5, unit = "nM", definition = "d", source_files = c("a.rds", "b.csv"))
  expect_identical(names(r), DERIVED_COLS)
  expect_identical(DERIVED_COLS, c("quantity", "value", "lo", "hi", "unit", "definition", "source_files"))
  expect_true(is.na(r$lo) && is.na(r$hi))
  expect_equal(r$source_files, "a.rds; b.csv")
  expect_identical(names(rbind(r, dq_row("y", 2, 1, 3, "pp", "e", "c.rds"))), DERIVED_COLS)
})

test_that("min_ess_tail returns the smallest posterior::ess_tail and its parameter", {
  set.seed(1); n <- 300; nch <- 4
  slow <- unlist(lapply(seq_len(nch), function(k) as.numeric(stats::arima.sim(list(ar = 0.97), n))))
  fit <- list(draws = cbind(a = stats::rnorm(n * nch), slow = slow, b = stats::rnorm(n * nch)), chain = rep(seq_len(nch), each = n))
  e <- min_ess_tail(fit)
  expect_equal(e$parameter, "slow")
  arr <- array(slow, c(n, nch))
  expect_equal(e$value, posterior::ess_tail(arr))
  expect_error(min_ess_tail(list(draws = fit$draws[-1, ], chain = fit$chain[-1])), "unequal")
})

test_that("thinned_index and mean_unbound_ntx follow the documented conventions", {
  i <- thinned_index(120000, 400)
  expect_length(i, 400); expect_equal(range(i), c(1, 120000))
  expect_equal(thinned_index(10, 400), 1:10)
  prof <- data.frame(hour = c(0, 6, 12, 18, 24), Cu_NTX = c(1, 2, 3, 4, 1))
  expect_equal(mean_unbound_ntx(prof), 2.5)          # the 24 h endpoint (= hour 0 at steady state) is excluded
})

test_that("pk_laplace_draw_summary back-transforms script 03's PK draws", {
  f <- file.path("..", "..", "output", "results", "pk_fit.rds")
  skip_if_not(file.exists(f), "no pk_fit.rds")
  pk <- readRDS(f); skip_if_not(isTRUE(pk$hessian_summary$pd) && "fp0" %in% pk$pars, "PK Hessian not PD or no fp0")
  s <- pk_laplace_draw_summary(pk, 20, seed = 7, dose_mg = 4.5)
  d <- laplace_draws(pk$theta, pk$hessian_summary$cov, 20, seed = 7)
  expect_equal(s$fp0, inv_logit(unname(d[, "fp0"])))
  expect_equal(s$F[1], pk_bioavailability(4.5, pk_from_theta(d[1, ], pk$par, pk$pars)))
  expect_true(all(s$F > 0 & s$F < 1))
})
