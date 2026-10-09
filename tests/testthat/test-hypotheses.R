priors <- data.frame(parameter = c("pbo_rel_mean", "pbo_logit_sd", "kpbo_per_day"),
                     value = c(0.12, 0.35, 0.1), sd = c(NA, NA, 0.5),
                     scale = c("relative_reduction", "logit", "per_day"))
flat <- list("0" = data.frame(hour = 0:24, Cu_NTX = 0, Cu_BN = 0, RO = 0),
             "4.5" = data.frame(hour = 0:24, Cu_NTX = 1, Cu_BN = 0, RO = 0.4))
arms <- data.frame(trial = "T", arm = c("pbo", "ldn"), schedule = c("0:0", "0:4.5"), day = 84,
                   outcome = "pct_reduction", baseline_nrs = 6, mean = c(12, 25), se = c(2, 2))

test_that("parameter vectors map to valid theta for every hypothesis", {
  for (h in c("H0", "H1", "H2", "H3")) {
    p <- init_par(h, "T"); th <- theta_from(p, h, list(IC50_TLR4 = 1000, Kd_OGFr = NA))
    expect_true(th$kout > 0 && th$kpbo > 0)
    if (h != "H0") expect_true(th$Emax > 0 && th$Emax < 1)
  }
})

test_that("log posterior is finite at the initial values", {
  for (h in c("H0", "H1", "H2", "H3")) {
    p <- init_par(h, "T")
    expect_true(is.finite(log_prior(p, h, priors, "T")))
    expect_true(is.finite(log_lik_arms(p, h, arms, flat, list(IC50_TLR4 = 1000, Kd_OGFr = NA))))
  }
})

test_that("pointwise log-likelihood has one entry per arm", {
  p <- init_par("H1", "T")
  expect_length(log_lik_arms(p, "H1", arms, flat, list(IC50_TLR4 = 1000, Kd_OGFr = NA), pointwise = TRUE), 2)
})

test_that("H3 with a fixed finite Kd uses it and needs no log10_Kd_OGFr", {
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = 50)
  p <- init_par("H3", "T"); p <- p[setdiff(names(p), "log10_Kd_OGFr")]
  expect_equal(theta_from(p, "H3", fx)$Kd_OGFr, 50)
  expect_true(is.finite(log_prior(p, "H3", priors, "T")))
  expect_true(is.finite(log_lik_arms(p, "H3", arms, flat, fx)))
  p2 <- init_par("H3", "T"); p2[["log10_Kd_OGFr"]] <- 2
  expect_equal(theta_from(p2, "H3", list(IC50_TLR4 = 1000, Kd_OGFr = NA))$Kd_OGFr, 100)
})

test_that("sigma_rel widens the likelihood: lowers a perfect fit, raises a bad fit", {
  p <- init_par("H0", "T")
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  pred <- arm_outcome(simulate_arm("H0", theta_from(p, "H0", fx), parse_schedule("0:0"), 84, flat, pmax_trial(p, "T")),
                      84, "pct_reduction", 6)
  perfect <- data.frame(trial = "T", arm = "pbo", schedule = "0:0", day = 84, outcome = "pct_reduction",
                        baseline_nrs = 6, mean = pred, se = 2)
  bad <- transform(perfect, mean = pred + 30)
  lo <- p; lo[["log_sigma_rel"]] <- log(0.001); hi <- p; hi[["log_sigma_rel"]] <- log(0.3)
  expect_lt(log_lik_arms(hi, "H0", perfect, flat, fx), log_lik_arms(lo, "H0", perfect, flat, fx))
  expect_gt(log_lik_arms(hi, "H0", bad, flat, fx), log_lik_arms(lo, "H0", bad, flat, fx))
})

test_that("outcome_scale returns 100 for pct_reduction and baseline for change_nrs", {
  expect_equal(outcome_scale("pct_reduction", 6), 100)
  expect_equal(outcome_scale("change_nrs", 6.5), 6.5)
  expect_equal(outcome_scale(c("pct_reduction", "change_nrs"), 7), c(100, 7))
})

test_that("log_prior errors clearly on missing priors", {
  expect_error(log_prior(init_par("H0", "T"), "H0", priors[-1, ], "T"), "pbo_rel_mean")
})

test_that("log_lik_arms errors on a non-finite se and is finite on the finite-se subset", {
  p <- init_par("H1", "T"); fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  bad <- arms; bad$se[2] <- NA
  expect_error(log_lik_arms(p, "H1", bad, flat, fx), "non-finite se")
  expect_true(is.finite(log_lik_arms(p, "H1", bad[is.finite(bad$se), ], flat, fx)))
})

# ---- R2 item 1 (review C2): the kpbo prior SE is on the per-day scale and must be converted for log_kpbo ----
byon <- data.frame(parameter = c("pbo_rel_mean", "pbo_logit_sd", "kpbo_per_day"),
                   value = c(0.1114, 0.1206, 0.106), sd = c(NA, NA, 0.00833),
                   scale = c("relative_reduction", "logit", "per_day"))

test_that("a per_day prior row is converted to the log scale by sd/value", {
  k <- log_scale_prior(byon, "kpbo_per_day")
  expect_true(k$converted)
  expect_equal(k$mean, log(0.106))
  expect_equal(k$sd, 0.00833 / 0.106)
  expect_equal(round(k$sd, 4), 0.0786)
})

test_that("a prior row already on the log scale is used as it stands", {
  lg <- byon; lg$scale[3] <- "log_per_day"; lg$value[3] <- log(0.106); lg$sd[3] <- 0.0786
  k <- log_scale_prior(lg, "kpbo_per_day")
  expect_false(k$converted)
  expect_equal(k$mean, log(0.106)); expect_equal(k$sd, 0.0786)
})

test_that("log_scale_prior refuses a row whose scale is absent, blank or unknown", {
  expect_error(log_scale_prior(byon[, setdiff(names(byon), "scale")], "kpbo_per_day"), "scale")
  blank <- byon; blank$scale[3] <- ""
  expect_error(log_scale_prior(blank, "kpbo_per_day"), "scale")
  odd <- byon; odd$scale[3] <- "per_week"
  expect_error(log_scale_prior(odd, "kpbo_per_day"), "per_week")
  dup <- rbind(byon, byon[3, ])
  expect_error(log_scale_prior(dup, "kpbo_per_day"), "exactly one")
})

test_that("log_prior applies the CONVERTED log-scale SD to log_kpbo, not the per-day SE", {
  p <- init_par("H0", "T")
  d <- function(x) { q <- p; q[["log_kpbo"]] <- x; log_prior(q, "H0", byon, "T") }
  # the only term that moves is the log_kpbo prior, so the difference is a difference of dnorm log-densities
  obs <- d(log(0.106) + 0.05) - d(log(0.106))
  expect_equal(obs, stats::dnorm(0.05, 0, 0.00833 / 0.106, log = TRUE) -
                      stats::dnorm(0, 0, 0.00833 / 0.106, log = TRUE))
  # with the uncorrected SE the same displacement would be ~9.4x further into the tail
  expect_gt(obs, stats::dnorm(0.05, 0, 0.00833, log = TRUE) - stats::dnorm(0, 0, 0.00833, log = TRUE))
})

test_that("the posterior for log_kpbo is wider than the uncorrected per-day SE of 0.00833", {
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  lp <- function(q) { pr <- log_prior(q, "H0", byon, "T"); if (!is.finite(pr)) return(-Inf)
                      pr + log_lik_arms(q, "H0", arms, flat, fx) }
  r <- am_sample(lp, init_par("H0", "T"), n_iter = 8000, n_burn = 3000, seed = 4)
  expect_gt(stats::sd(r$draws[, "log_kpbo"]), 0.02)
  expect_equal(stats::sd(r$draws[, "log_kpbo"]), 0.00833 / 0.106, tolerance = 0.5)
})

# ---- R2 item 3 (review I7): the model-structure residual is SHARED by the arms of one trial ----
test_that("outcome_shift carries the sign of d(outcome)/d(rel) and the magnitude of outcome_scale", {
  expect_equal(outcome_shift("pct_reduction", 6), -100)
  expect_equal(outcome_shift("change_nrs", 6.5), 6.5)
  expect_equal(abs(outcome_shift(c("pct_reduction", "change_nrs"), 7)), outcome_scale(c("pct_reduction", "change_nrs"), 7))
})

test_that("shared_resid_loglik equals a Monte Carlo integration over the shared term", {
  r <- c(1.2, -0.7, 0.4); se <- c(0.9, 1.3, 1.1); u <- c(0.8, 0.5, -0.6)
  set.seed(42); z <- stats::rnorm(400000)
  lw <- rowSums(vapply(seq_along(r), function(i) stats::dnorm(r[i], z * u[i], se[i], log = TRUE), numeric(length(z))))
  m <- max(lw); mc <- m + log(mean(exp(lw - m)))
  expect_equal(shared_resid_loglik(r, se, u), mc, tolerance = 1e-3)
})

test_that("shared_resid_loglik reduces to the independent formula for one arm and for sigma = 0", {
  expect_equal(shared_resid_loglik(1.1, 0.7, 0.4), stats::dnorm(1.1, 0, sqrt(0.7^2 + 0.4^2), log = TRUE))
  r <- c(1.2, -0.7); se <- c(0.9, 1.3)
  expect_equal(shared_resid_loglik(r, se, c(0, 0)), sum(stats::dnorm(r, 0, se, log = TRUE)))
})

# R2 fix round 1, IMPORTANT 2: PSIS-LOO for a non-factorising likelihood needs p(y_i | y_-i, theta) (Burkner,
# Gabry & Vehtari 2021), not chain-rule terms, which are a different quantity and depend on the row order.
dense_mvn <- function(r, S) { L <- t(chol(S)); w <- forwardsolve(L, r)
  -0.5 * length(r) * log(2 * pi) - sum(log(diag(L))) - 0.5 * sum(w^2) }
loo_ref <- function(r, se, u) {                     # p(y_i | y_-i) by conditioning a dense MVN, one i at a time
  S <- diag(se^2, length(r)) + tcrossprod(u)
  vapply(seq_along(r), function(i) {
    if (length(r) == 1) return(stats::dnorm(r, 0, sqrt(S[1, 1]), log = TRUE))
    m <- drop(S[i, -i] %*% solve(S[-i, -i], r[-i])); v <- S[i, i] - drop(S[i, -i] %*% solve(S[-i, -i], S[-i, i]))
    stats::dnorm(r[i], m, sqrt(v), log = TRUE) }, 0)
}

test_that("the pointwise terms are the leave-one-out conditionals p(y_i | y_-i)", {
  set.seed(5)
  for (k in 1:4) {
    r <- stats::rnorm(k, 0, 3); se <- stats::runif(k, 0.3, 4); u <- stats::rnorm(k, 0, 3)
    expect_equal(shared_resid_loglik_pointwise(r, se, u), loo_ref(r, se, u), tolerance = 1e-10)
  }
  expect_equal(shared_resid_loglik_pointwise(1.1, 0.7, 0.4), stats::dnorm(1.1, 0, sqrt(0.7^2 + 0.4^2), log = TRUE))
})

test_that("the pointwise terms do not depend on the order of the rows", {
  r <- c(1.2, -0.7, 0.4); se <- c(0.9, 1.3, 1.1); u <- c(0.8, 0.5, -0.6)
  pw <- shared_resid_loglik_pointwise(r, se, u)
  for (o in list(c(3, 2, 1), c(2, 3, 1), c(1, 3, 2)))
    expect_equal(shared_resid_loglik_pointwise(r[o], se[o], u[o]), pw[o], tolerance = 1e-12)
})

test_that("shared_resid_loglik equals a dense multivariate normal density", {
  set.seed(6)
  for (t in 1:50) {
    k <- sample(1:5, 1); r <- stats::rnorm(k, 0, 3); se <- stats::runif(k, 0.1, 5); u <- stats::rnorm(k, 0, 4)
    expect_equal(shared_resid_loglik(r, se, u), dense_mvn(r, diag(se^2, k) + tcrossprod(u)), tolerance = 1e-9)
  }
})

test_that("log_lik_arms equals a dense block-diagonal multivariate normal over trials", {
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  three <- rbind(transform(arms, trial = "T"), transform(arms[1, ], trial = "U", mean = 15, outcome = "change_nrs",
                                                         baseline_nrs = 6.5))
  p <- init_par("H1", c("T", "U")); p[["log_sigma_rel"]] <- log(0.15)
  th <- theta_from(p, "H1", fx); sig <- exp(p[["log_sigma_rel"]])
  m <- vapply(seq_len(nrow(three)), function(i) { a <- three[i, ]
    arm_outcome(simulate_arm("H1", th, parse_schedule(a$schedule), a$day, flat, pmax_trial(p, a$trial)),
                a$day, a$outcome, a$baseline_nrs) }, 0)
  u <- sig * outcome_shift(three$outcome, three$baseline_nrs)
  S <- diag(three$se^2)
  for (tr in unique(three$trial)) { ix <- which(three$trial == tr); S[ix, ix] <- S[ix, ix] + tcrossprod(u[ix]) }
  expect_equal(log_lik_arms(p, "H1", three, flat, fx), dense_mvn(three$mean - m, S), tolerance = 1e-9)
})

test_that("a shared term makes two trials independent but couples the arms within a trial", {
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  two <- rbind(transform(arms, trial = "T"), transform(arms, trial = "U"))
  p <- init_par("H0", c("T", "U")); p[["log_sigma_rel"]] <- log(0.1)
  expect_equal(log_lik_arms(p, "H0", two, flat, fx),
               log_lik_arms(p, "H0", two[two$trial == "T", ], flat, fx) +
                 log_lik_arms(p, "H0", two[two$trial == "U", ], flat, fx))
  # within one trial a COMMON displacement of both arms is cheaper than an OPPOSITE one of the same size;
  # under the old per-arm independent residual the two would score identically
  q <- init_par("H0", "T"); q[["log_sigma_rel"]] <- log(0.2)
  mfit <- vapply(1:2, function(i) arm_outcome(simulate_arm("H0", theta_from(q, "H0", fx), parse_schedule(arms$schedule[i]),
                   arms$day[i], flat, pmax_trial(q, "T")), arms$day[i], arms$outcome[i], arms$baseline_nrs[i]), 0)
  same <- transform(arms, mean = mfit - 10)       # both arms shifted the same way
  oppo <- transform(arms, mean = mfit + c(-10, 10))
  expect_gt(log_lik_arms(q, "H0", same, flat, fx), log_lik_arms(q, "H0", oppo, flat, fx))
})

test_that("log_lik_arms pointwise gives finite LOO conditionals that are invariant to row order", {
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  p <- init_par("H1", "T"); p[["log_sigma_rel"]] <- log(0.1)
  pw <- log_lik_arms(p, "H1", arms, flat, fx, pointwise = TRUE)
  expect_length(pw, 2); expect_true(all(is.finite(pw)))
  expect_equal(log_lik_arms(p, "H1", arms[2:1, ], flat, fx, pointwise = TRUE), rev(pw), tolerance = 1e-12)
  # the conditionals are NOT a decomposition of the joint (only the joint enters the posterior)
  expect_false(isTRUE(all.equal(sum(pw), log_lik_arms(p, "H1", arms, flat, fx))))
})

test_that("training_loo returns p_loo, the prior-determined fraction and a reliability flag", {
  fx <- list(IC50_TLR4 = 1000, Kd_OGFr = NA)
  lp <- function(q) { pr <- log_prior(q, "H0", byon, "T"); if (!is.finite(pr)) return(-Inf)
                      pr + log_lik_arms(q, "H0", arms, flat, fx) }
  f <- run_chains(lp, init_par("H0", "T"), n_chains = 2, n_iter = 1500, n_burn = 500)
  r <- suppressWarnings(training_loo(f, "H0", arms, flat, fx, n_keep = 400))
  expect_equal(r$n_params, ncol(f$draws))
  expect_equal(r$prior_determined_fraction, 1 - r$p_loo / r$n_params)
  expect_identical(r$p_loo_reliable, r$max_pareto_k <= PARETO_K_RELIABLE)
  expect_identical(r$ranking_valid, r$n_k_gt_0.7 == 0)
})

# ---- R2 fix round 2, ruling 3: the p_loo table says when the prior-determined fraction is indicative only ----
test_that("p_loo_table carries a per-row note tied to the Pareto k diagnostic", {
  lt <- data.frame(hypothesis = c("H0", "H1"), elpd_loo = c(-8, -7), se = c(1.8, 0.7), p_loo = c(1.4, 0.8),
                   n_params = c(6, 9), prior_determined_fraction = 1 - c(1.4, 0.8) / c(6, 9),
                   max_pareto_k = c(0.80, 0.68), n_k_gt_0.7 = c(1, 0), ranking_valid = c(FALSE, TRUE),
                   p_loo_reliable = c(FALSE, TRUE))
  t <- p_loo_table(lt, n_obs = 2)
  expect_true("note" %in% names(t))
  expect_match(t$note[1], "indicative only")
  expect_match(t$note[1], "Pareto k")
  expect_match(t$note[1], "0.80")
  expect_false(grepl("indicative only", t$note[2]))
  expect_match(t$note[2], "Pareto k")
  expect_equal(t$n_free_parameters, c(6, 9)); expect_equal(t$n_training_observations, c(2, 2))
})

# ---- R3a item 6, M2: both exclusion reasons for Younger 2009 ----
test_that("training_arms_table records the endpoint-mismatch reason from the data's own location column", {
  arms <- read_dataset(file.path("..", "..", "data", "training", "efficacy_training.csv"))
  u <- training_arms_table(arms)
  expect_equal(nrow(u), nrow(arms))
  y9 <- u[u$trial == "younger2009", ]
  expect_true(all(y9$status == "excluded"))
  expect_true(all(grepl("no dispersion reported in source", y9$reason)))
  expect_true(all(grepl("not a pain-specific score", y9$reason)))
  expect_true(all(grepl("overall fibromyalgia symptom severity", y9$reason)))
  y13 <- u[u$trial == "younger2013", ]
  expect_true(all(y13$status == "included" & y13$reason == ""))
})

# ---- R3a item 6, M7: Bruun-Plesner check with a model ED95 and a grid extended below 0.75 mg ----
test_that("dose_for_fraction interpolates on the rising limb and returns NA when the target is not crossed", {
  doses <- c(0, 1, 2, 4, 8); eff <- c(0, 0.5, 0.8, 0.9, 1)
  expect_equal(dose_for_fraction(doses, eff, 0.5), 1)
  expect_equal(dose_for_fraction(doses, eff, 0.95), 4 + (0.95 - 0.9) / (1 - 0.9) * 4)
  hump <- c(0, 1, 0.8, 0.6, 0.4)                     # peak at the lowest non-zero dose
  expect_equal(dose_for_fraction(doses, hump, 0.5), 0.5)
  expect_true(is.na(dose_for_fraction(doses, rep(0, 5), 0.5)))
})

test_that("bruun_check reports model ED50 and ED95, the peak and whether it sits at a grid edge", {
  cache_b <- list()
  for (d in c(0, 0.1, 0.5, 1, 2, 4, 6)) cache_b[[as.character(d)]] <- data.frame(hour = 0:24, Cu_NTX = d, Cu_BN = 0,
                                                                                RO = d / (d + 1))
  fx <- list(IC50_TLR4 = 1, Kd_OGFr = 1)
  mk <- function(h, extra) list(hyp = h, fixed = fx,
    draws = cbind(mu_p = 0, log_tau = 0, log_kpbo = log(0.1), log_kout = log(0.5), log_sigma_rel = log(0.01), t(extra)))
  fits <- list(H2 = mk("H2", c(logit_Emax = logit(0.5))), H3 = mk("H3", c(logit_Emax = logit(0.5), log10_Kd_OGFr = 0)))
  b <- bruun_check(fits, cache_b, doses = c(0, 0.1, 0.5, 1, 2, 4, 6), obs_ED50 = 3.88, obs_ED95 = 5.4, day = 120)
  expect_true(all(c("grid", "summary") %in% names(b)))
  s <- b$summary
  expect_true(all(c("hypothesis", "peak_dose_mg", "peak_position", "model_ED50_mg", "model_ED95_mg",
                    "observed_ED50_mg", "observed_ED95_mg", "lowest_nonzero_dose_mg") %in% names(s)))
  expect_equal(s$observed_ED95_mg, c(5.4, 5.4)); expect_equal(s$lowest_nonzero_dose_mg, c(0.1, 0.1))
  # occupancy-driven effect rises monotonically on this grid: peak at the top, ED50 < ED95 <= top dose
  expect_true(all(s$peak_position == "highest dose on grid"))
  # R3a fix round 1: an effect still rising at the grid top has no estimable EDx (it would move with the grid top)
  expect_true(all(is.na(s$model_ED50_mg) & is.na(s$model_ED95_mg)))
  expect_true(all(grepl("not estimable: effect still rising at grid top", s$ed_note, fixed = TRUE)))
  expect_true(all(grepl("quantal", s$comparison_note) & grepl("descriptive", s$comparison_note)))
  expect_equal(nrow(b$grid), 2 * 7)
  expect_error(bruun_check(fits, cache_b, doses = c(0, 0.25), obs_ED50 = 3.88, obs_ED95 = 5.4), "0.25")
})

test_that("the exposure-cache dose grid of script 02 reaches below 0.75 mg for the Bruun check", {
  s02 <- readLines(file.path("..", "..", "scripts", "02_binding_bottomup.R"))
  ln <- s02[startsWith(s02, "doses <- c(")]
  expect_length(ln, 1)
  doses <- eval(parse(text = sub("^doses <- ", "", ln)))
  expect_true(all(c(0.1, 0.25, 0.5) %in% doses))
  s04 <- paste(readLines(file.path("..", "..", "scripts", "04_fit_hypotheses.R")), collapse = "\n")
  expect_true(grepl("bruun_check(", s04, fixed = TRUE))
})


# ---- R3a fix round 1, IMPORTANT 2: EDx only where the peak is interior ----
test_that("bruun_ed gives EDx for an interior peak and NA with a reason at either grid edge", {
  doses <- c(0, 0.1, 0.5, 1, 2, 4, 6)
  hump <- c(0, 0.3, 0.8, 1, 0.9, 0.7, 0.6)                     # peak interior at 1 mg
  r <- bruun_ed(doses, hump)
  expect_equal(r$peak_position, "interior")
  expect_equal(r$ED50, dose_for_fraction(doses, hump, 0.5)); expect_equal(r$ED95, dose_for_fraction(doses, hump, 0.95))
  expect_true(is.finite(r$ED50) && is.finite(r$ED95)); expect_equal(r$ed_note, "")
  rising <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)
  e <- bruun_ed(doses, rising)
  expect_equal(e$peak_position, "highest dose on grid")
  expect_true(is.na(e$ED50) && is.na(e$ED95))
  expect_match(e$ed_note, "not estimable: effect still rising at grid top", fixed = TRUE)
  low <- c(0, 1, 0.8, 0.6, 0.4, 0.3, 0.2)
  l <- bruun_ed(doses, low)
  expect_true(is.na(l$ED50) && is.na(l$ED95)); expect_match(l$ed_note, "not estimable")
  z <- bruun_ed(doses, rep(0, 7)); expect_true(is.na(z$ED95)); expect_match(z$ed_note, "no effect")
})
