# R3a items 2-4 and M3: magnitude-guarded Sobol and power, responder scoring, SD sensitivity, two-sided p-values.
root <- normalizePath(file.path("..", ".."), winslash = "/")
s_cache <- list("0" = data.frame(hour = 0:24, Cu_NTX = 0, Cu_BN = 0, RO = 0),
                "4.5" = data.frame(hour = 0:24, Cu_NTX = 1, Cu_BN = 0, RO = 0.4))
s_fit <- function(hyp, emax = 0.3, n = 3, sigma = 0.05) {
  d <- cbind(mu_p = seq(-1.6, -1.4, length.out = n), log_tau = log(0.2), log_kpbo = log(0.1),
             log_kout = seq(log(0.08), log(0.12), length.out = n), log_sigma_rel = log(sigma))
  if (hyp != "H0") d <- cbind(d, logit_Emax = seq(logit(emax) - 0.1, logit(emax) + 0.1, length.out = n))
  if (hyp == "H3") d <- cbind(d, log10_Kd_OGFr = seq(-0.1, 0.1, length.out = n))
  list(hyp = hyp, fixed = list(IC50_TLR4 = 1e9, Kd_OGFr = NA), draws = d)
}

# ---- M3: one-sided and two-sided posterior predictive p-values ----
test_that("ppp_values gives both tails, the old direction-of-observation tail and twice the smaller tail", {
  sim <- c(-2, -1, 0, 1, 2, 3)
  p <- ppp_values(sim, 1)
  expect_equal(p$p_upper, mean(sim >= 1)); expect_equal(p$p_lower, mean(sim <= 1))
  expect_equal(p$p_one_sided_toward_sign_of_obs, mean(sim >= 1))     # obs >= 0: upper tail (old p_tail)
  expect_equal(ppp_values(sim, -1)$p_one_sided_toward_sign_of_obs, mean(sim <= -1))
  expect_equal(p$p_two_sided, 2 * min(mean(sim >= 1), mean(sim <= 1)))
  expect_equal(ppp_values(sim, 0.5)$p_two_sided, 1)                  # capped at 1
  expect_true(is.na(ppp_values(sim, NA)$p_two_sided))
})

# ---- item 3 (review I6): responder scoring ----
test_that("observed_responder_counts parses the curated (k of n) counts and cross-checks them", {
  tt <- read_dataset(file.path("..", "..", "data", "test", "trials_test.csv"))
  oc <- observed_responder_counts(tt, "final", 84)
  expect_equal(oc$cut, c(15, 30, 50))
  r30 <- oc[oc$cut == 30, ]
  expect_equal(c(r30$k_ldn, r30$n_ldn, r30$k_pbo, r30$n_pbo), c(20, 49, 13, 50))
  oi <- observed_responder_counts(tt, "innova", 90)
  expect_equal(oi$k_ldn, c(8, 5, 1)); expect_equal(oi$n_pbo, c(49, 49, 49))
  bad <- tt; i <- which(bad$trial == "final" & bad$endpoint == "resp30" & bad$arm == "placebo"); bad$mean[i] <- 0.40
  expect_error(observed_responder_counts(bad, "final", 84), "disagree")
})

test_that("responder_scores: binomial predictive probabilities, two-sided and one-sided p and the Fisher test", {
  S <- 4000; set.seed(4)
  rl <- stats::rbeta(S, 30, 70); rp <- stats::rbeta(S, 30, 70)        # no drug effect
  s <- responder_scores(rl, rp, k_ldn = 20, n_ldn = 49, k_pbo = 13, n_pbo = 50, seed = 1, n_binom = 5)
  expect_equal(s$pp_prob_count_ldn, mean(stats::dbinom(20, 49, rl)))
  expect_equal(s$pp_prob_counts_joint, mean(stats::dbinom(20, 49, rl) * stats::dbinom(13, 50, rp)))
  expect_equal(s$obs_diff, 20 / 49 - 13 / 50)
  expect_true(s$ppp_one_sided_upper > 0 && s$ppp_one_sided_upper < 0.5)
  expect_equal(s$ppp_two_sided, min(1, 2 * min(s$ppp_one_sided_upper, s$ppp_one_sided_lower)))
  expect_equal(s$fisher_p_two_sided, stats::fisher.test(matrix(c(20, 29, 13, 37), 2, byrow = TRUE))$p.value)
  expect_equal(s$fisher_label, "source-trial test, descriptive")
  expect_equal(s$n_pred_counts, S * 5)
  # with a large true difference the observed +15 pp is unremarkable in the upper tail
  big <- responder_scores(rl + 0.15, rp, 20, 49, 13, 50, seed = 1)
  expect_gt(big$ppp_one_sided_upper, s$ppp_one_sided_upper)
  # reproducible and RNG-neutral
  set.seed(9); a <- stats::runif(1); set.seed(9)
  s2 <- responder_scores(rl, rp, 20, 49, 13, 50, seed = 1, n_binom = 5)
  expect_identical(stats::runif(1), a); expect_identical(s, s2)
  expect_error(responder_scores(rl, rp[-1], 20, 49, 13, 50, seed = 1), "same length")
})

# ---- item 4 (review I10): SD sensitivity of the variability calibration ----
test_that("override_sd_pct replaces only the Younger 2013 sd_pct values, in memory", {
  arms <- read_dataset(file.path("..", "..", "data", "training", "efficacy_training.csv"))
  a2 <- override_sd_pct(arms, pbo = 14.6, ldn = 12.5)
  y <- grepl("younger2013", a2$trial)
  expect_equal(a2$sd_pct[y & a2$arm == "placebo"], 14.6); expect_equal(a2$sd_pct[y & a2$arm != "placebo"], 12.5)
  expect_equal(a2[!y, ], arms[!y, ])
  expect_equal(arms$sd_pct[y & arms$arm == "placebo"], 27.9)          # the curated table is untouched
  expect_error(override_sd_pct(arms[!y, ], 14.6, 12.5), "Younger")
})

test_that("calibrate_omegas reports the noise-term SD of its own simulated placebo cohort", {
  dr <- cbind(mu_p = logit(0.1114), log_kpbo = log(0.106), log_kout = log(0.1)); dr <- cbind(dr, `eta_younger2013` = 0)
  fit <- list(hyp = "H0", fixed = list(), draws = dr)
  arms <- data.frame(trial = "younger2013", arm = "placebo", schedule = "0:0", day = 28, baseline_nrs = 50.8,
                     scale = "0-100", sd_pct = 27.9)
  bt <- data.frame(trial = "younger2013", lower_nrs = 0.5, upper_nrs = 10)
  om <- suppressWarnings(calibrate_omegas(fit, arms, s_cache, sd_res = 0.48, n = 4000, bounds = bt, seed = 1001))
  # the same cohort by hand: H0 placebo draws pmax (n), then the truncated baseline, then the noise
  p <- match_trunc_parent(5.08, 1.5, 0.5, 10)
  set.seed(1001); invisible(stats::rnorm(4000)); x <- rtnorm_trunc(4000, p$mu, p$sigma, 0.5, 10); e <- stats::rnorm(4000, 0, 0.48)
  expect_equal(om$noise_sd_pbo, stats::sd(100 * e / x))
  expect_equal(om$min_baseline_pbo, min(x))
})

# ---- item 2 (review I2): magnitude guard, not a name guard ----
test_that("net_effect_at is zero for H0 and positive for a drug hypothesis; days come from one simulation", {
  expect_equal(net_effect_at(s_fit("H0"), s_cache, 84), 0)
  e <- net_effect_at(s_fit("H3"), s_cache, c(14, 84))
  expect_length(e, 2); expect_true(all(e > 0)); expect_lt(e[1], e[2])
  expect_equal(e[2], net_effect_at(s_fit("H3"), s_cache, 84))
})

test_that("effect_eligibility applies the threshold to the day-84 net effect and labels each hypothesis", {
  fits <- list(H0 = s_fit("H0"), H2 = s_fit("H3", emax = 1e-5), H3 = s_fit("H3"))
  fits$H2$hyp <- "H3"                                         # an effect of ~0 under a drug structure
  el <- effect_eligibility(fits, s_cache, day = 84, threshold_pp = 0.5)
  expect_equal(el$hypothesis, c("H0", "H2", "H3"))
  expect_equal(el$eligible, c(FALSE, FALSE, TRUE))
  expect_true(all(el$threshold_pp == 0.5))
  expect_lt(abs(el$net_effect_pp[el$hypothesis == "H2"]), 0.5)
  expect_match(el$power_role[el$hypothesis == "H0"], "type I error")
  expect_match(el$power_role[el$hypothesis == "H3"], "power")
  expect_match(el$power_role[el$hypothesis == "H2"], "not run")
})

test_that("sobol_hypothesis gives one row per parameter and day, and a cluster map gives identical numbers", {
  f <- s_fit("H3", n = 5)
  a <- suppressWarnings(sobol_hypothesis(f, s_cache, days = c(14, 84), N = 60, nboot = 20, seed = 3))   # boot warns at nboot = 20
  expect_equal(sort(unique(a$day)), c(14, 84))
  expect_setequal(unique(a$parameter), c("log_kpbo", "log_kout", "logit_Emax", "log10_Kd_OGFr"))
  expect_equal(nrow(a), 2 * 4)
  expect_true(all(c("hypothesis", "day", "S1", "S1_lo", "S1_hi", "ST", "ST_lo", "ST_hi", "net_effect_median_pp") %in% names(a)))
  b <- with_cluster(2, function(cl) suppressWarnings(sobol_hypothesis(f, s_cache, days = c(14, 84), N = 60, nboot = 20, seed = 3,
                                                    map = function(X, FUN) parallel::parLapplyLB(cl, X, FUN))), root = root)
  expect_identical(a, b)
})

test_that("power_curves labels roles, matches sequential on a cluster, and H0's rate is the type-I check", {
  fits <- list(H0 = s_fit("H0"), H3 = s_fit("H3", emax = 0.4))
  oms <- list(H0 = list(omega_p = 0.3, omega_e = 0, sd_res = 0.5), H3 = list(omega_p = 0.3, omega_e = 0.1, sd_res = 0.5))
  designs <- list(T = list(trial = "T", schedule_ldn = "0:4.5", day = 28, n_ldn = 20, n_pbo = 20, base_mean = 6, base_sd = 1.2))
  a <- power_curves(fits, oms, designs, s_cache, hyps = "H3", ns = c(20, 60), n_rep = 30, n_draws = 3)
  expect_setequal(a$hypothesis, c("H0", "H3"))
  expect_equal(nrow(a), 2 * 2)
  expect_match(a$role[a$hypothesis == "H0"][1], "type I error")
  expect_match(a$role[a$hypothesis == "H3"][1], "power")
  expect_true(all(a$rej_rate_with_model_unc >= 0 & a$rej_rate_with_model_unc <= 1))
  h3 <- a[a$hypothesis == "H3", ]; expect_gte(h3$rej_rate_with_model_unc[2], h3$rej_rate_with_model_unc[1])
  b <- with_cluster(3, function(cl) power_curves(fits, oms, designs, s_cache, hyps = "H3", ns = c(20, 60), n_rep = 30,
                                                n_draws = 3, map = function(X, FUN) parallel::parLapplyLB(cl, X, FUN)), root = root)
  expect_identical(a, b)
})

test_that("thin_fit keeps exactly the draws virtual_trials would pick, so results are unchanged", {
  f <- s_fit("H3", n = 40)
  design <- list(trial = "T", schedule_ldn = "0:4.5", day = 14, n_ldn = 10, n_pbo = 10, base_mean = 6, base_sd = 1)
  om <- list(omega_p = 0.3, omega_e = 0.1, sd_res = 0.5)
  a <- virtual_trials(f, design, s_cache, om, n_rep = 5, n_draws = 7, seed = 4)
  b <- virtual_trials(thin_fit(f, 7), design, s_cache, om, n_rep = 5, n_draws = 7, seed = 4)
  expect_equal(a[, setdiff(names(a), "draw")], b[, setdiff(names(b), "draw")])
})

test_that("run_vt_task pins the generator kind, so the result does not depend on the caller's RNG kind", {
  f <- s_fit("H0")
  task <- list(fit = thin_fit(f, 2), design = list(trial = "T", schedule_ldn = "0:4.5", day = 14, n_ldn = 10, n_pbo = 10,
               base_mean = 6, base_sd = 1), cache = s_cache, om = list(omega_p = 0.3, omega_e = 0, sd_res = 0.5),
               n_rep = 4, n_draws = 2, seed = 7)
  set.seed(1, kind = "Mersenne-Twister"); a <- run_vt_task(task)
  set.seed(1, kind = "L'Ecuyer-CMRG"); b <- run_vt_task(task)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG"); set.seed(1, kind = "Mersenne-Twister")
  expect_identical(a, b)
  expect_identical(run_vt_task(modifyList(task, list(keep = "diff"))), a$diff)
})

# Jansen's FIRST-order estimator is badly biased upward on the heavy-tailed net-effect outputs (an inert log_kpbo got
# S1 = 0.17 with ST = 2e-5 on the H1 posterior); S1 comes from the Sobol 2007 estimator and ST from Jansen, both on
# the SAME evaluations (the two designs are identical).
test_that("sobol_hypothesis takes S1 from sobol2007 and ST from Jansen, so an inert parameter has S1 near 0 and S1 <= ST", {
  f <- s_fit("H3", n = 9)
  f$draws[, "log_kpbo"] <- seq(log(0.05), log(0.2), length.out = 9)     # inert at day 84: the placebo term has saturated
  f$draws[, "logit_Emax"] <- seq(logit(0.05), logit(0.9), length.out = 9)   # wide, skewed output
  a <- suppressWarnings(sobol_hypothesis(f, s_cache, days = 84, N = 1000, nboot = 20, seed = 5))
  expect_true(all(c("S1_estimator", "ST_estimator") %in% names(a)))
  expect_match(a$S1_estimator[1], "sobol2007"); expect_match(a$ST_estimator[1], "Jansen")
  inert <- a[a$parameter %in% c("log_kpbo", "log_kout"), ]       # both saturated by day 84 in this fixture
  expect_true(all(abs(inert$S1) < 0.005))                         # Jansen's S1 gives about -0.018 for both here
  expect_true(all(abs(inert$S1_jansen) > abs(inert$S1)))
  expect_false("S1_jansen_biased" %in% names(a))
  expect_true(all(a$S1 <= a$ST + 0.03))
})

test_that("sobol_skipped rows bind with sobol_hypothesis rows (a below-threshold hypothesis gets a row, not silence)", {
  a <- suppressWarnings(sobol_hypothesis(s_fit("H3", n = 5), s_cache, days = c(14, 84), N = 40, nboot = 10, seed = 3))
  s <- sobol_skipped("H2", c(14, 84), c(0.0009, 0.0012))
  expect_identical(names(s), names(a))
  expect_equal(nrow(s), 2); expect_true(all(is.na(s$S1) & is.na(s$parameter)))
  expect_equal(s$net_effect_median_pp, c(0.0009, 0.0012))
  expect_silent(rbind(a, s))
})


# ---- R3a fix round 1 ----
test_that("power roles: eligible hypotheses are labelled assurance (posterior-averaged power)", {
  fits <- list(H0 = s_fit("H0"), H3 = s_fit("H3"))
  el <- effect_eligibility(fits, s_cache, day = 84, threshold_pp = 0.5)
  expect_identical(el$power_role[el$hypothesis == "H3"], "assurance (posterior-averaged power)")
  expect_identical(el$power_role[el$hypothesis == "H0"], "type I error check (no drug term)")
})

test_that("power_curves writes an explanatory row for a hypothesis that is in roles but not run", {
  fits <- list(H0 = s_fit("H0"), H3 = s_fit("H3", emax = 0.4))
  oms <- list(H0 = list(omega_p = 0.3, omega_e = 0, sd_res = 0.5), H3 = list(omega_p = 0.3, omega_e = 0.1, sd_res = 0.5))
  designs <- list(T = list(trial = "T", schedule_ldn = "0:4.5", day = 28, n_ldn = 20, n_pbo = 20, base_mean = 6, base_sd = 1.2))
  roles <- c(H0 = "type I error check (no drug term)", H2 = "not run: net effect 0.00122 pp is below the 0.5 pp threshold",
             H3 = "assurance (posterior-averaged power)")
  a <- power_curves(fits, oms, designs, s_cache, hyps = "H3", ns = 20, n_rep = 10, n_draws = 2, roles = roles)
  h2 <- a[a$hypothesis == "H2", ]
  expect_equal(nrow(h2), 1)
  expect_true(is.na(h2$rej_rate_with_model_unc) && is.na(h2$n_per_arm))
  expect_match(h2$role, "not run"); expect_match(h2$role, "effect_eligibility.csv", fixed = TRUE)
  expect_equal(a$role[a$hypothesis == "H3"], "assurance (posterior-averaged power)")
})

test_that("Sobol runs at N = 20000 by default and script 09 runs two seeds with a seed-agreement file", {
  expect_equal(formals(sobol_hypothesis)$N, 20000)
  s9 <- paste(readLines(file.path("..", "..", "scripts", "09_sobol.R")), collapse = "\n")
  expect_true(grepl("SOBOL_SEEDS <- c(3, 17)", s9, fixed = TRUE))
  expect_true(grepl("sobol_seed_agreement.csv", s9, fixed = TRUE))
  expect_false(grepl("N <- 2000;", s9, fixed = TRUE))
})

test_that("sobol_seed_agreement flags an index whose two seeds differ by more than the larger CI half-width", {
  mk <- function(s1, st, hw) data.frame(hypothesis = "H1", day = 14, parameter = c("a", "b"), S1 = s1, S1_lo = s1 - hw,
                                        S1_hi = s1 + hw, ST = st, ST_lo = st - hw, ST_hi = st + hw)
  a <- mk(c(0.10, 0.30), c(0.20, 0.50), 0.02); b <- mk(c(0.11, 0.20), c(0.21, 0.52), 0.03)
  g <- sobol_seed_agreement(a, b, seeds = c(3, 17))
  expect_equal(nrow(g), 4)                                       # 2 parameters x (S1, ST)
  expect_true(all(c("hypothesis", "day", "parameter", "index", "estimate_seed3", "estimate_seed17", "abs_diff",
                    "max_ci_half_width", "disagree") %in% names(g)))
  expect_equal(g$disagree[g$parameter == "a" & g$index == "S1"], FALSE)   # 0.01 <= 0.03
  expect_equal(g$disagree[g$parameter == "b" & g$index == "S1"], TRUE)    # 0.10 > 0.03
  expect_equal(g$disagree[g$parameter == "b" & g$index == "ST"], FALSE)   # 0.02 <= 0.03
  expect_error(sobol_seed_agreement(a, b[1, ], seeds = c(3, 17)), "same")
})
