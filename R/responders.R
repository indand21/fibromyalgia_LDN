# Posterior predictive p-values, responder scoring and the SD sensitivity input (R3a items 3 and 4, review M3, I6, I10).

# Posterior predictive p-values for an observed value against simulated values (review M3). The old p_tail was
# one-sided with a direction that flipped with the sign of the observation, so FINAL (obs < 0) and INNOVA
# (obs > 0) values were not comparable. Both tails are returned, the old quantity under an explicit name, and
# the two-sided value 2 x the smaller tail (capped at 1).
ppp_values <- function(sim, obs) {
  if (length(obs) != 1L || !is.finite(obs))
    return(list(p_upper = NA_real_, p_lower = NA_real_, p_one_sided_toward_sign_of_obs = NA_real_, p_two_sided = NA_real_))
  up <- mean(sim >= obs); lo <- mean(sim <= obs)
  list(p_upper = up, p_lower = lo, p_one_sided_toward_sign_of_obs = if (obs >= 0) up else lo,
       p_two_sided = min(1, 2 * min(up, lo)))
}

# Observed responder COUNTS per arm and cut for one trial and day. The curated rows hold proportions in `mean` and
# the counts in `location` as "(k of n)"; the count is parsed and cross-checked against mean x n and the row's n.
observed_responder_counts <- function(tt, trial, day, cuts = c(15, 30, 50)) {
  one <- function(cut, pbo) {
    r <- tt[tolower(tt$trial) == tolower(trial) & tt$day == day & tt$endpoint == paste0("resp", cut) &
              ((tt$arm == "placebo") == pbo), ]
    if (nrow(r) != 1) stop("observed_responder_counts: need exactly one resp", cut, " ", if (pbo) "placebo" else "LDN",
                           " row for ", trial, " day ", day, " (found ", nrow(r), ")")
    m <- regmatches(r$location, regexec("\\((\\d+) of (\\d+)\\)", r$location))[[1]]
    if (length(m) != 3) stop("observed_responder_counts: no '(k of n)' count in the location of ", trial, " resp", cut)
    k <- as.integer(m[2]); n <- as.integer(m[3])
    p <- if (r$mean > 1) r$mean / 100 else r$mean
    if (n != r$n || abs(k / n - p) > 0.01)
      stop(sprintf("observed_responder_counts: %s resp%d %s: count %d of %d and proportion %g (n %g) disagree",
                   trial, cut, if (pbo) "placebo" else "LDN", k, n, r$mean, r$n))
    c(k = k, n = n)
  }
  do.call(rbind, lapply(cuts, function(cut) { a <- one(cut, FALSE); b <- one(cut, TRUE)
    data.frame(trial = trial, day = day, cut = cut, k_ldn = a[["k"]], n_ldn = a[["n"]], k_pbo = b[["k"]], n_pbo = b[["n"]]) }))
}

# Score observed responder counts against model-predicted responder RATES (review I6).
# rate_ldn / rate_pbo: paired vectors of population-level responder probabilities, one pair per posterior draw and
# replicate trial (the two arms of one replicate share its placebo ceiling and model-structure residual).
# Each arm's count is sampled Binomial(n_arm, rate) at the trial's ACTUAL arm sizes, n_binom times per pair.
#   pp_prob_count_*      P(count = observed) per arm, and jointly (arms conditionally independent given the
#                        rates); the exact binomial expectation (mean of dbinom), i.e. the Rao-Blackwellised
#                        version of the sampled indicator
#   ppp_two_sided_count_* 2 x the smaller tail of the sampled count at the observed count, per arm
#   ppp_one_sided_upper  P(simulated LDN-minus-placebo difference in proportions >= observed): the direction of an
#                        LDN benefit, fixed for every trial so FINAL and INNOVA are comparable (M3)
#   ppp_two_sided        2 x the smaller tail of that difference (capped at 1)
# The source trial's own Fisher exact test on the observed 2 x 2 table is attached as descriptive context only.
responder_scores <- function(rate_ldn, rate_pbo, k_ldn, n_ldn, k_pbo, n_pbo, seed, n_binom = 20) {
  if (length(rate_ldn) != length(rate_pbo)) stop("responder_scores: rate_ldn and rate_pbo must have the same length (paired draws)")
  S <- length(rate_ldn)
  K <- with_mt_seed(seed, list(l = stats::rbinom(S * n_binom, n_ldn, rep(rate_ldn, each = n_binom)),
                               p = stats::rbinom(S * n_binom, n_pbo, rep(rate_pbo, each = n_binom))))
  dl <- stats::dbinom(k_ldn, n_ldn, rate_ldn); dp <- stats::dbinom(k_pbo, n_pbo, rate_pbo)
  two <- function(x, o) min(1, 2 * min(mean(x >= o), mean(x <= o)))
  d_sim <- K$l / n_ldn - K$p / n_pbo; d_obs <- k_ldn / n_ldn - k_pbo / n_pbo
  pv <- ppp_values(d_sim, d_obs)
  ft <- stats::fisher.test(matrix(c(k_ldn, n_ldn - k_ldn, k_pbo, n_pbo - k_pbo), 2, byrow = TRUE))
  data.frame(k_ldn = k_ldn, n_ldn = n_ldn, k_pbo = k_pbo, n_pbo = n_pbo,
             obs_rate_ldn = k_ldn / n_ldn, obs_rate_pbo = k_pbo / n_pbo, obs_diff = d_obs,
             pred_rate_ldn_median = stats::median(rate_ldn), pred_rate_pbo_median = stats::median(rate_pbo),
             pred_diff_median = stats::median(d_sim), pred_diff_lo90 = unname(stats::quantile(d_sim, 0.05)),
             pred_diff_hi90 = unname(stats::quantile(d_sim, 0.95)),
             pp_prob_count_ldn = mean(dl), pp_prob_count_pbo = mean(dp), pp_prob_counts_joint = mean(dl * dp),
             ppp_two_sided_count_ldn = two(K$l, k_ldn), ppp_two_sided_count_pbo = two(K$p, k_pbo),
             ppp_one_sided_upper = pv$p_upper, ppp_one_sided_lower = pv$p_lower, ppp_two_sided = pv$p_two_sided,
             n_rate_pairs = S, n_pred_counts = S * n_binom,
             fisher_p_two_sided = ft$p.value, fisher_odds_ratio = unname(ft$estimate),
             fisher_label = "source-trial test, descriptive")
}

# SD sensitivity input (review I10): Bested et al. 2023 quotes Younger 2013's statistics with SDs of 14.6 (placebo)
# and 12.5 (LDN), against the 27.9 / 24.0 derived from the reported 95% CI. Returns a COPY of the training table
# with those sd_pct values on the Younger 2013 rows; the curated file is never changed.
override_sd_pct <- function(arms, pbo, ldn, trial_pattern = "Younger.*2013") {
  y <- grepl(trial_pattern, arms$trial, ignore.case = TRUE)
  if (!any(y & arms$arm == "placebo") || !any(y & arms$arm != "placebo"))
    stop("override_sd_pct: no Younger 2013 placebo and LDN rows to override")
  arms$sd_pct[y & arms$arm == "placebo"] <- pbo
  arms$sd_pct[y & arms$arm != "placebo"] <- ldn
  arms
}
