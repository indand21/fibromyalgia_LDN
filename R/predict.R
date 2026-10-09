lpd_normal <- function(obs, pred_draws, se) {
  d <- stats::dnorm(obs, pred_draws, se, log = TRUE); m <- max(d)
  m + log(mean(exp(d - m)))
}

# Monte Carlo standard error of a held-out elpd (review I8). Each lpd is log(mean_s w_is) over a FINITE set of
# S posterior-predictive draws, w_is = dnorm(obs_i, pred_is, se_i), so elpd = sum_i log(m_i), m_i = mean_s w_is,
# carries draw error of its own. The delta method gives
#   Var(elpd) ~= sum_ij Cov(m_i, m_j) / (m_i m_j) = (1/S) Var_s( sum_i w_is / m_i ),
# i.e. MCSE = sd_s(sum_i w_is/m_i) / sqrt(S), which for one observation is sd(w)/(mean(w) sqrt(S)) and which
# correctly adds the covariance when observations are scored on the SAME draws.
# This is predictive-DRAW error only. It is NOT the observation-level sampling SE of an elpd contrast (which with
# n = 2 held-out endpoints is the quantity already reported as se_diff), and it does not cover a refit with other
# MCMC chains or other inputs.
# `pred` is an S x n matrix (or a list of n equal-length draw vectors).
elpd_mcse <- function(obs, pred, se) {
  q <- elpd_ratio_draws(obs, pred, se)
  list(elpd = attr(q, "elpd"), mcse = stats::sd(q) / sqrt(length(q)), n_draws = length(q), n_obs = attr(q, "n_obs"))
}

# Per-draw influence q_s = sum_i w_is / m_i of an elpd (the quantity whose SD/sqrt(S) is the delta-method MCSE),
# with the elpd itself attached as an attribute.
elpd_ratio_draws <- function(obs, pred, se) {
  pm <- if (is.matrix(pred)) pred else if (is.list(pred)) do.call(cbind, pred) else matrix(pred, ncol = 1)
  n <- ncol(pm); S <- nrow(pm)
  if (length(obs) != n) stop("elpd_mcse: obs has ", length(obs), " entries for ", n, " prediction column(s)")
  if (length(se) == 1L) se <- rep(se, n)
  if (length(se) != n) stop("elpd_mcse: se has ", length(se), " entries for ", n, " prediction column(s)")
  if (S < 2L) stop("elpd_mcse: need at least 2 draws to estimate a Monte Carlo standard error")
  ratio <- vapply(seq_len(n), function(i) {
    d <- stats::dnorm(obs[i], pm[, i], se[i], log = TRUE)
    e <- exp(d - max(d))               # the shift cancels in w / mean(w)
    e / mean(e)
  }, numeric(S))
  q <- rowSums(matrix(ratio, nrow = S))
  attr(q, "elpd") <- sum(vapply(seq_len(n), function(i) lpd_normal(obs[i], pm[, i], se[i]), 0))
  attr(q, "n_obs") <- n
  q
}

# Monte Carlo SE of an elpd CONTRAST (R2 fix round 1, Minor 3). Hypotheses are predicted with the same per-trial
# predictive seeds, so their draw errors are positively correlated and sqrt(mcse_a^2 + mcse_b^2) overstates the
# error of a difference. By the delta method on PAIRED draws, Var(elpd_a - elpd_b) ~= Var_s(q_as - q_bs) / S,
# so the contrast MCSE is sd_s(q_a - q_b) / sqrt(S). The two prediction matrices must be paired draw by draw.
elpd_mcse_diff <- function(obs, pred_a, pred_b, se) {
  qa <- elpd_ratio_draws(obs, pred_a, se); qb <- elpd_ratio_draws(obs, pred_b, se)
  if (length(qa) != length(qb)) stop("elpd_mcse_diff: the two prediction sets need the same number of draws")
  d <- as.vector(qa) - as.vector(qb)
  list(diff = attr(qa, "elpd") - attr(qb, "elpd"), mcse = stats::sd(d) / sqrt(length(d)), n_draws = length(d))
}

# rows: arms of ONE new trial (arm-level endpoints change_nrs); a new trial-level placebo draw is shared by its arms.
# The model-structure residual is also a TRIAL-level term (review I7): ONE z ~ N(0, 1) per draw multiplies
# sigma_rel in every arm of the trial, so a misspecified disease course moves both arms together. Each arm's
# absolute level still carries the full residual N(0, sigma_rel * baseline_nrs), while it cancels exactly in the
# difference between two arms of the same trial. Independent trials (separate calls) get independent terms.
# sigma_rel is 0 if the fit has no log_sigma_rel.
predict_trial <- function(fit, rows, cache, n_draws = 1000, seed = 7) {
  set.seed(seed)
  idx <- round(seq(1, nrow(fit$draws), length.out = n_draws))
  out <- lapply(seq_along(idx), function(k) {
    par <- fit$draws[idx[k], ]
    th <- theta_from(par, fit$hyp, fit$fixed)
    sigma_rel <- if ("log_sigma_rel" %in% names(par)) exp(par[["log_sigma_rel"]]) else 0
    pmax_new <- inv_logit(par[["mu_p"]] + stats::rnorm(1, 0, exp(par[["log_tau"]])))
    z <- stats::rnorm(1)    # ONE shared model-structure residual for this trial-draw; always drawn so the
                            # RNG stream is identical whatever sigma_rel is
    vapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      sim <- simulate_arm(fit$hyp, th, parse_schedule(r$schedule), r$day, cache, pmax_new)
      arm_outcome(sim, r$day, "change_nrs", r$baseline_nrs) + z * sigma_rel * outcome_shift("change_nrs", r$baseline_nrs)
    }, 0)
  })
  m <- do.call(rbind, out)
  data.frame(row = rep(seq_len(nrow(rows)), each = length(idx)), draw = rep(seq_along(idx), nrow(rows)), pred = as.vector(m))
}

# Held-out intervals (R3a item 1, controller ruling after R2). Since R2 shares the model-structure residual within a
# trial, the draws of a predicted DIFFERENCE describe the model's uncertainty about the TRUE difference; they carry
# none of the trial's own sampling error. An observed value must be compared with the posterior predictive
# interval for the OBSERVED value, formed by adding N(0, obs_se) observation error to every draw:
#   true_lo90/true_hi90  5-95% quantiles of the model draws (TRUE value)
#   pred_lo90/pred_hi90  5-95% quantiles of draws + N(0, obs_se) (OBSERVED value); NA when no obs_se exists
#   inside90             observed value inside the PREDICTIVE interval (NA when there is none)
# The observation-error draws use their own Mersenne-Twister stream from `seed` (with_mt_seed), so the result is
# reproducible and the caller's RNG stream is untouched. The lpd scoring (lpd_normal) already uses obs_se.
HELDOUT_INTERVAL_DEF <- paste(
  "pred_lo90/pred_hi90 = 90% posterior predictive interval for the OBSERVED value (model draws plus N(0, obs_se)",
  "observation error per draw); true_lo90/true_hi90 = 90% interval for the TRUE value (model uncertainty only);",
  "inside90 compares the observed value with the predictive interval")
HELDOUT_INTERVAL_DEF_NOSE <- paste(
  "no observation SE for this row (not reported and not derivable), so no predictive interval was formed:",
  "pred_lo90/pred_hi90 and inside90 are NA; true_lo90/true_hi90 = 90% interval for the TRUE value (model uncertainty only)")
heldout_intervals <- function(obs, draws, obs_se, seed) {
  tq <- unname(stats::quantile(draws, c(0.05, 0.95)))
  has_se <- length(obs_se) == 1L && is.finite(obs_se) && obs_se > 0
  if (has_se) {
    noisy <- with_mt_seed(seed, draws + stats::rnorm(length(draws), 0, obs_se))
    pq <- unname(stats::quantile(noisy, c(0.05, 0.95)))
    inside <- if (is.finite(obs)) obs >= pq[1] & obs <= pq[2] else NA
  } else { pq <- c(NA_real_, NA_real_); inside <- NA }
  data.frame(true_lo90 = tq[1], true_hi90 = tq[2], pred_lo90 = pq[1], pred_hi90 = pq[2], inside90 = inside,
             interval_definition = if (has_se) HELDOUT_INTERVAL_DEF else HELDOUT_INTERVAL_DEF_NOSE)
}
