lpd_normal <- function(obs, pred_draws, se) {
  d <- stats::dnorm(obs, pred_draws, se, log = TRUE); m <- max(d)
  m + log(mean(exp(d - m)))
}

# rows: arms of ONE new trial (arm-level endpoints change_nrs); a new trial-level placebo draw is shared by its arms.
# Each arm prediction also carries a model-structure residual N(0, sigma_rel * baseline_nrs) (0 if the fit has no log_sigma_rel).
predict_trial <- function(fit, rows, cache, n_draws = 1000, seed = 7) {
  set.seed(seed)
  idx <- round(seq(1, nrow(fit$draws), length.out = n_draws))
  out <- lapply(seq_along(idx), function(k) {
    par <- fit$draws[idx[k], ]
    th <- theta_from(par, fit$hyp, fit$fixed)
    sigma_rel <- if ("log_sigma_rel" %in% names(par)) exp(par[["log_sigma_rel"]]) else 0
    pmax_new <- inv_logit(par[["mu_p"]] + stats::rnorm(1, 0, exp(par[["log_tau"]])))
    vapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      sim <- simulate_arm(fit$hyp, th, parse_schedule(r$schedule), r$day, cache, pmax_new)
      z <- stats::rnorm(1)  # always drawn so the RNG stream is identical whatever sigma is
      arm_outcome(sim, r$day, "change_nrs", r$baseline_nrs) + z * sigma_rel * outcome_scale("change_nrs", r$baseline_nrs)
    }, 0)
  })
  m <- do.call(rbind, out)
  data.frame(row = rep(seq_len(nrow(rows)), each = length(idx)), draw = rep(seq_along(idx), nrow(rows)), pred = as.vector(m))
}
